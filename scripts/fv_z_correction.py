"""
FV correction + z-correction for S1-only alpha analysis.

Requires peak-level position data.  Two sources (in preference order):
  1. peak_s1_positions_cnn  — built-in straxen, S1 CNN (x,y,z in cm)
  2. naive_peak_positions   — custom plugin, PMT-area-weighted

The full chain:
  1. Load peak_basics + position data → compute r²
  2. FV correction:  area_corrected = area × ref / linearFunc(r²)
  3. z-correction:   area_corrected = area × ref / exp_bump_model(AFT)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
import pandas as pd

from naive_peak_positions import (
    merge_position_into_basics,
    load_positions_for_df,
    register_plugin,
)


# ═══════════════════════════════════════════════════════════════════════
# Notebook integration — data loading
# ═══════════════════════════════════════════════════════════════════════

# After creating the straxen context:
#
#   st = straxen.contexts.xenonnt_online(...)
#   st.storage.append(...)
#
# Modify the loading loop to also load peak_s1_positions_cnn:
#
#   for chunk in st.get_iter(
#       this_run_name,
#       ["peak_basics", "peak_s1_positions_cnn"],   # ← add position
#       selection_str=cut_str,
#       progress_bar=False
#   ):
#       if len(chunk.data) > 0:
#           merged = merge_position_into_basics(      # ← merge + compute r2
#               chunk.data, columns_to_keep, pos_source='cnn'
#           )
#           run_chunks.append(merged)
#
# This gives you columns: area, area_fraction_top, ..., x, y, r2.


# ═══════════════════════════════════════════════════════════════════════
# FV correction — fit & correct area vs r²
# ═══════════════════════════════════════════════════════════════════════

def linear_func(x, a, b):
    """Linear function: a*x + b.  x = r² in cm²."""
    return a * x + b


def plot_r2_dependence(df, r2_col='r2', aft_range=(0.12, 0.55)):
    """
    Diagnostic: show the radial dependence that FV correction removes.
    Area vs r², colored by AFT.
    """
    mask = (df['area_fraction_top'] >= aft_range[0]) & \
           (df['area_fraction_top'] <= aft_range[1]) & \
           (df[r2_col] > 0) & (df[r2_col] < 6e3)
    d = df.loc[mask]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    h = ax.hist2d(d[r2_col], d['area'], bins=100,
                  range=[[0, 6e3], [20e3, 100e3]],
                  cmin=1, norm=matplotlib.colors.LogNorm())
    plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('r² [cm²]')
    ax.set_ylabel('Area [PE]')
    ax.set_title('Area vs r² (before FV correction)')

    # Median in r² slices
    ax = axes[1]
    r2_edges = np.linspace(0, 6e3, 31)
    d['r2_bin'] = pd.cut(d[r2_col], bins=r2_edges)
    med = d.groupby('r2_bin', observed=True)['area'].median()
    cents = med.index.map(lambda iv: iv.mid).astype(float).values
    ax.plot(cents, med.values, 'o-', color='red', lw=2, label='Median area')
    ax.set_xlabel('r² [cm²]')
    ax.set_ylabel('Median Area [PE]')
    ax.set_title('Median area vs r² (AFT slice)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    return fig


def apply_fv_correction(df, r2_col='r2', aft_range=(0.12, 0.55),
                        r2_range=(0, 5e3), n_r2_bins=30,
                        plot_fit=True):
    """
    Fit area vs r² in a narrow AFT slice, then divide out the r² dependence.

    The fit is done on a narrow AFT band so the z-dependence doesn't
    contaminate the r²-dependence measurement.

    Returns
    -------
    df : pd.DataFrame
        With new column 'area_fv_corrected'.
    popt : tuple
        (slope, intercept) of linear fit.
    """
    df = df.copy()

    mask = (
        (df['area_fraction_top'] >= aft_range[0]) &
        (df['area_fraction_top'] <= aft_range[1]) &
        (df[r2_col] >= r2_range[0]) &
        (df[r2_col] <= r2_range[1])
    )
    df_slice = df.loc[mask].copy()

    r2_edges = np.linspace(r2_range[0], r2_range[1], n_r2_bins + 1)
    df_slice['r2_bin'] = pd.cut(df_slice[r2_col], bins=r2_edges)

    medians = df_slice.groupby('r2_bin', observed=True)['area'].median()
    r2_centers = medians.index.map(lambda iv: iv.mid).astype(float).values
    median_areas = medians.values

    counts = df_slice.groupby('r2_bin', observed=True)['area'].count()
    good = counts.values >= 10
    r2_centers = r2_centers[good]
    median_areas = median_areas[good]

    if len(r2_centers) < 5:
        print("FV correction: too few bins, skipping")
        df['area_fv_corrected'] = df['area']
        return df, (0, np.median(df['area']))

    p0 = [0.005, np.median(median_areas)]
    try:
        popt, pcov = curve_fit(linear_func, r2_centers, median_areas,
                               p0=p0, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        print(f"FV fit:  slope = {popt[0]:.6f} ± {perr[0]:.6f}  PE/cm²")
        print(f"         intercept = {popt[1]:.0f} ± {perr[1]:.0f}  PE")
    except Exception as e:
        print(f"FV fit failed: {e}, using constant")
        popt = [0, np.median(median_areas)]
        perr = [0, 0]

    area_ref_fv = linear_func(0, *popt)
    r2_values = df[r2_col].values
    correction_factor = area_ref_fv / linear_func(r2_values, *popt)
    correction_factor = np.clip(correction_factor, 0.3, 3.0)

    df['area_fv_corrected'] = df['area'] * correction_factor
    df['fv_correction_factor'] = correction_factor

    print(f"Ref area (FV, r²=0) = {area_ref_fv:.0f} PE")
    print(f"Correction range: [{correction_factor.min():.3f}, {correction_factor.max():.3f}]")

    if plot_fit:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.errorbar(r2_centers, median_areas,
                    yerr=np.sqrt(median_areas) / np.sqrt(counts[good].values),
                    fmt='o', capsize=3, label='Median area')
        r2_fine = np.linspace(0, r2_range[1], 200)
        ax.plot(r2_fine, linear_func(r2_fine, *popt), 'r-', lw=2, label='Fit')
        ax.axhline(area_ref_fv, color='gray', ls='--', label=f'Ref = {area_ref_fv:.0f} PE')
        ax.set_xlabel('r² [cm²]')
        ax.set_ylabel('Area [PE]')
        ax.set_title('FV correction fit')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return df, popt


# ═══════════════════════════════════════════════════════════════════════
# z-correction (AFT correction)
# ═══════════════════════════════════════════════════════════════════════

def exp_bump_model(x, a, b, c, d, e):
    """
    Florian's empirical model for the alpha band in Area-vs-AFT space.

        f(AFT) = a × exp((AFT + b)/c) × (1 − exp((AFT + b)/d)) + e
    """
    return a * np.exp((x + b) / c) * (1 - np.exp((x + b) / d)) + e


def apply_aft_correction(df, area_col='area_fv_corrected',
                         aft_range=(0.12, 0.55), aft_ref=0.30,
                         n_aft_bins=50,
                         band_upper_pars=None, band_lower_pars=None,
                         plot_fit=True):
    """
    Fit the alpha band centerline in area-vs-AFT, then divide out the
    AFT-dependent light collection variation.

    Apply AFTER FV correction.

    Returns
    -------
    df : pd.DataFrame
        With new column 'area_z_corrected'.
    popt : tuple
        Fitted (a, b, c, d, e).
    """
    df = df.copy()

    mask_fv = (
        (df['area_fraction_top'] >= aft_range[0]) &
        (df['area_fraction_top'] <= aft_range[1])
    )
    df_fv = df.loc[mask_fv].copy()

    if band_upper_pars is not None and band_lower_pars is not None:
        df_fv['in_rn_band'] = (
            (df_fv[area_col] <
             exp_bump_model(df_fv['area_fraction_top'], *band_upper_pars)) &
            (df_fv[area_col] >
             exp_bump_model(df_fv['area_fraction_top'], *band_lower_pars))
        )
    else:
        df_fv['in_rn_band'] = (
            (df_fv[area_col] > 30e3) & (df_fv[area_col] < 55e3)
        )

    n_rn = df_fv['in_rn_band'].sum()
    print(f"Events in Rn band: {n_rn}")
    if n_rn < 100:
        print("WARNING: very few Rn band events — fit unreliable")

    aft_edges = np.linspace(aft_range[0], aft_range[1], n_aft_bins + 1)
    df_fv['aft_bin'] = pd.cut(df_fv['area_fraction_top'], bins=aft_edges)

    medians_rn = (df_fv[df_fv['in_rn_band']]
                  .groupby('aft_bin', observed=True)[area_col]
                  .median())
    aft_centers = medians_rn.index.map(lambda iv: iv.mid).astype(float).values
    median_values = medians_rn.values

    good = ~np.isnan(median_values)
    aft_centers_fit = aft_centers[good]
    median_values_fit = median_values[good]

    print(f"Fitting AFT with {len(aft_centers_fit)} bins")
    if len(aft_centers_fit) < 8:
        raise ValueError("Too few AFT bins")

    p0 = [4.5e4, -0.05, -0.1, -0.05, 35e3]
    popt, pcov = curve_fit(
        exp_bump_model, aft_centers_fit, median_values_fit,
        p0=p0,
        bounds=([0, -np.inf, -np.inf, -np.inf, 0],
                [np.inf, 0, 0, 0, np.inf]),
        maxfev=20000
    )
    perr = np.sqrt(np.diag(pcov))
    names = ['a(ampl)', 'b(shift)', 'c(scale₁)', 'd(scale₂)', 'e(offset)']
    print("AFT fit:")
    for n, v, e in zip(names, popt, perr):
        print(f"  {n:12s} = {v:.6e} ± {e:.6e}")

    area_ref_aft = exp_bump_model(aft_ref, *popt)
    print(f"Ref area (AFT={aft_ref}) = {area_ref_aft:.0f} PE")

    aft_values = df['area_fraction_top'].values
    aft_factor = area_ref_aft / exp_bump_model(aft_values, *popt)
    aft_factor = np.clip(aft_factor, 0.3, 3.0)

    df['area_z_corrected'] = df[area_col] * aft_factor
    df['aft_correction_factor'] = aft_factor

    print(f"AFT correction range: [{aft_factor.min():.3f}, {aft_factor.max():.3f}]")

    if plot_fit:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(median_values_fit, aft_centers_fit, 'o', color='lime',
                markersize=5, label='Medians')
        aft_fine = np.linspace(aft_range[0], aft_range[1], 200)
        ax.plot(exp_bump_model(aft_fine, *popt), aft_fine, 'r-', lw=2,
                label='Fit')
        ax.axvline(area_ref_aft, color='gray', ls='--',
                   label=f'Ref = {area_ref_aft:.0f} PE')
        ax.set_xlabel('Area [PE]')
        ax.set_ylabel('Area Fraction Top')
        ax.set_title('AFT correction fit')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return df, popt


# ═══════════════════════════════════════════════════════════════════════
# Comparison plots
# ═══════════════════════════════════════════════════════════════════════

def plot_correction_comparison(df, area_cols=None):
    """Side-by-side Area vs AFT at each correction stage."""
    if area_cols is None:
        area_cols = {}
        for label, col in [('Raw', 'area'),
                           ('FV-corrected', 'area_fv_corrected'),
                           ('FV+Z-corrected', 'area_z_corrected')]:
            if col in df.columns:
                area_cols[label] = col

    n = len(area_cols)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5.5), squeeze=False)

    if 'area' in df.columns:
        h_test, _, _ = np.histogram2d(
            df['area'], df['area_fraction_top'],
            bins=[200, 200], range=[[20e3, 100e3], [0.1, 0.65]]
        )
        vmax = h_test.max()
    else:
        vmax = None

    for ax, (label, col) in zip(axes[0], area_cols.items()):
        h = ax.hist2d(
            df[col], df['area_fraction_top'],
            bins=[200, 200],
            range=[[20e3, 100e3], [0.1, 0.65]],
            cmin=1,
            norm=matplotlib.colors.LogNorm(vmax=vmax)
        )
        plt.colorbar(h[3], ax=ax, label='Counts')
        ax.set_xlabel('Area [PE]')
        ax.set_ylabel('Area Fraction Top')
        ax.set_title(label)
        ax.grid(True, alpha=0.2)

    fig.suptitle('Area vs AFT — correction stages', fontsize=14, y=1.02)
    plt.tight_layout()
    return fig


def quantify_band_width(df, area_cols=None, aft_range=(0.15, 0.46),
                        area_range=(35e3, 50e3)):
    """
    Compute band width metrics at each correction stage.
    Narrower = better correction.
    """
    if area_cols is None:
        area_cols = {
            'Raw': 'area',
            'FV-corrected': 'area_fv_corrected',
            'FV+Z-corrected': 'area_z_corrected',
        }

    mask = (
        (df['area_fraction_top'] >= aft_range[0]) &
        (df['area_fraction_top'] <= aft_range[1])
    )

    results = {}
    print(f"{'Stage':20s}  {'μ':>8s}  {'σ':>8s}  "
          f"{'IQR':>8s}  {'FWHM_est':>10s}  {'n':>8s}")
    print("-" * 70)

    for label, col in area_cols.items():
        if col not in df.columns:
            continue
        data = df.loc[mask, col].dropna()
        data = data[(data >= area_range[0]) & (data <= area_range[1])]

        r = {
            'mean': data.mean(),
            'std': data.std(),
            'median': data.median(),
            'iqr': data.quantile(0.75) - data.quantile(0.25),
            'fwhm_est': 2.355 * data.std(),
            'n': len(data),
        }
        results[label] = r

        print(f"{label:20s}  {r['mean']:8.0f}  {r['std']:8.0f}  "
              f"{r['iqr']:8.0f}  {r['fwhm_est']:10.0f}  {r['n']:8d}")

    return results
