import json

with open(r'C:\Users\13696\Desktop\xenon\Rn222\rn222_s1s2.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

# New markdown cell — inserted after cell 15 (df display), before cell 16 (Alpha Event Selection Cuts)
new_md_cell = {
    "cell_type": "markdown",
    "metadata": {},
    "source": [
        "## Raw Data Visualization — Before Cuts\n",
        "\n",
        "Visualize the raw `event_info` data (pre-selected with `s1_area > 15000` PE)\n",
        "to **manually tune each alpha selection cut** for the SR3 dataset.\n",
        "\n",
        "Each panel corresponds to one cut stage in `apply_alpha_cuts()`:\n",
        "1. S2 vs S1 — gamma/electron removal threshold\n",
        "2. S1 Area vs AFT — alpha band identification (key plot)\n",
        "3. S1 Area vs Drift Time — DT-based z proxy\n",
        "4. Drift Time vs AFT — S1-S2 pair envelope\n",
        "5. Drift Time vs S2 Width — unphysical event removal\n",
        "6. S1 Area vs r² — radial uniformity / FV r² cut\n",
        "7. z Distribution — FV z cut\n",
        "8. AFT Distribution — BiPo AFT lower bound\n",
        "9. alt_s1_area vs s1_area — multi-scatter rejection\n",
        "\n",
        "**Use these plots to adjust the cut parameters in the next cell.**"
    ]
}

# New code cell
plot_code = r'''# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  RAW DATA VISUALIZATION — 3×3 diagnostic grid for cut tuning               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# Use these plots to manually adjust cut parameters in apply_alpha_cuts() below.
# Current cut values are overlaid as dashed red lines where applicable.

if df is not None:
    fig, axes = plt.subplots(3, 3, figsize=(22, 18))
    CMAP = 'plasma'
    SCATTER_ALPHA = 0.3
    SCATTER_S = 1

    # ── Panel 1: S2 vs S1 Area (gamma/electron removal) ──
    ax = axes[0, 0]
    if 's2_area' in df.columns:
        h = ax.hist2d(df['s1_area'], df['s2_area'],
                      bins=[200, 200], range=[[1.5e4, 1.2e5], [0, 5e4]],
                      cmin=1, norm=LogNorm(), cmap=CMAP)
        # Gamma cut curve overlay: s2 = 3000 * 10^(s1/20000)
        s1_grid = np.linspace(1.5e4, 1.2e5, 500)
        gamma_line = 3.0e3 * np.power(10, s1_grid / 2.0e4)
        ax.plot(s1_grid, np.clip(gamma_line, 0, 5e4), 'r--', lw=1.5, alpha=0.7, label='gamma cut')
        ax.axhline(y=1.2e4, color='r', ls=':', lw=1, alpha=0.5)
        plt.colorbar(h[3], ax=ax, label='Counts')
        ax.legend(fontsize=7)
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('S2 Area [PE]')
    ax.set_title('S2 vs S1 — Gamma/Alpha Separation')
    ax.set_yscale('log')

    # ── Panel 2: S1 Area vs AFT (THE key diagnostic — alpha bands) ──
    ax = axes[0, 1]
    if 's1_area_fraction_top' in df.columns:
        h = ax.hist2d(df['s1_area'], df['s1_area_fraction_top'],
                      bins=[300, 200], range=[[1.5e4, 1.2e5], [0.0, 0.55]],
                      cmin=1, norm=LogNorm(), cmap=CMAP)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('Area Fraction Top')
    ax.set_title('S1 Area vs AFT — Alpha Band Structure')

    # ── Panel 3: S1 Area vs Drift Time ──
    ax = axes[0, 2]
    if 'drift_time' in df.columns:
        h = ax.hist2d(df['s1_area'], df['drift_time'] / 1e3,
                      bins=[200, 200], range=[[1.5e4, 1.2e5], [0, 500]],
                      cmin=1, norm=LogNorm(), cmap=CMAP)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('Drift Time [μs]')
    ax.set_title('S1 Area vs Drift Time')

    # ── Panel 4: Drift Time vs AFT (S1-S2 pair envelope) ──
    ax = axes[1, 0]
    if 'drift_time' in df.columns and 's1_area_fraction_top' in df.columns:
        h = ax.hist2d(df['s1_area_fraction_top'], df['drift_time'] / 1e3,
                      bins=[150, 200], range=[[0.0, 0.55], [0, 500]],
                      cmin=1, norm=LogNorm(), cmap=CMAP)
        # Overlay current polyFunc envelope bounds
        aft_grid = np.linspace(0.0, 0.55, 200)
        dt_lower = polyFunc(aft_grid, [-0.05, 2.0e6, -4.0e6, 0.0]) / 1e3
        dt_upper = polyFunc(aft_grid, [0.18, 2.0e6, -4.0e6, 0.0]) / 1e3
        ax.plot(aft_grid, dt_lower, 'r--', lw=1.5, alpha=0.7, label='lower bound')
        ax.plot(aft_grid, dt_upper, 'r--', lw=1.5, alpha=0.7, label='upper bound')
        ax.legend(fontsize=7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('Area Fraction Top')
    ax.set_ylabel('Drift Time [μs]')
    ax.set_title('Drift Time vs AFT — S1-S2 Pair Envelope')

    # ── Panel 5: Drift Time vs S2 Width (unphysical event removal) ──
    ax = axes[1, 1]
    if 'drift_time' in df.columns and 's2_range_50p_area' in df.columns:
        h = ax.hist2d(df['s2_range_50p_area'] / 1e3, df['drift_time'] / 1e3,
                      bins=[200, 200], range=[[0, 15], [0, 500]],
                      cmin=1, norm=LogNorm(), cmap=CMAP)
        # Overlay current unphysical lower bound: dt > polyFunc(s2w, [-2000, -50000, 0, 0.025])
        s2w_grid = np.linspace(0, 15e3, 200)
        dt_lower_s2w = polyFunc(s2w_grid, [-2.0e3, -0.5e5, 0.0, 0.025]) / 1e3
        ax.plot(s2w_grid / 1e3, dt_lower_s2w, 'r--', lw=1.5, alpha=0.7, label='lower bound')
        ax.legend(fontsize=7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S2 Range 50% Area [μs]')
    ax.set_ylabel('Drift Time [μs]')
    ax.set_title('Drift Time vs S2 Width — Unphysical Removal')

    # ── Panel 6: S1 Area vs r² (radial uniformity / FV r² cut) ──
    ax = axes[1, 2]
    if 'r2' in df.columns:
        h = ax.hist2d(df['s1_area'], df['r2'],
                      bins=[200, 200], range=[[1.5e4, 1.2e5], [0, 6000]],
                      cmin=1, norm=LogNorm(), cmap=CMAP)
        ax.axhline(y=FV_R2_MAX, color='red', ls='--', lw=1.5, alpha=0.7, label=f'r² < {FV_R2_MAX}')
        ax.legend(fontsize=7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('r² [cm²]')
    ax.set_title('S1 Area vs r² — Radial Uniformity')

    # ── Panel 7: z Distribution (FV z cut tuning) ──
    ax = axes[2, 0]
    if 'z' in df.columns:
        ax.hist(df['z'], bins=120, range=[-150, 10],
                histtype='stepfilled', color='steelblue', alpha=0.6, edgecolor='k', lw=0.5)
        ax.axvline(x=FV_Z_MIN, color='red', ls='--', lw=1.5, alpha=0.7, label=f'z > {FV_Z_MIN}')
        ax.axvline(x=FV_Z_MAX, color='red', ls='--', lw=1.5, alpha=0.7, label=f'z < {FV_Z_MAX}')
        ax.legend(fontsize=7)
    ax.set_xlabel('z [cm]')
    ax.set_ylabel('Counts')
    ax.set_yscale('log')
    ax.set_title('z Distribution — FV z Cut')

    # ── Panel 8: AFT Distribution (BiPo rejection tuning) ──
    ax = axes[2, 1]
    if 's1_area_fraction_top' in df.columns:
        ax.hist(df['s1_area_fraction_top'], bins=100, range=[0.0, 0.55],
                histtype='stepfilled', color='darkorange', alpha=0.6, edgecolor='k', lw=0.5)
        ax.axvline(x=0.12, color='red', ls='--', lw=1.5, alpha=0.7, label='AFT > 0.12')
        ax.legend(fontsize=7)
    ax.set_xlabel('Area Fraction Top')
    ax.set_ylabel('Counts')
    ax.set_yscale('log')
    ax.set_title('AFT Distribution — BiPo Rejection')

    # ── Panel 9: alt_s1_area vs s1_area (multi-scatter rejection) ──
    ax = axes[2, 2]
    if all(c in df.columns for c in ['alt_s1_area', 's2_area', 's1_area']):
        # Only show events with non-NaN alt_s1_area > 0
        ms_mask = df['alt_s1_area'].notna() & (df['alt_s1_area'] > 0)
        n_ms = ms_mask.sum()
        if n_ms > 0:
            # Color by s2_area to see the odd multi-scatter region (s2: 6e4-2e6)
            scatter = ax.scatter(df.loc[ms_mask, 's1_area'].values,
                       df.loc[ms_mask, 'alt_s1_area'].values,
                       c=df.loc[ms_mask, 's2_area'].values,
                       s=SCATTER_S, alpha=SCATTER_ALPHA, norm=LogNorm(), cmap='viridis')
            plt.colorbar(scatter, ax=ax, label='s2_area [PE]')
            # Overlay multi-scatter rejection region: s1 0-5.9e4, alt_s1 > 70, s2 6e4-2e6
            ms_rect = plt.Rectangle((0, 70), 5.9e4, 500, fill=False, edgecolor='red',
                                    ls='--', lw=1.5, alpha=0.7, label='MS rejection box')
            ax.add_patch(ms_rect)
        else:
            ax.text(0.5, 0.5, f'No events with alt_s1_area > 0\n(total: {len(df):,})',
                    ha='center', va='center', transform=ax.transAxes)
        ax.legend(fontsize=7)
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('alt_s1_area [PE]')
    ax.set_title('alt_s1_area vs S1 Area — Multi-Scatter')

    plt.tight_layout()
    plt.savefig('/project2/lgrandi/jietang/sr3_commissioning/s1s2_raw_data_overview.png',
                dpi=150, bbox_inches='tight')
    plt.show()
else:
    print('ERROR: df is None — cannot generate plots. Check data loading.')
'''

new_code_cell = {
    "cell_type": "code",
    "metadata": {},
    "source": [plot_code]
}

# Insert after cell 15 (the 'df' display cell), before cell 16 (Alpha Event Selection Cuts)
nb['cells'].insert(16, new_md_cell)
nb['cells'].insert(17, new_code_cell)

with open(r'C:\Users\13696\Desktop\xenon\Rn222\rn222_s1s2.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print(f'Done! Inserted 2 cells after cell 15.')
print(f'Total cells now: {len(nb["cells"])}')
