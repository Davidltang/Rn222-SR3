//[[xenon:xenonnt:analysis:analysts_overview_page:jie_tang|David Tang]], jietang@link.cuhk.edu.cn//\\
20260728, created

==== Environment ====
<code>
Environment TAG XENONnT_el7.2025.07.2
 module version                                                                             path                     git
 python  3.11.8                                /opt/XENONnT/anaconda/envs/XENONnT_el7.2025.07.2/bin/python3                    None
  strax   2.2.1                                /opt/XENONnT/anaconda/envs/XENONnT_el7.2025.07.2/lib/python3.11/site-packages/strax
straxen   3.2.8                                                  /home/jietang/straxen/straxen branch:daq-vanilla
  cutax   2.3.1                                  /project2/lgrandi/xenonnt/software/cutax/v2.3.1/cutax                    None
</code>

==== Data ====
<code>
# S1-only alpha candidates (pre-correction)
PICKLE_PATH = "/project2/lgrandi/jietang/sr3_commissioning/sr3_s1only_pre_corr.pkl"

# S1+S2 alpha candidates (pre-correction)
PICKLE_PATH = "/project2/lgrandi/jietang/sr3_commissioning/sr3_s1s2_pre_corr.pkl"

# Fully corrected events for fitting (per-partition corrections)
PICKLE_PATH = "/project2/lgrandi/jietang/sr3_commissioning/sr3_s1only_corrected_for_fit_pp.pkl"
PICKLE_PATH = "/project2/lgrandi/jietang/sr3_commissioning/sr3_s1s2_corrected_for_fit_pp.pkl"
</code>

==== Link to the Code ====

The analysis is performed in Jupyter notebooks:

  * ''rn222_s1only-analysis-new.ipynb'' — Main S1-only analysis: run selection, data loading, 3-step correction pipeline, multi-Gaussian fitting, and per-partition rate extraction.
  * ''rn222_s1s2.ipynb'' — S1+S2 paired analysis: 6-stage alpha selection, 3-step correction, per-run rate fitting.
  * ''rn222_fit_test.ipynb'' — Per-run and partition-based fitting: semi-constrained vs free Gaussian models, 4-partition rate extraction, and time-trend plots.
  * ''rn222_correction_test.ipynb'' — Per-partition AFT centerline diagnosis across the 4 purifier periods.
  * ''plot_full_rates.ipynb'' — Reads the final rate CSV and produces publication-quality time-trend plots.

Supporting scripts:

  * ''scripts/naive_peak_positions.py'' — Strax plugin for PMT-weighted (x, y) position reconstruction. Registered via ''register_plugin(st)'', provides ''naive_peak_positions'' with ''x_naive'', ''y_naive'' (cm).
  * ''scripts/fv_z_correction.py'' — Correction pipeline: ''apply_fv_correction()'', ''apply_aft_correction()'', diagnostic plotting.

The production reference code from the SR0/SR1 analysis is ''Reference/radon_analysis_v2_s1only.py''.

--------------------------

===== Overview =====

The Rn-222 alpha emitter is a dominant background source in liquid xenon (LXe) dark matter detectors, continuously emanated from detector materials. Its decay chain includes several mono-energetic alpha emitters (Po-218, Rn-220, Po-216, Po-214) in the high-energy region. Since these alpha decays are physically linked to the beta-emitter Pb-214 in the same chain, measuring the alpha rate provides constraints on the Pb-214 background in the low-energy region of interest.

The goal of this analysis is to measure the Rn-222 alpha rate using S1-only data taken during SR3 commissioning, covering periods when the detector was in various purifier configurations. The S1-only mode operates without the gate grid, enabling data-taking even when the gate was not commissioned. This analysis also cross-validates against the S1+S2 reconstruction mode in overlapping runs.

<figure cut>
|{{:xenon:xenonnt:ting:rndecaychain_.png?800|}}|
^Fig.0: Rn-222 decay chain with alpha energies.^
</figure>

--------------------------

===== Run Selection =====

Runs are selected from the straxen run database:

<code>
st = straxen.contexts.xenonnt_online()

# All SR3 tpc_bkg runs
runs_sr3 = st.select_runs(
    run_mode="tpc_bkg",
    exclude_tags=["abandon", "messy"]
)
# → 1272 runs total

# S1-only runs (gate grid off)
runs_s1only = st.select_runs(
    run_mode="tpc_bkg",
    include_tags=["s1_only"],
    exclude_tags=["abandon", "messy"]
)
# → 652 runs (078417–083466, May 3 – Jul 2)

# S1+S2 runs (standard operation)
runs_s1s2 = st.select_runs(
    run_mode="tpc_bkg",
    exclude_tags=["abandon", "messy", "s1_only"]
)
# → 620 runs (078399–083533)
</code>

==== S1-only Data Loading ====

Data is loaded from ''peak_basics'' with an alpha ROI pre-selection:

<code>
columns_to_keep = ['time', 'endtime', 'area', 'area_fraction_top',
                   'range_50p_area', 'type']

selection_str = '(type == 1) & (area > 5e3) & (area < 150e3) & (area > 0) & (range_50p_area > 40) & (range_50p_area < 200)'
</code>

XY position is obtained using the ''naive_peak_positions'' plugin, which computes PMT-area-weighted (x, y) positions. The plugin is registered before data loading:

<code>
from scripts.naive_peak_positions import register_plugin
register_plugin(st)
</code>

The loaded data (with naive positions merged) is saved as a pickle file for reuse: ''/project2/lgrandi/jietang/sr3_commissioning/sr3_s1only_pre_corr.pkl''.

==== S1+S2 Data Loading ====

For cross-validation, S1+S2 data is loaded from ''event_info'' with a parallel pipeline:

<code>
columns_to_keep = ['time', 'endtime', 's1_area', 's2_area',
                   's1_area_fraction_top', 'drift_time',
                   'alt_s1_area', 'x', 'y', 'z', 'cs1', 'cs2']

# Post-load: s1_area > 15000 PE
</code>

The S1+S2 data is saved as ''/project2/lgrandi/jietang/sr3_commissioning/sr3_s1s2_pre_corr.pkl''.

--------------------------

===== Area Correction Pipeline (S1-only) =====

The S1 scintillation light collection efficiency varies with event position in the TPC due to geometric and optical effects. Three empirical corrections are applied sequentially to remove these position-dependent biases: Area vs AFT, AFT vs r², r² vs Area.

==== Correction 1: Area vs AFT (Depth Correction) ====

The correction function:

$f_{\rm bump}(x; a, b, c, d, e) = a \cdot \exp\!\left(\frac{x + b}{c}\right) \cdot \left[\,1 - \exp\!\left(\frac{x + b}{d}\right)\,\right] + e$

<code>
def exp_bump_model(x, a, b, c, d, e):
    """Florian's empirical model for the alpha band: area vs AFT."""
    return a * np.exp((x + b) / c) * (1 - np.exp((x + b) / d)) + e
</code>

**Procedure:**

- Isolate the Rn-222 band with manual //exp_bump_model// boundaries:
  * Upper: ''[8.86e4, -0.06, -0.09, -0.06, 68000]''
  * Lower: ''[7.30e4, -0.06, -0.09, -0.06, 58000]''

- Fit ''exp_bump_model'' to the median area in 50 AFT bins within AFT ∈ [0.12, 0.55].

- Apply initial correction: ''area_aft_corrected = area × exp_bump(0.30) / exp_bump(AFT)''.

- Redefine the Rn band as a **vertical cut** (±1σ) in area_aft_corrected space.

- Re-fit the centerline on the improved selection.

- Apply the refined correction to the full dataset.

<figure>
|{{:jie_tang:corr1_p4_s1s2.png?direct|}}|
^Fig.1: Left: Raw Area vs AFT with manual boundaries and refined centerline. Right: area_aft_corrected vs AFT showing the flattened band.^
</figure>

==== Correction 2: AFT vs r² (Removing AFT–r² Coupling) ====

The correction function:

$r^2(AFT) = \sum_{i=0}^{N} c_i \cdot AFT^i \quad (N = 10)$

<code>
def polynomial(x, *coeffs):
    """Polynomial: coeffs[0] + coeffs[1]*x + coeffs[2]*x² + …"""
    return sum(c * x**i for i, c in enumerate(coeffs))
</code>

**Procedure:**

- Select Rn-band events in correction 1's vertical band: ''area_aft_corrected ∈ [40000, 45000] PE''.

- In the fit range AFT ∈ [0.12, 0.55], r² ∈ [0, 750] cm², bin by AFT (40 bins) and compute median r².

- Fit a degree-10 polynomial: $r^2(AFT)$ using bins with ≥ 20 events.

- Correct all events: ''r2_aft_corrected = r² × r²(AFT=0.30) / r²(AFT)''.

<figure>
|{{:jie_tang:corr2_fit_band_s1only.png?direct|}}|{{:jie_tang:corr2_bef_aft_s1only.png?direct&600|}}|
^Fig.2a: AFT vs r² with degree-10 polynomial fit.^Fig.2b: Before/after — r² vs AFT showing removal of AFT dependence.^
</figure>

==== Correction 3: r² vs Area (Residual Radial Correction) ====

The correction function uses the same ''polynomial(x, *coeffs)'' with degree 6.

**Procedure:**

- Define linear band boundaries in r2_aft_corrected vs area_aft_corrected space: ''area = k × r² + intercept'' (k = −110).

- Select band events, bin in r2_aft_corrected, compute median area_aft_corrected in each bin.

- Fit degree-6 polynomial to $area(r^2)$ for r² < 55 cm².

- Correct all events: ''area_2_corrected = area_aft_corrected × poly(0) / poly(r2_aft_corrected)''.

- Clip correction factors to [0.5, 2.0].

<figure>
|{{:jie_tang:corr3_fit_band_s1only.png?direct&600|}}|{{:jie_tang:corr3_bef_aft_s1only.png?direct&600|}}|
^Fig.3a: r2_aft_corrected vs area_aft_corrected with band and polynomial fit.^Fig.3b: Before/after Correction 3 on the full dataset.^
</figure>

==== Correction Stages Overview ====

<figure>
|{{:jie_tang:corr_stages_s1only.png?direct&850|}}|
|{{:jie_tang:corr_stages_s1s2.png?direct&850|}}|
^Fig.4: Area vs AFT at each correction stage. Upper: S1-only. Lower: S1+S2.^
</figure>

--------------------------

===== Pre-Correction Quality Cuts =====

The following cuts are applied before the correction pipeline:

  * Alpha ROI pre-selection in the strax selection string: ''type == 1'', ''area ∈ [5e3, 150e3] PE'', ''range_50p_area ∈ [40, 200] ns''.
  * ''r² < 3800 cm²'' — loose radial cut that removes < 0.01% of events.

Note: The r² cut is not strictly necessary but is retained to exclude any events with clearly unphysical positions at the TPC edge. The AFT fiducial cut is deferred until after corrections.

--------------------------

===== Fiducial Volume Cut (Post-Correction) =====

After all three area corrections are applied, the fiducial volume is selected:

  * AFT ∈ [0.13, 0.46]
  * r2_aft_corrected < 3800 cm² ( < 75 cm² for S1+S2 data)

<figure>
|{{:jie_tang:fv_cut_s1s2.png?direct&800|}}|
^Fig.5: FV cut diagnostic: AFT 1D, r² distribution, Area vs AFT, and AFT vs r2_aft_corrected with FV box.^
</figure>

--------------------------

===== Multi-Gaussian Spectrum Fitting =====

The fully corrected S1 area spectrum (''area_2_corrected'') within the fiducial volume is fitted with a sum of 5 Gaussians plus a constant background:

$$
f(x) = \sum_{i=1}^{5} A_i \exp\!\left(-\frac{(x - \mu_i)^2}{2\sigma_i^2}\right) + \text{bg}
$$

The five alpha-emitting isotopes and their initial guess PE positions (scaled from the Rn-222 anchor by known alpha energies) are:

| Isotope | Eα (MeV) | S1-only Guess [PE] | S1+S2 Guess [PE] |
|---------|----------------|---------------------|-------------------|
| Rn-222 | 5.49 | 52,000 | 50,000 |
| Po-218 | 6.00 | 56,830 | 54,645 |
| Rn-220 | 6.29 | 59,577 | 57,285 |
| Po-216 | 6.78 | 64,218 | 61,748 |
| Po-214 | 7.69 | 72,837 | 70,036 |

The fit is performed independently for each run using ''scipy.optimize.curve_fit'' with adaptive binning (min(300, max(80, n_events/5)) bins) and progressive mu bounds (±400 → ±700 → ±1200 PE) to maximize convergence.

<figure>
|{{:jie_tang:fit_spectra_s1only.png?direct&700|}}|
^Fig.6: 5-Gaussian fit with component peaks and total fit overlay (S1-only data).^
</figure>

==== Fitted μ Evolution ====

The per-run fitted μ (mean) for each isotope is monitored as a function of time to track light yield stability:

<figure>
|{{:jie_tang:fitted_mean_evolution_full_0707.png?direct&800|}}|
^Fig.7: Per-run fitted μ for all 5 isotopes. S1+S2 fits show instability during May (shift in μ), indicating light yield drift in that period.^
</figure>

--------------------------

===== Rate Calculation =====

The event rate for each isotope is computed from the per-run fit results:

$$
R = \frac{N_{\rm events}}{\Delta t \times M_{\rm FV}}
$$

where events are counted within ±2σ of the fitted peak center position for each run. The fiducial mass is $M_{\rm FV} = 3660$ kg.

Error contributions include Poisson counting uncertainty, the uncertainty from the Gaussian fit, and the binomial uncertainty from comparing S1-only to S1+S2 fiducial volume selection (see the SR0/SR1 analysis note for the FV estimation method).

==== Partition-Based Rate Summary ====

The full SR3 data-taking period is divided into 4 partitions based on the Q5 purifier status:

| Partition | Dates | Key Event | Rn-222 Rate [μBq/kg] | Notes |
|-----------|-------|-----------|----------------------|-------|
| P1: OLD Q5 | May 1–11 | Original Q5 purifier | 49.6 ± 1.7 | Highest and most stable rate |
| P2: removed | May 11–26 | Q5 removed, ST707 cleaning | 19.9 ± 2.5 | Rate decaying, large spread |
| P3: cleaned | May 26–Jun 12 | Q5 Rn-depleted | 4.5 ± 3.4 | Lowest rate, gain unstable |
| P4: FRESH Q5 | Jun 12–Jul 2 | Fresh Q5 + Kr distillation | 40.4 ± 2.0 | Lower than P1 (Kr distillation effect?) |

--------------------------

===== Time Evolution of Rn Rates =====

The Rn-222 rate trend is consistent with the known Q5 purifier transitions: highest with the old Q5, decaying after its removal, reaching a minimum after Rn cleaning, and rising again when fresh Q5 pellets come online. The Rn-220 (thoron) rate remains stable throughout, serving as an internal calibration reference.

A notable observation is that S1+S2 rates and light yield (fitted μ) were unstable during May, with the fitted peak positions shifting by ~5%. This is visible in Fig.7 and suggests a gain drift in that period.

<figure>
|{{:jie_tang:rn_rate_evolution_full_0707_s1only_s1s2.png?direct|}}|
^Fig.8: Errorbar scatter of all 5 isotope rates vs time (μBq/kg). S1-only (steelblue) and S1+S2 (darkorange) overlaid.^
</figure>

<figure>
|{{:jie_tang:rn_rate_evolution_full_0707_step.png?direct|}}|
^Fig.9: Rn-222 rate evolution step plot for the full SR3 period.^
</figure>

<figure>
|{{:jie_tang:image_1_.png?direct|}}|{{:jie_tang:step_plot_rn_rates_p4.png?direct|}}|
^Fig.10a: P3 zoom — without Q5 (Rn-depleted period).^Fig.10b: P4 zoom — with Q5 in place.^
</figure>

--------------------------

===== Systematic Uncertainties =====

  * **AFT as z-proxy**: AFT provides only an estimate of light collection, not true geometric depth. The bump model captures band curvature but residual position-dependent biases may remain.
  * **r²-AFT polynomial fit:** The degree-10 polynomial is chosen empirically. In the production reference code, this step is replaced by quantile interpolation for improved stability.
  * **Peak overlap (Rn-220 vs Po-218):** These peaks are separated by only ~3% in energy. The deconvolution is sensitive to the choice of fit model (common vs independent sigma).
  * **Light yield drift:** The fitted mu values shifted by ~5% across partitions (see Fig.7). Partition-dependent anchor points partially absorb this effect.

===== Reference =====
<table reference>
^Note^Author(s)^Date^
|[[xenon:xenonnt:analysis:dominick:xenonnt_sr0_radon_corrections|XENONnT SR0 radon analysis — signal corrections]]|Dominick|Mar 30th, 2022|
|[[xenon:xenon1t:jorg:xent_lxe_alpha_study_4|Closer look into alpha rate evolution in XEnT]]|Florian|Dec 17th, 2020|
|[[xenon:xenonnt:lin:sr1:240717|SR0+SR1 Radon-222 Analysis]]|Ting|Jul 17th, 2024|
</table>
