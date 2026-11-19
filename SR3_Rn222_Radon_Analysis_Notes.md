# XENONnT SR3 Commissioning — Rn-222 Alpha Background Rate Analysis Notes

**Author:** Jie Tang  
**Date:** 2026-07-07 (last update)  
**Computing environment:** University of Chicago Midway3 cluster  
**Software:** Python 3.11.8, strax 2.2.1, straxen 3.2.8 (daq-vanilla branch), cutax 2.3.1  

---

## 1. Motivation and Physics Goals

Measure the Rn-222 alpha background rate in the XENONnT TPC during SR3 commissioning. Rn-222 is a radioactive background source in liquid xenon (LXe), continuously emanated from detector materials. Its alpha decays produce characteristic signals in the TPC, and tracking its rate as a function of time and detector conditions is essential for:

- Monitoring the internal radon level evolution during commissioning
- Evaluating the effectiveness of the radon distillation system
- Providing input for the background model of the science runs
- Cross-validating two independent reconstruction modes (S1-only vs S1+S2)

Five alpha-emitting isotopes in the radon chain are analyzed simultaneously:

| Isotope | Alpha Energy (MeV) | Primary Decay |
|---------|---------------------|---------------|
| **Rn-222** | 5.49 | Direct alpha decay |
| **Po-218** | 6.00 | Rn-222 daughter |
| **Rn-220** | 6.29 | Thoron chain |
| **Po-216** | 6.78 | Rn-220 daughter |
| **Po-214** | 7.69 | Rn-222 granddaughter |

---

## 2. Reconstruction Modes

### 2.1 S1-only Mode (gate grid off)

Used when the gate grid is not operational during commissioning. Only S1 (primary scintillation) signals are recorded. S2 (ionization) signals are absent, which means:

- **No z-position from drift time** — depth is inferred from AFT (Area Fraction Top), the fraction of S1 light collected by the top PMT array
- **No gamma/electron rejection via S2/S1 ratio**
- **No BiPo tagging**

**Advantage:** Maximizes available data-taking periods, including those where the gate was not yet commissioned.

### 2.2 S1+S2 Mode (standard operation)

Full S1+S2 event reconstruction with drift time, providing:

- Precise z-position from drift time → better fiducial volume (FV) definition
- S2/S1 ratio cuts for gamma/electron background rejection
- Alt-S1 for multi-scatter rejection
- Time differences for BiPo (Bi-212 → Po-212, Bi-214 → Po-214) tagging

**Advantage:** Better background rejection, lower systematic uncertainties on the rate.

---

## 3. Data and Run Selection

### 3.1 Run Database Query

All runs queried from the straxen run database with criteria:

```python
run_mode = "tpc_bkg"           # TPC background runs
exclude_tags = ["abandon", "messy"]  # quality cuts
```

- **S1-only runs:** `include_tags = ["s1_only"]` → 652 runs
- **S1+S2 runs:** `exclude_tags = ["s1_only"]` → 620 runs

### 3.2 SR3 Run Coverage

| Mode | Run Range | Dates | # Runs |
|------|-----------|-------|--------|
| S1-only | 078417–083466 | 2026-05-03 to 2026-07-06 | 652 |
| S1+S2 | 078399–083533 | ~2026-05-12 to 2026-07-04 | 620 |

Total SR3 runs: **1272** (= 652 S1-only + 620 S1+S2)

### 3.3 Data Types Loaded

**S1-only analysis (`peak_basics`):**
- Columns: `time`, `endtime`, `area`, `area_fraction_top`, `range_50p_area`, `type`
- Pre-filter (alpha ROI): `type == 1`, `area ∈ [5e3, 150e3]` PE, `range_50p_area ∈ [40, 200]` ns
- Position: naive (PMT-weighted) or CNN (`peak_s1_positions_cnn`)
- S1-only analysis full dataset (post time cut): ~2.26M events across 374 runs

**S1+S2 analysis (`event_info`):**
- Columns: `s1_area`, `s2_area`, `s1_area_fraction_top`, `drift_time`, `alt_s1_area`, `x`, `y`, `z`, `cs1`, `cs2`, `run_number`
- Pre-filter: `s1_area > 15000` PE
- Alpha selection via 6 dedicated cuts (see §4.2)

### 3.4 Temporal Division

The analysis is split at **2026-06-08** (UTC) into two periods:

- **Before (2026-06-08):** Pre-injection/stable period
  - S1-only: 582 runs, ~57,064 Rn-222 events
  - Rn-222 rate: ~12–21 μBq/kg (consistent with SR0/SR1 baseline)
- **After (2026-06-08):** Post-Rn-injection period
  - S1-only: 373 runs, ~173,651 Rn-222 events
  - Rn-222 rate: ~36–42 μBq/kg (elevated, indicating Rn injection/distillation tests)

---

## 4. Event Selection Cuts

### 4.1 S1-only Cuts (pre-correction)

| # | Cut | Parameter | Purpose |
|---|------|-----------|---------|
| 1 | Alpha ROI | `type == 1`, `area ∈ [5e3, 150e3]` PE | Select S1-like peaks |
| 2 | S1 pulse width | `range_50p_area ∈ [40, 200]` ns | Reject muons/electronic noise |
| 3 | r² FV | `r2 < 3800 cm²` | Radial fiducial volume |
| 4 | AFT fiducial | `AFT ∈ [0.1, 0.55]` (deferred to after corrections) | Z-proxy fiducial |

**Note:** The AFT cut is **deferred** until after the correction pipeline so that the correction fits have access to the full AFT range for accurate centerline determination.

Cut efficiency: ~46.8M → ~32.8M events (70.1% pass pre-correction quality cuts).

### 4.2 S1+S2 Cuts (event-level)

A six-stage cut sequence, each self-contained with tunable parameters:

| # | Cut | Parameter | Purpose |
|---|------|-----------|---------|
| 1 | Gamma/electron removal | `s2 < 3e3 × 10^(s1/2e4)` AND `s2 < 1.2e4` PE | Reject ER backgrounds via S2/S1 ratio |
| 2 | Physical S1-S2 pairs | `drift_time` enveloped by `polyFunc(AFT, pars_low)` and `polyFunc(AFT, pars_high)` | Select events with consistent drift time vs AFT |
| 3 | BiPo rejection | `AFT > 0.12` | Remove BiPo events near top PMTs |
| 4 | Unphysical events | `drift_time > polyFunc(s2_range_50p_area, pars)` | Remove anomalous drift time vs S2 width |
| 5 | Multi-scatter rejection | `alt_s1 > 70` PE in S2 [60k, 2000k] & s1 [0, 59k] region | Tag multi-scatter candidates |
| 6 | Fiducial volume | `z ∈ [-120, -10] cm`, `r² < 3800 cm²` | Standard FV cut |

---

## 5. S1 Area Correction Pipeline (S1-only)

The S1 scintillation light collection efficiency varies with event position in the TPC due to geometric and optical effects. These position-dependent biases are corrected in a **three-step iterative procedure** applied to the S1-only data.

### 5.1 Correction Functions

**The bump model** (empirical, from Florian Jörg):
$$
\DeclareMathOperator{\AFT}{AFT}
f_{\rm bump}(\AFT; a,b,c,d,e) = a \cdot \exp\!\left(\frac{\AFT + b}{c}\right) \cdot \left[\,1 - \exp\!\left(\frac{\AFT + b}{d}\right)\,\right] + e
$$
where:
- `e` = high-AFT asymptote (area at top of TPC)
- `c` = tilt/curvature of the band (dominant shape parameter)
- `a` = overall stretch factor
- `d` = low-AFT knee sharpness
- `b` = horizontal shift

**Linear model** (for r² dependence):
$$
f_{\rm linear}(r^2; a, b) = a \cdot r^2 + b
$$

### 5.2 Correction Steps

#### Correction 1 — AFT (Depth) Correction

**Goal:** Remove the curvature of the Rn-222 band in Area-vs-AFT space, normalizing all event areas to a reference depth AFT = 0.30.

**Procedure:**
1. **Initial band isolation:** Use manual `exp_bump_model` boundaries to bracket the Rn-222 band:
   - Upper pars: `[8.86e4, -0.06, -0.09, -0.06, 68000]`
   - Lower pars: `[7.30e4, -0.06, -0.09, -0.06, 58000]`
   - Fit range: AFT ∈ [0.12, 0.55]
2. **Initial centerline fit:** Fit `exp_bump_model` to the median area in each AFT bin within the band
3. **Refinement:** Apply the initial correction, then fit a vertical band (±1σ) in corrected space, and re-fit
4. **Final correction:**
   ```
   area_aft_corrected = area × exp_bump_model(0.30, p_opt) / exp_bump_model(AFT, p_opt)
   ```

**Refined fit parameters (post-April-22 data):**
- `a = 4.61 × 10⁴`, `b = -7.84 × 10⁻³`, `c = -0.112`, `d = -0.414`, `e = 6.21 × 10⁴`
- Reference area at AFT = 0.30: **~62,100 PE** (depends on data subset)

#### Correction 2 — AFT → r² Correction

**Goal:** Remove AFT-dependent radial bias. Higher AFT (events near top of TPC) systematically maps to different r² distributions. We normalize all r² values to the expected r² at AFT = 0.30.

**Procedure:**
1. Select Rn-band events via vertical cut in `area_aft_corrected ∈ [40000, 45000]` PE
2. Bin by AFT (40 bins), compute **median r² per AFT bin**
3. Fit **degree-10 polynomial**: `r²(AFT) = c₀ + c₁·AFT + … + c₁₀·AFT¹⁰`
4. Apply correction:
   ```
   r2_aft_corrected = r² × poly(0.30) / poly(AFT)
   ```

**Fit range:** AFT ∈ [0.12, 0.55], r² ∈ [0, 750] cm²

#### Correction 3 — r² → Area Correction

**Goal:** Remove residual r² dependence in the AFT-corrected area. Even after depth correction, light collection varies radially.

**Procedure:**
1. Draw linear band boundaries `area = k × r2_aft_corrected + intercept` (k = -110)
2. Fit **degree-6 polynomial**: `area(r2_aft_corrected)` for r² < 55 cm²
3. Apply correction:
   ```
   area_2_corrected = area_aft_corrected × poly(0) / poly(r2_aft_corrected)
   ```
4. Clip correction factors to [0.5, 2.0] to prevent extreme values

### 5.3 Post-Correction Fiducial Cuts

After corrections, the AFT and FV fiducial cuts are applied:

| Cut | Parameter | Events passing |
|-----|-----------|----------------|
| AFT cut | `AFT ∈ [0.1, 0.55]` | 12.1M |
| FV cut | `AFT ∈ [0.13, 0.46]`, `r2_aft_corrected < 38` (or 50) | ~1.5M → 517k |

Final cut survival rate: ~3.2% of input events (dominated by the narrow Rn-222 band selection).

---

## 6. Rate Extraction

### 6.1 Multi-Gaussian Spectrum Fitting

The fully corrected S1 area spectrum (`area_2_corrected`) in the fiducial volume is fitted with a sum of **N Gaussians + background**:

$$
f(x) = \sum_{i=1}^{N} \frac{A_i}{\sqrt{2\pi\sigma_i^2}} \exp\!\left(-\frac{(x - \mu_i)^2}{2\sigma_i^2}\right) + \text{bg}
$$

where `i` indexes isotopes: Rn-222, Po-218, Rn-220, Po-216, Po-214 (and optionally Po-212).

**Two fit variants are used:**

| Variant | Description |
|---------|-------------|
| **Linked Gaussians** (production code) | Common sigma σ for all peaks, µ constrained by known energy ratios |
| **Free Gaussians** (notebook) | Independent µ and σ for each peak, full flexibility |

### 6.2 Rate Calculation

- Counting window: **±2σ** (or ±3σ) around fitted peak center
- Rate = (events in window) / (total livetime × fiducial mass)
- Fiducial mass: **MFV = 3660 kg** (liquid xenon)
- Overlapping peaks (Rn-220/Po-218/Bi-212) are deconvolved by Gaussian tail subtraction

### 6.3 Rate Results Summary (Full SR3 Dataset)

Combined S1-only + S1+S2 per-run rates, 955 runs total:

| Isotope | Before 2026-06-08 (μBq/kg) | After 2026-06-08 (μBq/kg) | Ratio A/B |
|---------|----------------------------|---------------------------|-----------|
| **Rn-222** | 21.2 ± 18.8 | 41.7 ± 8.1 | 1.97× |
| **Po-218** | 16.9 ± 14.5 | 38.2 ± 8.0 | 2.26× |
| **Rn-220** | 8.5 ± 15.2 | 9.4 ± 8.0 | 1.11× ∼ stable |
| **Po-216** | 0.9 ± 1.6 | 1.9 ± 1.0 | 2.11× |
| **Po-214** | 8.1 ± 7.5 | 16.2 ± 5.1 | 2.00× |

**Key observations:**
1. Rn-220 (thoron chain) rate is **stable** across the injection event (~8–9 μBq/kg), acting as an internal calibration reference
2. Rn-222 rate **nearly doubles** after 2026-06-08, indicating radon injection and/or LXe distillation tests
3. Po-218 follows its parent Rn-222 closely (2.26× increase), consistent with secular equilibrium expectation
4. Po-214 shows a similar increase (2.00×), though with larger scatter due to lower statistics

---

## 7. Correction Functions — S1+S2 Mode

The S1+S2 correction pipeline differs from S1-only due to the availability of drift time:

### 7.1 AFT Correction
Fit `exp2Func(AFT; a, b, c, d) = a × exp(AFT − b)/c + d` to s1_area in Area-vs-AFT space. Reference depth set at AFT = 0.30.

### 7.2 Drift Time Correction
Fit `exp2Func(drift_time)` to s1_area. Reference at drift_time = 300,000 ns (center of TPC).

### 7.3 Radial (r²) Correction
Fit `linearFunc(r²; a, b) = a × (1 + b × r²)` to s1_corrF. Reference at r² = 0.

Applied sequentially: AFT → r² → final `s1_corrF`.

### 7.4 BiPo Tagging (S1+S2 specific)

| BiPo Chain | S1 Selection [PE] | Time Diff [ns] | Physics |
|------------|-------------------|----------------|---------|
| BiPo-212 | s1_corrF ∈ [70k, 100k] | ∈ [600, 2700] | Bi-212 (β) → Po-212 (α, τ = 299 ns) |
| BiPo-214 | s1_corrF ∈ [60k, 70k] | ∈ [700, 400k] | Bi-214 (β) → Po-214 (α, τ = 164 µs) |

Lifecycle efficiency is computed from exponential fits to time differences:
- `expFunc(t; a, b, c, d) = a × exp(−(t − c)/b) + d`

---

## 8. Detector Parameters

| Parameter | Value | Note |
|-----------|-------|------|
| LXe density | 2.862 kg/cm³ | At TPC operating temperature |
| Fiducial mass (MFV) | 3660 kg | Standard FV: z ∈ [-120, -10] cm, r² < 3800 cm² |
| Fiducial volume | ~1.25 × 10⁶ cm³ = ~1250 L | |
| Runs | 1272 total (652 S1-only + 620 S1+S2) | |
| Run duration | ~10 min (standard) or ~1 h (long runs) | |
| Analysis cutoff date | 2026-06-08 UTC | Separates pre/post-Rn-injection periods |

---

## 9. Data Products

### 9.1 Intermediate Pickles (on Midway3)

| File | Description |
|------|-------------|
| `sr3_alpha_s1only.pkl` | S1-only alpha candidates with naive positions |
| `sr3_alpha_s1only_merged_0706.pkl` | Merged S1-only + S1+S2 alpha peaks |
| `sr3_s1s2_event_info.pkl` | S1+S2 event-level data |
| `sr3_s1only_corrected_for_fit.pkl` | Fully corrected S1-only events ready for fitting |

### 9.2 Final Rate Output

**`sr3_full_rn_rates_0707.csv`** — 4775 rows × 10 columns:
- `run`, `start_s`, `end_s`, `livetime_s`
- `isotope` ∈ {Rn-222, Rn-220, Po-216, Po-218, Po-214}
- `n_events`, `rate_ukgb` [μBq/kg], `rate_err_ukgb`
- `period` ∈ {before, after}
- `source` ∈ {S1-only, S1+S2}

### 9.3 Plots

- `rn_rate_evolution_full_0707_bef_and_aft.png` — Before/after 2D corrected area spectra
- `rn_rate_evolution_full_0707_s1only_and_s1s2.png` — S1-only vs S1+S2 comparison
- `time_evolution_rn_rates_0601.png` — Time trend plot (older version)
- `step_plot_rn_rates_0601.png` — Step plot of all 5 isotopes
- `1d_hist_bef_aft_0608.png` — 1D comparison before/after
- `Partition_fit_*_rates.png` — Partition-based fits
- `plots/3peaks/`, `plots/6peaks/` directories

---

## 10. Code Structure and Workflow

### 10.1 Notebooks (Analysis)

| Notebook | Purpose |
|----------|---------|
| `rn222_s1only-analysis-full-0707.ipynb` | **Main S1-only analysis** (post-May-2026 data). Full 3-step correction + multi-Gaussian rate fit + time variation |
| `rn222_s1only-analysis-new.ipynb` | S1-only analysis (pre-April-22 data). Same pipeline, tuned for earlier runs |
| `rn222_s1s2.ipynb` | **S1+S2 paired analysis** (100 latest runs). 6-stage alpha selection + 3-step correction + per-chunk rate fit |
| `rn222_fit_test.ipynb` | Per-run and partition-based fitting. Tests ratio-constrained vs free Gaussian models, 5 time partitions |
| `plot_full_rates.ipynb` | **Local plotting notebook.** Reads `sr3_full_rn_rates_0707.csv` and creates publication-quality time-trend plots (5-panel + step overlay) |
| `notebooks/rn222_s1only-0.ipynb` | Early data loader + single-Gaussian rate extraction |
| `notebooks/rn222_s1only-data_loader.ipynb` | Data loader only: registers plugin, streams `peak_basics`, merges with `naive_peak_positions` |
| `notebooks/rn222_s1only.ipynb` | Original S1-only (first version). Per-cut tuning interface |
| `notebooks/tune_bump_model.ipynb` | Interactive `exp_bump_model` parameter tuning tool |

### 10.2 Scripts

| Script | Purpose |
|--------|---------|
| `scripts/naive_peak_positions.py` | Strax plugin for PMT-weighted (x, y) position reconstruction (`naive_peak_positions`) |
| `scripts/fv_z_correction.py` | Correction pipeline: `apply_fv_correction()`, `apply_aft_correction()`, diagnostic plotting |
| `scripts/load_data.py` | Convenience data loader: `load_s1_data()` |
| `scripts/_curr_fit_ipy.py` | Per-run 5-Gaussian fitting (`fit_one_run()`, `fit_all_runs()`) |
| `scripts/_curr_plot_ipy.py` | Single-run fit visualization |
| `scripts/jupyter_exec.py` | WebSocket-based remote Jupyter kernel executor (for CLI interaction with notebooks) |
| `scripts/gen_download_cmds.py` | Bootstrap script to auto-download peaklet data for 284 S1-only runs |

### 10.3 Production Code (Reference)

| File | Lines | Purpose |
|------|-------|---------|
| `Reference/radon_analysis_v2_s1s2.py` | 2460 | Base framework: CLI, data loading, S1+S2 3-step correction + rate extraction |
| `Reference/radon_analysis_v2_s1only.py` | 3573 | S1-only extension: 4-step iterative correction, fiducial mass estimation, surface rate |

**Correction flow in production code (4 steps):**
1. AFT correction (exp2Func on area vs AFT) → `s1_corr1`
2. r² correction (quantile interpolation of AFT → r² profile) → `r2_corr`
3. r² → area correction (linearFunc) → `s1_corr2`
4. AFT re-correction on r²-corrected areas → `s1_corrF` (final)

### Key parameters tuned from r² polynomial fit (in notebooks) to quantile interpolation (in production code; documented reason: "polynomial fits get too bad for this application").

---

## 11. Weekly Commissioning Meeting Reports

11 weekly presentations documented in `nT Meeting/` (PDFs, 2026-04-30 through 2026-07-23), covering:

- Run-by-run rate evolution plots
- S1-only vs S1+S2 rate comparisons
- Pre/post-April-22 detector status
- Rn injection tests and distillation system performance
- BiPo tagging and detector stability checks

---

## 12. Key References

1. **SR0/SR1 Rn-222 analysis** — `Reference/sr0+sr1_radon-222_analysis.pdf`: Previous science-run analysis methodology (3-step correction framework for S1+S2 mode)
2. **Florian Jörg's S1-only alpha analysis** — `Reference/S1-only radon event selection.pdf`: Original alpha selection strategy for S1-only data
3. **SR3 S1-only alpha analysis overview** — `Reference/SR3-s1-only_alpha_analysis_in_xenonnt.pdf`: Preliminary SR3 S1-only results
4. **Production code** — `Reference/radon_analysis_v2_s1only.py` and `radon_analysis_v2_s1s2.py`: Full CLI-driven pipeline used for final rate measurement
5. **CLAUDE.md** — `CLAUDE.md`: Comprehensive project documentation including software stack, data paths, plugin details, and physics context (see inline for complete reference)

---

## 13. Appendix: Physics Constants

| Quantity | Symbol | Value | Units |
|----------|--------|-------|-------|
| Rn-222 alpha energy | Eα(Rn-222) | 5.49 / 5.59 | MeV |
| Po-218 alpha energy | Eα(Po-218) | 6.00 / 6.002 | MeV |
| Rn-220 alpha energy | Eα(Rn-220) | 6.29 / 6.288 | MeV |
| Po-216 alpha energy | Eα(Po-216) | 6.78 / 6.778 | MeV |
| Po-214 alpha energy | Eα(Po-214) | 7.69 / 7.687 | MeV |
| LXe density | ρXe | 2.862 | g/cm³ |
| Fiducial mass | MFV | 3660 | kg |
| Bi-212 half-life | τ(Bi-212) | 60.55 | min |
| Po-212 half-life | τ(Po-212) | 299 | ns |
| Bi-214 half-life | τ(Bi-214) | 19.9 | min |
| Po-214 half-life | τ(Po-214) | 164 | µs |

---

## 14. Systematic Uncertainties and Known Issues

1. **AFT as z-proxy in S1-only:** AFT provides only an approximate z position. The `exp_bump_model` fit accounts for band curvature but assumes the band centerline tracks light collection efficiency, not a change in the alpha energy spectrum itself.

2. **r²-AAFT polynomial fit quality:** The polynomial degree (5–10) is chosen empirically. In the production code, this was replaced by **quantile interpolation** because polynomial fits showed poor convergence behavior in some configurations.

3. **Peak overlap (Rn-220 vs Po-218):** These two peaks have only ~3% separation in energy. The multi-Gaussian fit deconvolves them, but systematic biases can arise depending on the choice of common sigma vs independent sigma.

4. **S1-only fiducial mass uncertainty:** The S1-only fiducial mass is estimated by comparing overlapping S1-only and S1+S2 runs. The `validS1S2Ratio` factor compensates for differences in efficiency, but this ratio depends on the stability of the comparison.

5. **Deadtime correction:** Applied from SCADA slow-control data. Not accounted for in the notebook-level analysis (handled in production code).

6. **Light yield drift:** The production code tracks light yield evolution via `lightYieldCorrPars`. In the notebook analysis, this is partially accounted for by fitting separate before/after periods.

---

*Document prepared for SR3 commissioning analysis archiving and reporting.*  
*Last compiled: 2026-07-07*
