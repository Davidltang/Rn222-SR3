# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## ⚠️ Current status (updated 2026-11-19) — read first

- **Only the S1-only pipeline exists and works.** The S1+S2 method is NOT in place yet; development starts 2026-11-19. `rn222_s1s2.ipynb` is exploratory. Do not suggest "just use S1+S2 / event_info x,y" as an existing option.
- The r² correction (Correction 3) needs `peaks` (`area_per_channel` → naive positions). Recent runs (e.g. ~091000) have no peaks/peaklets, so only `peak_basics` → AFT-only correction is possible there.
- Open question from colleagues (Olga Razuvaev, Caio Ishikawa, Nov 2026): fast Rn-222 monitoring without XY correction. Caio's AFT-only result is ~6 µBq/kg vs ~4 expected.
- Findings so far (scripts/fast_noxy/toy_aft_only.py, perbin_export npz):
  - The bump at ~0.92×Rn-222 is wall events (99% at r²>40). Ratio wall/Rn-222 rises from 0.15 (P1, P4) to 0.90 (P3, low Rn), so it biases low-rate periods most.
  - Rn-222, Po-218 and Po-214 share the same r² distribution.
  - Po-218/Rn-222 = 0.92–0.95 in steady state (P3–P6), stable within stats. Po-214/Rn-222 is NOT stable (0.34–0.50).
  - Toy: upper-side counting [μ, μ+1.5σ] works (few %) only if the AFT-only Rn-222 σ stays ≲ ~1.2k PE. If σ is ~2k, every spectral method fails.
  - `max_pmt` / `max_pmt_area` are in `peak_basics` and give a rough radius without peaks. This is an option that stays inside the S1-only pipeline.

## Project overview

Analysis of Rn-222 alpha background rates in the XENONnT dark matter experiment TPC, using SR3 commissioning data. The analysis exists in two Jupyter notebooks representing two reconstruction modes.

## Computing environment

This code runs on the **University of Chicago Midway3 cluster**, not locally. The software environment is:

```
/opt/XENONnT/anaconda/envs/XENONnT_el7.2025.07.2/
```

Key data paths:
- Scratch output: `/scratch/midway3/jietang/strax_data`
- Shared commissioning data: `/project/lgrandi-shared/commissioning/` (read-only)
- Personal work: `/project2/lgrandi/jietang/sr3_commissioning/`
- Custom plotting style: `/home/jietang/xenon_plot_style`

To activate the environment on Midway3:
```bash
source /opt/XENONnT/anaconda/envs/XENONnT_el7.2025.07.2/bin/activate
```

## Software stack (locked versions)

| Package | Version | Notes |
|---------|---------|-------|
| python  | 3.11.8  |       |
| strax   | 2.2.1   | Core streaming data processor |
| straxen | 3.2.8   | XENONnT plugins; custom branch `daq-vanilla` at `/home/jietang/straxen/straxen` |
| cutax   | 2.3.1   | Cut/selection framework at `/project2/lgrandi/xenonnt/software/cutax/v2.3.1/` |

Also: numpy, scipy, pandas, matplotlib, xgboost.

## Notebooks

**⚠️ Before any notebook edit:** Always re-read the file (or the relevant cells via grep/Python extraction) immediately before making changes. The user may have edited it since the last read, and stale reads cause silent overwrites of their changes. Use grep or Python JSON extraction to read only the target cells — never re-read the full notebook (it's a ~28k-token single-line JSON).

### `rn222_s1only.ipynb` — S1-only alpha analysis

Purpose: Measure Rn-222 rate using only S1 (scintillation) signals, without requiring S2 (ionization) coincidence.

Workflow:
1. **Run selection** — Queries `straxen.contexts.xenonnt_online` for `tpc_bkg` runs tagged `s1_only`, excluding `abandon`/`messy`. Filters to SR3 run range 078760–079371.
2. **Data loading** — Streams `peak_basics` with a strax-level selection string cutting on `area`, `range_50p_area`, and `type==1` to select alpha-like S1 signals. Saves as pickle at `/project2/lgrandi/jietang/sr3_commissioning/sr3_alpha_s1only.pkl`.
3. **Band isolation** — Uses Florian's empirical `exp_bump_model(a, b, c, d, e)` to manually bracket the Rn-222 band in area-vs-AFT space, then fits the band centerline via `scipy.optimize.curve_fit`.
4. **AFT correction** — Normalizes all event areas to a reference AFT=0.30, removing the z-dependent light collection bias.
5. **Rate calculation** — Fits a Gaussian to the corrected 1D area spectrum and counts events within ±3σ to compute the Rn-222 rate.

Key columns used: `time`, `endtime`, `area`, `area_fraction_top`, `range_50p_area`, `type`.

### `rn222_s1s2.ipynb` — S1+S2 paired analysis

Purpose: Rn-222 analysis using full S1+S2 event reconstruction for improved background rejection.

Workflow:
1. **Run selection** — Broader run selection (078399–078718), excluding `s1_only` tagged runs.
2. **Data loading** — Loads `peak_basics` for all peak types, filtering `type>0` post-load.
3. **Alpha selection** — Implements multi-stage cuts defined in `apply_alpha_cuts()`:
   - Gamma removal via S2 vs S1 area threshold
   - Drift-time vs AFT envelope (physical S1-S2 pair selection)
   - BiPo rejection (AFT lower bound at 0.12)
   - Unphysical event removal (drift-time vs S2 width)
   - Multi-scatter rejection
   - Fiducial volume (z in [-120, -10] cm, r < ~61.6 cm)
4. **Exploratory plotting** — 6-panel overview (area spectrum, area vs AFT, area vs pulse width, type distribution, S1-only and S2-only zoomed views).

Key columns used (in addition to S1-only columns): `range_90p_area`, `rise_time`, `tight_coincidence`, `n_channels`, `max_pmt`, `max_pmt_area`.

## Key physics concepts

- **AFT (Area Fraction Top)**: Proxy for z-position. Fraction of S1 light collected by the top PMT array. Higher AFT = interaction closer to the top of the TPC. Alpha bands appear as curved structures in area-vs-AFT because light collection efficiency varies with depth.
- **S1 vs S2**: S1 is prompt scintillation light, S2 is delayed electroluminescence from drifted ionization electrons. The S1-only mode skips S2 reconstruction, enabling runs where the gate grid isn't operating.
- **Rn-222**: A radioactive background source in liquid xenon. Its alpha decays produce a characteristic peak around ~41,000 PE (after AFT correction), distinct from Po-218 and Po-214 alphas at different energies.
- **The bump model** (`exp_bump_model`): An empirical function `a * exp((x+b)/c) * (1 - exp((x+b)/d)) + e` that traces the curved alpha band in area-vs-AFT space, used both for band selection and AFT correction.

## Run selection convention

Runs are selected from the straxen run database by:
- `run_mode = "tpc_bkg"` — TPC background runs
- `exclude_tags = "abandon", "messy"` — quality cuts
- `include_tags` or `exclude_tags` with `"s1_only"` depending on the analysis mode
- Run IDs are 6-digit zero-padded strings (e.g., `"078760"`)

## Remote JupyterLab access

The analysis runs on Midway3 via JupyterLab. Access it through an SSH tunnel:

```bash
ssh -N -f -L localhost:16812:10.50.250.134:16812 jietang@midway3.rcc.uchicago.edu
```

This requires interactive password + Duo 2FA (must be run in a proper terminal, NOT via `!` in Claude Code). Once the tunnel is up, the JupyterLab REST API at `http://localhost:16812` can be used to read/edit notebooks, start kernels, and execute code. The Jupyter instance uses password authentication (not token-based). See `../.claude/memory/jupyter-access.md` for auth details.

The remote workspace contains:
- `sr3_commissioning/` — Rn-222 analysis notebooks + production scripts
- `S2_only/`, `-S2xy/` — Other analysis directories
- `xenon_plot_style/`, `straxen/` — Custom code repos
- `strax_data/` — Local strax data cache

## Production analysis framework

The remote workspace at `sr3_commissioning/radon_analysis_v2_s1only.py` (261 KB) is a full production analysis pipeline with CLI interface. Its companion config is `radon_analysis_config_v2_s1only_monitor.py` (20 KB).

### Three-step S1 area correction (from production code)

The production framework applies corrections in this specific order:

1. **z-correction (CNN position)**: `s1_corr1 = area * ref / interpFunc(z_s1_cnn)`
2. **FV radial correction**: `s1_corr2 = area * ref / linearFunc(r²_corr)` — removes radial dependence with a linear function of r²
3. **AFT/z-proxy correction**: `s1_corrF = s1_corr2 * ref / exp2Func(area_fraction_top)` — removes residual AFT dependence

The correction parameters are derived by `correction_derivation()` and saved as pickle files, then loaded by downstream analysis steps.

### Peak-level data types (from straxen docs)

| Data type | Has `area_per_channel`? | Has position? | Weight |
|-----------|------------------------|---------------|--------|
| `peaks` | Yes (494 PMTs) | Only `loc_x_som`, `loc_y_som` (SOM grid, int16) | Heavy — includes waveforms |
| `peak_basics` | No | Only `loc_x_som`, `loc_y_som` (SOM grid) | Light — summary columns |
| `peak_s1_positions_cnn` | No | **`x_s1_cnn`, `y_s1_cnn`, `z_s1_cnn` [cm, float32]** — S1 CNN 3D position | Light |
| `peak_positions` | No | `x`, `y` [cm] — S2-based position (NOT available in S1-only runs) | Light |

**For S1-only FV correction, use `peak_s1_positions_cnn`** — it's a built-in straxen plugin (no registration needed) that provides x, y, z in cm directly from the S1 hit pattern via CNN.

Also available as fallback: `naive_peak_positions` (custom plugin using `area_per_channel` × PMT positions).

Key config parameters:
- `S1AFT_cutRange = [0.15, 0.455]` — AFT fiducial range (3660 kg FV)
- `r2_corr_cutRange = [0.0, 2.5E3]` — r² cut (cm²) for FV
- `z_cutRange = [-140.0, -10.0]` — z cut (cm)
- `s1_corrF_cutRange = [4.2E4, 5.6E4]` — Rn-222 peak window (PE)

## Local analysis code

### `naive_peak_positions.py`

A strax plugin (`class naive_peak_positions(strax.Plugin)`) that computes naive (x, y) position from `area_per_channel` weighted by PMT positions. This is the **prerequisite for FV correction** — without it, there is no r² to correct against.

Functions:
- `register_plugin(st)` — registers the plugin with a straxen context
- `merge_naive_position(chunk_data, columns_to_keep)` — merges `peak_basics` and `naive_peak_positions` from a strax chunk, computing r2_naive
- `load_naive_positions(st, run_ids, df_existing)` — loads positions separately and merges with an existing DataFrame (useful if you already have a saved pickle)
- `plot_position_check(df)` — diagnostic XY map and area-vs-r² plot

The plugin depends on `peaks` (needs `area_per_channel`) and provides `naive_peak_positions` with fields `x_naive`, `y_naive` (cm). Register it before data loading: `register_plugin(st)`.

### `fv_z_correction.py`

Imports from `naive_peak_positions`. Provides the correction pipeline:
- `load_naive_positions(st, run_ids, df_existing)` — re-exported for convenience
- `apply_fv_correction(df, r2_col='r2_naive')` — fits `linear_func(r²)` to area in a narrow AFT slice, divides out radial dependence
- `apply_aft_correction(df, area_col='area_fv_corrected')` — fits `exp_bump_model(AFT)` to the alpha band centerline, divides out depth dependence
- `plot_correction_comparison(df)` — side-by-side Area vs AFT at each correction stage
- `quantify_band_width(df)` — computes σ, IQR, FWHM to quantify band narrowing

**Apply corrections in this order:**
1. Register plugin + load data with positions → get `r2_naive`
2. `apply_fv_correction()` → `area_fv_corrected`
3. `apply_aft_correction()` → `area_z_corrected`
4. `plot_correction_comparison()` → verify the band gets narrower

## Plotting conventions

- **Default colormap for 2D histograms:** `viridis` (unless explicitly specified otherwise). Do not use `Blues`, `Reds`, or other single-hue colormaps for 2D plots unless the user asks for them.
- Use `matplotlib.colors.LogNorm()` for 2D histograms with `cmin=1`.
- Overlay cut boundaries as dashed red lines (`'r--'`, `lw=2`).

## Data persistence

Intermediate results are pickled to `/project2/lgrandi/jietang/sr3_commissioning/`:
- `sr3_alpha_s1only.pkl` — S1-only alpha candidates DataFrame
- `sr3_s1s2_alpha_peaks.pkl` — S1+S2 peak DataFrame
