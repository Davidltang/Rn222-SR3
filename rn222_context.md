# Rn-222 Analysis — Global Context (Updated Jul 27)

## Data Paths

| File | Description |
|------|-------------|
| `/project2/lgrandi/jietang/sr3_commissioning/sr3_s1only_corrected_for_fit_pp.pkl` | S1-only corrected events (per-partition correction) |
| `/project2/lgrandi/jietang/sr3_commissioning/sr3_s1s2_corrected_for_fit_pp.pkl` | S1+S2 corrected events (per-partition correction) |
| `/project2/lgrandi/jietang/sr3_commissioning/sr3_s1only_pre_corr.pkl` | S1-only pre-correction (df_cut state) |
| `/project2/lgrandi/jietang/sr3_commissioning/sr3_s1s2_pre_corr.pkl` | S1+S2 pre-correction (df_cut state) |

## Current Cuts & Corrections (as of Jul 27)

### Common (both S1-only and S1+S2)
- Pre-correction: `range_50p_area ∈ [40, 200]`, `r2 < 3800`
- AFT fit range: `[0.12, 0.55]`, AFT reference: `0.30`
- Band selection: `RN_UPPER = [8.86e+04, -0.0600, -0.0900, -0.0600, 68000]`, `RN_LOWER = [7.30e+04, -0.0600, -0.0900, -0.0600, 58000]`
- 3-step correction: AFT → r²/AFT → r²/Area

### S1-only specific
- AFT cut: `[0.10, 0.55]`
- FV cut: `AFT ∈ [0.15, 0.455]`, `r2_aft_corrected < 3800`

### S1+S2 specific
- AFT cut: `[0.10, 0.55]`
- FV cut: `AFT ∈ [0.15, 0.455]`, `r2_aft_corrected < 75`

## Partitions (4-partition, Q5-driven)

```
P1: OLD Q5               May 1–11    original Q5 purifier (small sample for S1-only)
P2: ST707, Q5 removed    May 11–26   Q5 removed, ST707 cleaning, Rn rate decaying
P3: ST707, Q5 Rn cleaned May 26–Jun12 ST707, Q5 Rn depleted, lowest Rn rate (~4.5 uBq/kg)
P4: FRESH Q5             Jun 12–Jul2 new Q5 pellets online, Kr distillation active
```

## Purifier Operations Timeline

| Date | Event |
|------|-------|
| Apr 29 | Q5 first installed in LXePUR loop |
| **May 11** | Q5 removed, ST707 installed |
| May 21 | Second ST707 installed, Q5 vessel emptied |
| May 28 | New Q5 material shipped |
| Jun 5 | Online Kr distillation started (~18:00) |
| Jun 8-9 | Q5 refilled with new pellets, regenerated |
| Jun 11 | Q5 vessel reinstalled, cooldown |
| **Jun 12 ~18:00** | Fresh Q5 online in LXePUR loop |
| Jun 20 | LXe level 4→5 mm |
| Jun 25 | LXe level 5→5.7 mm |
| Jun 28 | DAQ room AC died, electrodes OFF → S1-only data |
| Jul 3 | Q5 removed, new ST707 installed (beyond data range) |

## Fit Model

- **Semi-constrained**: Rn-222/Po-218/Rn-220 ratio-constrained by energy, Po-216/Po-214 free
- Linear background, 16 params
- Asymmetric counting windows: Po-218 right ≤ Rn-220 left − 1σ, Po-218 ≥ +1.5σ

## Current Results Summary

| Partition | Rn-222 Rate | Fitted Mu | Notes |
|-----------|:----------:|:---------:|-------|
| P1 OLD Q5 | 49.6 | 52,790 ± 1,677 | Highest rate, most stable |
| P2 removed | 19.9 | 53,426 ± 2,541 | Rate decaying, large spread |
| P3 cleaned | 4.5 | 50,258 ± 3,442 | Lowest rate, gain unstable |
| P4 FRESH Q5 | 40.4 | 51,722 ± 2,000 | Lower than P1 (Kr distillation?) |

## Notebooks

| File | Purpose |
|------|---------|
| `rn222_fit_test.ipynb` | Main fit + rate analysis + all trend plots |
| `rn222_correction_test.ipynb` | Per-partition AFT centerline diagnosis |
| `rn222_s1only-analysis-new.ipynb` | S1-only correction pipeline |
| `rn222_s1only-analysis-full-0707.ipynb` | S1+S2 correction pipeline |
| `rn222_correction_test_context.md` | Context for correction test session |

## Plotting Conventions
- Source colors: S1-only=steelblue, S1+S2=darkorange
- Partition vlines: **black solid** (Q5 replaced, Q5 Rn cleaned, Q5 installed)
- Event vlines: **gray dashed** (Kr distillation, LXe level)
- capsize=0 (no error bar caps), legend fontsize scaled 1.5x from default
- Step plots: single continuous line, no source separation
