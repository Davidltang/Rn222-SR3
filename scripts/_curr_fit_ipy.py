# ═══════════════════════════════════════════════════════════════
# PER-RUN FIT — fit a single run + batch fit all runs
# ═══════════════════════════════════════════════════════════════
# Returns success/failure label for each run.
# Depends on: multi_gaussian(), N_PEAKS, PEAK_NAMES (from fit cell above).

from scipy.optimize import curve_fit

GUESS_SIGMA = 800  # initial sigma guess (PE)

def fit_one_run(events, guess_mu, fit_min=50000, label=''):
    """
    Fit 5-Gaussian model to a single run's spectrum.

    Parameters
    ----------
    events : DataFrame with 'area_2_corrected' column
    guess_mu : array of 5 float — initial μ per peak
    fit_min : float — exclude < fit_min PE from fit
    label : str — run name for printout

    Returns dict:
        success   : bool
        status    : str — success_tight / success_medium / success_wide
                          / fit_failed / too_few_events
        p_fit     : array or None — fitted params [amp0,mu0,sig0, ..., bg]
        n_events  : int
        n_bins    : int — adaptive bin count used
    """
    n_ev = len(events)
    if n_ev < 50:
        print(f'  [{label:8s}] SKIP  — only {n_ev} events')
        return {'success': False, 'status': f'too_few_events({n_ev})',
                'p_fit': None, 'n_events': n_ev, 'n_bins': 0}

    # Adaptive binning — fewer bins for low-stat runs
    n_bins = min(300, max(80, n_ev // 5))

    counts, edges = np.histogram(
        events['area_2_corrected'], bins=n_bins, range=[45e3, 80e3]
    )
    centers = 0.5 * (edges[:-1] + edges[1:])

    m = (centers >= fit_min) & (centers <= 80e3)
    x_fit, y_fit = centers[m], counts[m]
    y_err = np.sqrt(np.maximum(y_fit, 1))

    # Try progressively wider mu bounds
    for width, tag in [(400, 'tight'), (700, 'medium'), (1200, 'wide')]:
        p0, lo, hi = [], [], []
        for mu0 in guess_mu:
            idx = np.argmin(np.abs(x_fit - mu0))
            amp0 = max(y_fit[idx] * 0.4, 3)
            p0.extend([amp0, mu0, GUESS_SIGMA])
            lo.extend([0, mu0 - width, 100])
            hi.extend([np.inf, mu0 + width, 3500])
        p0.append(1.0); lo.append(0); hi.append(np.inf)

        try:
            p_fit, _ = curve_fit(
                multi_gaussian, x_fit, y_fit, p0=p0, sigma=y_err,
                bounds=(lo, hi), maxfev=30000
            )
            # Sanity checks — reject if any mu/sigma hit bounds
            ok = True
            for i in range(N_PEAKS):
                mu_f  = p_fit[3*i + 1]
                sig_f = p_fit[3*i + 2]
                if abs(mu_f - guess_mu[i]) > width * 0.95:
                    ok = False
                if sig_f < 50 or sig_f > 4000:
                    ok = False
            if ok:
                print(f'  [{label:8s}] OK    ({tag:6s}, {n_ev:4d} evts, {n_bins:3d} bins)')
                return {'success': True, 'status': f'success_{tag}',
                        'p_fit': p_fit, 'n_events': n_ev, 'n_bins': n_bins}
        except Exception:
            continue

    print(f'  [{label:8s}] FAIL  ({n_ev:4d} evts)')
    return {'success': False, 'status': 'fit_failed',
            'p_fit': None, 'n_events': n_ev, 'n_bins': n_bins}


def fit_all_runs(runs_df, events_df, cutoff_s,
                 guess_mu_before=GUESS_MU_BEFORE,
                 guess_mu_after=GUESS_MU_AFTER,
                 fit_min_before=50000, fit_min_after=48000):
    """
    Fit every run. Returns DataFrame — one row per run.

    Columns: run, period, n_events, success, status,
             fit_mu_<isotope>, fit_sigma_<isotope>  (5 peaks each)
    """
    rows = []
    for _, run in runs_df.iterrows():
        name = run['name']
        t0   = int(run['start_s'] * 1e9)
        t1   = int(float(run['end'].timestamp()) * 1e9)
        evts = events_df[(events_df['time'] >= t0) & (events_df['time'] < t1)]

        if run['start_s'] < cutoff_s:
            gm, fm, per = guess_mu_before, fit_min_before, 'before'
        else:
            gm, fm, per = guess_mu_after, fit_min_after, 'after'

        res = fit_one_run(evts, gm, fm, label=name)

        row = {'run': name, 'period': per, 'n_events': res['n_events'],
               'success': res['success'], 'status': res['status']}
        if res['p_fit'] is not None:
            for i, iso in enumerate(PEAK_NAMES):
                row[f'fit_mu_{iso}']    = res['p_fit'][3*i + 1]
                row[f'fit_sigma_{iso}'] = res['p_fit'][3*i + 2]
        else:
            for i, iso in enumerate(PEAK_NAMES):
                row[f'fit_mu_{iso}']    = np.nan
                row[f'fit_sigma_{iso}'] = np.nan
        rows.append(row)

    df = pd.DataFrame(rows)
    n_ok = df['success'].sum()
    print(f'\n{"="*60}')
    print(f'  FIT SUMMARY: {n_ok}/{len(df)} runs successful ({100*n_ok/len(df):.1f}%)')
    print(f'{"="*60}')
    for st in sorted(df['status'].unique()):
        n = (df['status'] == st).sum()
        print(f'  {st:25s} {n:4d} runs')
    print(f'{"="*60}')
    return df


print('fit_one_run() and fit_all_runs() defined.')
print('Usage: df_fits = fit_all_runs(runs_sr3, rn_events, cutoff_s)')
