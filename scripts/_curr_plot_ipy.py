# ═══════════════════════════════════════════════════════════════
# PLOT FIT SPECTRUM — inspect any run's fit
# ═══════════════════════════════════════════════════════════════
# Re-fits the run and plots histogram + multi-Gaussian components.

def plot_run_fit(run_name, events_df, runs_df, cutoff_s,
                 guess_mu_before=GUESS_MU_BEFORE,
                 guess_mu_after=GUESS_MU_AFTER,
                 fit_min_before=50000, fit_min_after=48000):
    """
    Plot histogram + 5-Gaussian fit for one run.

    Parameters
    ----------
    run_name : str — e.g. '078900'
    events_df : DataFrame with 'time', 'area_2_corrected'
    runs_df : DataFrame from st.select_runs() with 'start_s' column
    cutoff_s : unix timestamp for before/after boundary
    """
    rrow = runs_df[runs_df['name'] == run_name]
    if len(rrow) == 0:
        print(f'Run {run_name} not found in runs_df')
        return None

    rrow = rrow.iloc[0]
    t0 = int(rrow['start_s'] * 1e9)
    t1 = int(float(rrow['end'].timestamp()) * 1e9)
    evts = events_df[(events_df['time'] >= t0) & (events_df['time'] < t1)]

    period = 'before' if rrow['start_s'] < cutoff_s else 'after'
    gm = guess_mu_before if period == 'before' else guess_mu_after
    fm = fit_min_before if period == 'before' else fit_min_after

    # Fit this run
    res = fit_one_run(evts, gm, fm, label=run_name)
    p_fit = res['p_fit']

    # ── Plot ──
    n_bins = res['n_bins'] if res['n_bins'] > 0 else 200
    hist_color = 'steelblue' if period == 'before' else 'darkorange'

    fig, ax = plt.subplots(figsize=(16, 7))

    ax.hist(evts['area_2_corrected'], bins=n_bins, range=[45e3, 80e3],
            color=hist_color, alpha=0.5, label=f'Data ({len(evts)} evts)')

    x_fine = np.linspace(45e3, 80e3, 2000)
    peak_colors = plt.cm.tab10(np.linspace(0, 1, N_PEAKS))

    if p_fit is not None:
        # Individual Gaussians
        for i in range(N_PEAKS):
            amp, mu, sig = p_fit[3*i], p_fit[3*i+1], p_fit[3*i+2]
            g = amp * np.exp(-0.5 * ((x_fine - mu) / sig)**2)
            ax.plot(x_fine, g, '--', color=peak_colors[i], lw=1.5, alpha=0.7,
                    label=f'{PEAK_NAMES[i]}  mu={mu:.0f}  sig={sig:.0f}')
            # +-2sigma window shading
            ax.axvspan(mu - 2*sig, mu + 2*sig, color=peak_colors[i], alpha=0.06)

        # Total fit
        y_tot = multi_gaussian(x_fine, *p_fit)
        ax.plot(x_fine, y_tot, 'r-', lw=2.5, label='Total fit')
        ax.axhline(y=p_fit[-1], color='gray', ls=':', lw=1.5,
                   label=f'bg = {p_fit[-1]:.1f}')

    # Initial guess markers
    for mu0, nm in zip(gm, PEAK_NAMES):
        ax.axvline(x=mu0, color='gray', ls=':', lw=0.8, alpha=0.4)

    ax.set_xlabel('AFT + r^2 Corrected Area [PE]', fontsize=12)
    ax.set_ylabel('Counts per bin', fontsize=12)
    ok_str = 'OK' if res['success'] else 'FAILED'
    ax.set_title(f'Run {run_name}  [{period}]  '
                 f'fit: {ok_str} ({res["status"]})  '
                 f'{len(evts)} events  |  {n_bins} bins',
                 fontsize=13)
    ax.legend(fontsize=7, ncol=2, loc='upper right')
    ax.grid(True, alpha=0.3)
    safe_tight(fig)
    plt.show()

    # ── Print table ──
    if p_fit is not None:
        print(f'\n{"="*70}')
        print(f'  Run {run_name} [{period}] — Fitted Parameters')
        print(f'{"="*70}')
        print(f'{"Peak":10s} {"mu [PE]":>12s} {"sig [PE]":>10s} {"Amp":>10s} {"sig/mu":>8s}')
        print(f'{"-"*70}')
        for i in range(N_PEAKS):
            amp, mu, sig = p_fit[3*i], p_fit[3*i+1], p_fit[3*i+2]
            print(f'{PEAK_NAMES[i]:10s} {mu:10.0f} {sig:10.0f} {amp:10.0f} {100*sig/mu:7.2f}%')
        print(f'{"bg":10s} {"-":>12s} {"-":>10s} {p_fit[-1]:10.1f}')
        print(f'{"="*70}')

    return res


print('plot_run_fit() defined.')
print('Usage: plot_run_fit("078900", rn_events, runs_sr3, cutoff_s)')
print()
print('# Quick-look — pick a successful run:')
print('# ok = df_fits[df_fits["success"]]["run"].values')
print('# plot_run_fit(np.random.choice(ok), rn_events, runs_sr3, cutoff_s)')
