import json

with open(r'C:\Users\13696\Desktop\xenon\Rn222\rn222_s1s2.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

new_source = '''# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  RAW DATA VISUALIZATION — 3×3 diagnostic grid for cut tuning               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# Use these plots to manually adjust cut parameters in apply_alpha_cuts() below.
# Current cut boundaries are overlaid as dashed red lines where applicable.

if df is not None:
    fig, axes = plt.subplots(3, 3, figsize=(22, 18))

    # ── Panel 1: s1_area vs s2_area (gamma/electron removal) ──
    ax = axes[0, 0]
    if 's2_area' in df.columns:
        h = ax.hist2d(df['s1_area'], df['s2_area'],
                      bins=[200, 200], range=[[1.5e4, 1.2e5], [0, 5e4]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        # Gamma cut curve: s2 = 3000 * 10^(s1/20000)
        s1_grid = np.linspace(1.5e4, 1.2e5, 500)
        gamma_line = 3.0e3 * np.power(10, s1_grid / 2.0e4)
        ax.plot(s1_grid, np.clip(gamma_line, 0, 5e4), 'r--', lw=1.5, alpha=0.7)
        ax.axhline(y=1.2e4, color='r', ls=':', lw=1, alpha=0.5)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('S2 Area [PE]')
    ax.set_title('s1_area vs s2_area')
    ax.set_yscale('log')

    # ── Panel 2: drift_time vs s1_area_fraction_top (S1-S2 pair envelope) ──
    ax = axes[0, 1]
    if 'drift_time' in df.columns and 's1_area_fraction_top' in df.columns:
        h = ax.hist2d(df['s1_area_fraction_top'], df['drift_time'] / 1e3,
                      bins=[150, 200], range=[[0.0, 0.55], [0, 500]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        # Overlay polyFunc envelope bounds
        aft_grid = np.linspace(0.0, 0.55, 200)
        dt_lower = polyFunc(aft_grid, [-0.05, 2.0e6, -4.0e6, 0.0]) / 1e3
        dt_upper = polyFunc(aft_grid, [0.18, 2.0e6, -4.0e6, 0.0]) / 1e3
        ax.plot(aft_grid, dt_lower, 'r--', lw=1.5, alpha=0.7)
        ax.plot(aft_grid, dt_upper, 'r--', lw=1.5, alpha=0.7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('Area Fraction Top')
    ax.set_ylabel('Drift Time [μs]')
    ax.set_title('drift_time vs s1_area_fraction_top')

    # ── Panel 3: x vs y (spatial distribution) ──
    ax = axes[0, 2]
    if 'x' in df.columns and 'y' in df.columns:
        h = ax.hist2d(df['x'], df['y'],
                      bins=[200, 200], range=[[-70, 70], [-70, 70]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        # FV r boundary circle
        r_fv = np.sqrt(FV_R2_MAX)
        circle = plt.Circle((0, 0), r_fv, fill=False, edgecolor='red',
                            ls='--', lw=1.5, alpha=0.7)
        ax.add_patch(circle)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_title('x vs y')
    ax.set_aspect('equal')

    # ── Panel 4: s1_area vs s1_area_fraction_top (alpha band structure) ──
    ax = axes[1, 0]
    if 's1_area_fraction_top' in df.columns:
        h = ax.hist2d(df['s1_area'], df['s1_area_fraction_top'],
                      bins=[300, 200], range=[[1.5e4, 1.2e5], [0.0, 0.55]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('Area Fraction Top')
    ax.set_title('s1_area vs s1_area_fraction_top')

    # ── Panel 5: drift_time vs s2_range_50p_area (unphysical event removal) ──
    ax = axes[1, 1]
    if 'drift_time' in df.columns and 's2_range_50p_area' in df.columns:
        h = ax.hist2d(df['s2_range_50p_area'] / 1e3, df['drift_time'] / 1e3,
                      bins=[200, 200], range=[[0, 15], [0, 500]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        # Unphysical lower bound: dt > polyFunc(s2w, [-2000, -50000, 0, 0.025])
        s2w_grid = np.linspace(0, 15e3, 200)
        dt_lower_s2w = polyFunc(s2w_grid, [-2.0e3, -0.5e5, 0.0, 0.025]) / 1e3
        ax.plot(s2w_grid / 1e3, dt_lower_s2w, 'r--', lw=1.5, alpha=0.7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S2 Range 50% Area [μs]')
    ax.set_ylabel('Drift Time [μs]')
    ax.set_title('drift_time vs s2_range_50p_area')

    # ── Panel 6: r vs z (fiducial volume) ──
    ax = axes[1, 2]
    if 'z' in df.columns and 'r' in df.columns:
        h = ax.hist2d(df['z'], df['r'],
                      bins=[200, 200], range=[[-150, 10], [0, 70]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        # FV boundaries
        ax.axvline(x=FV_Z_MIN, color='red', ls='--', lw=1.5, alpha=0.7)
        ax.axvline(x=FV_Z_MAX, color='red', ls='--', lw=1.5, alpha=0.7)
        ax.axhline(y=np.sqrt(FV_R2_MAX), color='red', ls='--', lw=1.5, alpha=0.7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('z [cm]')
    ax.set_ylabel('r [cm]')
    ax.set_title('r vs z')

    # ── Panel 7: s1_area vs drift_time ──
    ax = axes[2, 0]
    if 'drift_time' in df.columns:
        h = ax.hist2d(df['s1_area'], df['drift_time'] / 1e3,
                      bins=[200, 200], range=[[1.5e4, 1.2e5], [0, 500]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('Drift Time [μs]')
    ax.set_title('s1_area vs drift_time')

    # ── Panel 8: S1 area spectrum (1D histogram) ──
    ax = axes[2, 1]
    ax.hist(df['s1_area'], bins=400, range=[1.5e4, 1.2e5],
            histtype='step', color='steelblue', lw=1.5)
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('Counts')
    ax.set_yscale('log')
    ax.set_title('S1 Area Spectrum')

    # ── Panel 9: s1_area vs r² (radial uniformity) ──
    ax = axes[2, 2]
    if 'r2' in df.columns:
        h = ax.hist2d(df['s1_area'], df['r2'],
                      bins=[200, 200], range=[[1.5e4, 1.2e5], [0, 6000]],
                      cmin=1, norm=LogNorm(), cmap='viridis')
        ax.axhline(y=FV_R2_MAX, color='red', ls='--', lw=1.5, alpha=0.7)
        plt.colorbar(h[3], ax=ax, label='Counts')
    ax.set_xlabel('S1 Area [PE]')
    ax.set_ylabel('r² [cm²]')
    ax.set_title('s1_area vs r²')

    plt.tight_layout()
    plt.savefig('/project2/lgrandi/jietang/sr3_commissioning/s1s2_raw_data_overview.png',
                dpi=150, bbox_inches='tight')
    plt.show()
else:
    print('ERROR: df is None — cannot generate plots. Check data loading.')
'''

nb['cells'][17]['source'] = [new_source]

with open(r'C:\Users\13696\Desktop\xenon\Rn222\rn222_s1s2.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print('Cell 17 rewritten successfully.')
print(f'No legend calls: {"ax.legend(" not in new_source}')
print(f'Uses viridis: {"viridis" in new_source}')
