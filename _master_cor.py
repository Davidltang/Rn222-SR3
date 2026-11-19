import json

NB = "rn222_analysis_master.ipynb"
nb = json.load(open(NB, encoding="utf-8"))
C = nb["cells"]

def setcell(i, src):
    C[i]["source"] = src.splitlines(keepends=True)
    C[i]["outputs"] = []
    C[i]["execution_count"] = None

def replace(i, old, new):
    s = "".join(C[i]["source"])
    assert old in s, f"cell {i}: missing:\n{old[:120]}"
    C[i]["source"] = s.replace(old, new).splitlines(keepends=True)

# ---------------- cell 8 : peak list (one line per peak, # = drop) ----------------
setcell(8, '''\
# ---- Fit components: ONE PEAK PER LINE. Comment a line (#) to drop that peak. ----
# tuple: (name, Q [MeV], mu0, (mu_lo, mu_hi), sig0, (sig_lo, sig_hi))
FIT_PEAKS = [
    ('Rn-222', 5.590, 51000, (49500, 53500),  800, (300, 2000)),
    ('Po-218', 6.115, 56000, (54500, 59000),  900, (300, 2000)),
    ('Rn-220', 6.288, 58500, (57000, 59500),  800, (250, 2000)),
    ('Po-216', 6.778, 63000, (61000, 65000),  900, (250, 2000)),
    ('Po-214', 7.834, 72000, (69000, 75000), 1100, (300, 2500)),
]

FIT_COMPONENTS    = [p[0] for p in FIT_PEAKS]
FIT_ENERGIES      = np.array([p[1] for p in FIT_PEAKS])
N_PEAKS           = len(FIT_COMPONENTS)
FIT_COMPONENT_CFG = {p[0]: {'mu0': p[2], 'mu_b': p[3], 'sig0': p[4], 'sig_b': p[5]}
                     for p in FIT_PEAKS}
E_RN222 = {p[0]: p[1] for p in FIT_PEAKS}.get('Rn-222', 5.590)

# ---- Signal peaks used for the RATE (must be present in FIT_PEAKS) ----
SIGNAL_NAMES = ['Rn-222', 'Po-218', 'Po-214']
SIGNAL_IDX   = [FIT_COMPONENTS.index(s) for s in SIGNAL_NAMES]

# ---- Rate constants ----
FMASS   = 3660.0    # fiducial mass [kg]
N_SIGMA = 2.0       # counting window half-width (x sigma)

print('Fit components:', FIT_COMPONENTS)
print('Signal peaks  :', SIGNAL_NAMES, '-> idx', SIGNAL_IDX)
print('FMASS=%.0f kg, N_SIGMA=%.1f' % (FMASS, N_SIGMA))
''')

# ---------------- cell 63 : drop the local FIT_COMPONENT_CFG (now from cell 8) ---
replace(63,
"""RN222_REF = 52000

# Per-component initial guess + bounds (mu, sigma in PE).
# Rn-220 / Po-216 get tighter mu bounds so they cannot collapse onto the
# signal peaks in low-statistics partitions.
FIT_COMPONENT_CFG = {
    'Rn-222': {'mu0': 51000, 'mu_b': (49500, 53500), 'sig0':  800, 'sig_b': (300, 2000)},
    'Po-218': {'mu0': 56000, 'mu_b': (54500, 59000), 'sig0':  900, 'sig_b': (300, 2000)},
    'Rn-220': {'mu0': 58500, 'mu_b': (57000, 59500), 'sig0':  800, 'sig_b': (250, 2000)},
    'Po-216': {'mu0': 63000, 'mu_b': (61000, 65000), 'sig0':  900, 'sig_b': (250, 2000)},
    'Po-214': {'mu0': 72000, 'mu_b': (69000, 75000), 'sig0': 1100, 'sig_b': (300, 2500)},
}""",
"""RN222_REF = 52000

# FIT_COMPONENT_CFG comes from section 0.3 (built from FIT_PEAKS) so the peak
# list is a single switch: comment a line in FIT_PEAKS to drop that peak.""")

# ---------------- cell 35 : corr1 params with per-partition override -----------
setcell(35, '''\
# Correction 1: exp_bump_model band [a, b, c, d, e] per partition.
# Global default, then per-partition overrides (write only the partitions you change).
_CORR1_DEFAULT = {'upper': [8.86e+04, -0.0600, -0.0900, -0.0600, 68000],
                  'lower': [7.30e+04, -0.0600, -0.0900, -0.0600, 58000],
                  'n_sigma': 2.0}

# Per-partition overrides. Keys: 'upper'/'lower' (full [a,b,c,d,e]),
# or just 'e_upper'/'e_lower' to change the 5th number only, or 'n_sigma'.
_CORR1_OVERRIDE = {
    # 'P5a: 0715-0725': {'e_upper': 68000, 'e_lower': 58000},
    # 'P5c: 0808-0822': {'e_upper': 68000, 'e_lower': 58000},
}

def _merge_corr1(default, ov):
    cfg = {'upper': list(default['upper']), 'lower': list(default['lower']),
           'n_sigma': default['n_sigma']}
    ov = ov or {}
    if 'upper' in ov: cfg['upper'] = list(ov['upper'])
    if 'lower' in ov: cfg['lower'] = list(ov['lower'])
    if 'e_upper' in ov: cfg['upper'][4] = ov['e_upper']
    if 'e_lower' in ov: cfg['lower'][4] = ov['e_lower']
    if 'n_sigma' in ov: cfg['n_sigma'] = ov['n_sigma']
    return cfg

CORR1_PARAMS = {p[0]: _merge_corr1(_CORR1_DEFAULT, _CORR1_OVERRIDE.get(p[0]))
                for p in PARTITIONS}
''')

# ---------------- cell 40 : corr2 params with per-partition override ----------
setcell(40, '''\
_CORR2_DEFAULT = {'area_lo': 40000, 'area_hi': 45000, 'poly_deg': 10,
                  'r2_min': 0, 'r2_max': 750, 'n_bins': 40}

# Per-partition overrides (any subset of keys).
_CORR2_OVERRIDE = {
    # 'P5a: 0715-0725': {'area_lo': 40000, 'area_hi': 45000},
}

CORR2_PARAMS = {p[0]: {**_CORR2_DEFAULT, **_CORR2_OVERRIDE.get(p[0], {})}
                for p in PARTITIONS}
''')

# ---------------- cell 45 : corr3 params with per-partition override ----------
setcell(45, '''\
# Band line:  area = k_slope * r2 + intercept   (upper/lower).
_CORR3_DEFAULT = {'k_slope': -110, 'intercept_upper': 75000, 'intercept_lower': 68000,
                  'r2_min': 0, 'r2_max': 100, 'r2_fit_max': 55, 'poly_deg': 6, 'n_bins': 30}

# Per-partition overrides (any subset of keys).
_CORR3_OVERRIDE = {
    # 'P5a: 0715-0725': {'k_slope': -110, 'intercept_upper': 75000, 'intercept_lower': 68000},
}

CORR3_PARAMS = {p[0]: {**_CORR3_DEFAULT, **_CORR3_OVERRIDE.get(p[0], {})}
                for p in PARTITIONS}
''')

json.dump(nb, open(NB, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("patched master: peak switch (cell 8), corr1/2/3 overrides (35/40/45), cell 63 de-dup")
