import json

NB = 'rn222_analysis_master.ipynb'
with open(NB, encoding='utf-8') as f:
    nb = json.load(f)


def get_src(i):
    return ''.join(nb['cells'][i]['source'])


def set_src(i, text):
    nb['cells'][i]['source'] = text.splitlines(keepends=True)


def replace(i, old, new):
    s = get_src(i)
    assert old in s, 'cell %d: anchor not found: %r' % (i, old[:80])
    set_src(i, s.replace(old, new, 1))


# corr2 after panel (band subset)
replace(41,
        "    after_hi = np.nanpercentile(d_fit['r2_aft_corrected'], 99)\n"
        "    if not np.isfinite(after_hi) or after_hi <= 0:\n"
        "        after_hi = 3.0 * r2_ref\n"
        "    after_hi = max(after_hi, 1.5 * r2_ref, 1.0)",
        "    after_hi = np.nanpercentile(d_fit['r2_aft_corrected'], 99)\n"
        "    if not np.isfinite(after_hi) or after_hi <= 0:\n"
        "        after_hi = 1.5 * r2_ref\n"
        "    after_hi = float(np.clip(after_hi, 1.5 * r2_ref, 3.0 * r2_ref))")

# corr2 full-ROI after panel
replace(41,
        "    _r2a_hi_full = np.nanpercentile(d_all['r2_aft_corrected'], 99.5)\n"
        "    if not np.isfinite(_r2a_hi_full) or _r2a_hi_full <= 0:\n"
        "        _r2a_hi_full = 3.0 * r2_ref",
        "    _r2a_hi_full = np.nanpercentile(d_all['r2_aft_corrected'], 99.5)\n"
        "    if not np.isfinite(_r2a_hi_full) or _r2a_hi_full <= 0:\n"
        "        _r2a_hi_full = 1.5 * r2_ref\n"
        "    _r2a_hi_full = float(np.clip(_r2a_hi_full, 1.5 * r2_ref, 3.0 * r2_ref))")

# FV ax1 x-range cap
replace(51,
        "_r2_hi = np.nanpercentile(df['r2_aft_corrected'], 99.5)\n"
        "if not np.isfinite(_r2_hi) or _r2_hi <= 0:\n"
        "    _r2_hi = 1.5 * R2_FV_MAX\n"
        "_r2_hi = max(_r2_hi, 1.2 * R2_FV_MAX, 1.0)",
        "_r2_hi = np.nanpercentile(df['r2_aft_corrected'], 99.5)\n"
        "if not np.isfinite(_r2_hi) or _r2_hi <= 0:\n"
        "    _r2_hi = 1.2 * R2_FV_MAX\n"
        "_r2_hi = float(np.clip(_r2_hi, 1.2 * R2_FV_MAX, 2.0 * R2_FV_MAX))")

# FV radial x-range cap
replace(51,
        "_r99 = np.nanpercentile(_r_all, 99.5) if len(_r_all) else 0\n"
        "if not np.isfinite(_r99) or _r99 <= 0:\n"
        "    _r99 = 1.5 * _r_fv\n"
        "_edges_r = np.linspace(0, max(1.5*_r_fv, _r99, 1.0), 80)",
        "_r99 = np.nanpercentile(_r_all, 99.5) if len(_r_all) else 0\n"
        "if not np.isfinite(_r99) or _r99 <= 0:\n"
        "    _r99 = 1.5 * _r_fv\n"
        "_r99 = float(np.clip(_r99, 1.5 * _r_fv, 2.0 * _r_fv))\n"
        "_edges_r = np.linspace(0, max(1.5*_r_fv, _r99, 1.0), 80)")

with open(NB, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
print('range caps added')
