"""
Peak-level position plugins for S1-only analysis.

TWO options, in order of preference:
  1. peak_s1_positions_cnn  — built-in straxen, S1 CNN (x,y,z in cm)
  2. naive_peak_positions   — custom plugin using area_per_channel × PMT positions

According to straxen docs (datastructure.html):
  - peak_s1_positions_cnn provides: x_s1_cnn, y_s1_cnn, z_s1_cnn [cm], time, endtime
  - peaks provides: area_per_channel (494 PMTs), loc_x_som, loc_y_som
  - peak_basics provides: area, area_fraction_top, range_50p_area, type, etc.
    (NO area_per_channel, NO x/y position)

For FV correction, prefer peak_s1_positions_cnn.  Fall back to naive_peak_positions
if the CNN is unavailable (e.g. old processing).
"""

import numpy as np
import strax
import straxen


# ═══════════════════════════════════════════════════════════════════════
# Option 1: built-in S1 CNN positions (PREFERRED)
# ═══════════════════════════════════════════════════════════════════════
# Just load 'peak_s1_positions_cnn' alongside 'peak_basics'.
# No plugin registration needed — it's standard straxen.
#
#   st.get_iter(run, ['peak_basics', 'peak_s1_positions_cnn'], ...)
#
# Columns: x_s1_cnn, y_s1_cnn, z_s1_cnn (float32, cm)


# ═══════════════════════════════════════════════════════════════════════
# Option 2: naive PMT-weighted position (FALLBACK)
# ═══════════════════════════════════════════════════════════════════════

class naive_peak_positions(strax.Plugin):
    """
    Compute naive (x, y) position for each peak using area_per_channel
    as weights on the known PMT positions.

    This is simpler than the CNN but available even without tensorflow.
    """
    depends_on = ('peaks',)
    provides = 'naive_peak_positions'
    data_kind = 'naive_peak_positions'

    dtype = [
        ('x_naive', np.float32, 'Naive X position (cm)'),
        ('y_naive', np.float32, 'Naive Y position (cm)'),
    ] + strax.time_fields

    __version__ = '0.0.1'

    def compute(self, peaks):
        weight_map = np.abs(peaks['area_per_channel']).astype(np.float64)
        row_sums = np.sum(weight_map, axis=1)
        row_sums = np.where(row_sums == 0, 1.0, row_sums)
        weight_map /= row_sums[:, None]

        result = np.ones(len(peaks), dtype=self.dtype)
        result['time'] = peaks['time']
        result['endtime'] = strax.endtime(peaks)

        pmt_x = straxen.pmt_positions().x.values.astype(np.float64)
        pmt_y = straxen.pmt_positions().y.values.astype(np.float64)

        result['x_naive'] = np.sum(weight_map * pmt_x[None, :], axis=1)
        result['y_naive'] = np.sum(weight_map * pmt_y[None, :], axis=1)

        return result


def register_plugin(st):
    """Register naive_peak_positions with a straxen context."""
    try:
        st.register(naive_peak_positions)
        print("Registered naive_peak_positions plugin.")
    except Exception as e:
        print(f"Plugin registration skipped: {e}")


# ═══════════════════════════════════════════════════════════════════════
# Merge helper — works with EITHER position source
# ═══════════════════════════════════════════════════════════════════════

def merge_position_into_basics(chunk_data, columns_to_keep, pos_source='cnn'):
    """
    Given a strax chunk dict with 'peak_basics' plus a position data type,
    merge position fields into the basics array and compute r2.

    Parameters
    ----------
    chunk_data : dict of str → np.ndarray
        Output from st.get_iter(run, ['peak_basics', 'peak_s1_positions_cnn'])
        or st.get_iter(run, ['peak_basics', 'naive_peak_positions']).
    columns_to_keep : list of str
        Which peak_basics columns to retain.
    pos_source : str
        'cnn' (uses x_s1_cnn, y_s1_cnn) or 'naive' (uses x_naive, y_naive).

    Returns
    -------
    merged : np.ndarray
        Structured array with basics columns + x, y, r2.
    """
    basics = chunk_data['peak_basics']

    if pos_source == 'cnn':
        pos = chunk_data['peak_s1_positions_cnn']
        xcol, ycol = 'x_s1_cnn', 'y_s1_cnn'
    else:
        pos = chunk_data['naive_peak_positions']
        xcol, ycol = 'x_naive', 'y_naive'

    # Verify time alignment
    if not (np.array_equal(basics['time'], pos['time']) and
            np.array_equal(basics['endtime'], pos['endtime'])):
        raise ValueError(f"peak_basics and position time fields disagree")

    # Build merged dtype
    available = [c for c in columns_to_keep if c in basics.dtype.names]
    keep_dtype = np.dtype([(c, basics.dtype[c]) for c in available])

    extra_dtype = np.dtype([
        ('x', np.float32),
        ('y', np.float32),
        ('r2', np.float32),
    ])

    merged_dtype = np.dtype(keep_dtype.descr + extra_dtype.descr)
    merged = np.zeros(len(basics), dtype=merged_dtype)

    for name in keep_dtype.names:
        merged[name] = basics[name]

    merged['x'] = pos[xcol]
    merged['y'] = pos[ycol]
    merged['r2'] = pos[xcol]**2 + pos[ycol]**2

    return merged


def load_positions_for_df(st, run_ids, df_existing, prefer_cnn=True):
    """
    Load positions for runs and merge with an existing peak_basics DataFrame.

    Tries peak_s1_positions_cnn first (if prefer_cnn=True), falls back to
    naive_peak_positions (registering the plugin if needed).

    Parameters
    ----------
    st : strax.Context
    run_ids : list of str
    df_existing : pd.DataFrame
        Must have time, endtime, run_name columns.
    prefer_cnn : bool
        Try the built-in S1 CNN positions first.

    Returns
    -------
    df : pd.DataFrame
        df_existing with x, y, r2 columns added.
    """
    import pandas as pd

    # Decide which position source to try
    sources = []
    if prefer_cnn:
        sources.append(('cnn', 'peak_s1_positions_cnn', 'x_s1_cnn', 'y_s1_cnn'))
    sources.append(('naive', 'naive_peak_positions', 'x_naive', 'y_naive'))

    for source_name, data_type, xcol, ycol in sources:
        if source_name == 'naive':
            register_plugin(st)

        pos_frames = []
        success = True
        for run_id in run_ids:
            run_str = run_id if isinstance(run_id, str) else f"{run_id:06d}"
            try:
                pos_arr = st.get_df(run_str, data_type, progress_bar=False)
                pos_df = pd.DataFrame(pos_arr)
                pos_df['run_name'] = run_str
                pos_frames.append(pos_df)
            except Exception as e:
                print(f"Position load failed for {run_str} ({source_name}): {e}")
                success = False
                break

        if not success:
            continue

        df_pos = pd.concat(pos_frames, ignore_index=True)
        df = df_existing.merge(
            df_pos[['time', 'endtime', xcol, ycol, 'run_name']],
            on=['time', 'endtime', 'run_name'],
            how='inner'
        )
        df['x'] = df[xcol]
        df['y'] = df[ycol]
        df['r2'] = df['x']**2 + df['y']**2
        print(f"Loaded positions via {source_name}: {len(df)} events (from {len(df_existing)} input)")
        return df

    print("ERROR: No position source available. FV correction will not work.")
    return df_existing
