"""
Minimal data loader for S1-only Rn-222 analysis.
Loads peak_basics + naive_peak_positions, returns a DataFrame with r2.

Usage on Midway3:
    from load_data import load_s1_data
    df = load_s1_data(st, run_ids, max_runs=10)
"""

import numpy as np
import pandas as pd
import strax
import straxen
from tqdm.notebook import tqdm


class naive_peak_positions(strax.Plugin):
    """
    Naive (x,y) from area_per_channel weighted by PMT positions.
    Does NOT set data_kind → inherits from depends_on ('peaks'),
    so strax can auto-merge with other peak-kind data types.
    """
    depends_on = ('peaks',)
    provides = 'naive_peak_positions'
    # NO data_kind override — inherit from peaks (tutorial §3)

    dtype = [
        ('x_naive', np.float32, 'Naive X position (cm)'),
        ('y_naive', np.float32, 'Naive Y position (cm)'),
    ] + strax.time_fields
    __version__ = '0.0.2'

    def compute(self, peaks):
        weight_map = np.abs(peaks['area_per_channel']).astype(np.float64)
        row_sums = np.sum(weight_map, axis=1)
        row_sums = np.where(row_sums == 0, 1.0, row_sums)
        weight_map /= row_sums[:, None]
        pmt_x = straxen.pmt_positions().x.values.astype(np.float64)
        pmt_y = straxen.pmt_positions().y.values.astype(np.float64)
        result = np.ones(len(peaks), dtype=self.dtype)
        result['time'] = peaks['time']
        result['endtime'] = strax.endtime(peaks)
        result['x_naive'] = np.sum(weight_map * pmt_x[None, :], axis=1)
        result['y_naive'] = np.sum(weight_map * pmt_y[None, :], axis=1)
        return result


def load_s1_data(st, run_ids, max_runs=None, progress_bar=True):
    """
    Load peak_basics + naive_peak_positions, return DataFrame with r2.

    The plugin inherits data_kind from 'peaks', same as peak_basics,
    so strax auto-merges them — no allow_multiple needed.
    """
    # Register (idempotent — safe to call multiple times)
    try:
        st.register(naive_peak_positions)
    except Exception:
        pass  # already registered

    keep_cols = ['time', 'endtime', 'area', 'area_fraction_top',
                 'range_50p_area', 'type']

    if max_runs:
        run_ids = run_ids[:max_runs]

    frames = []
    n_loaded, n_failed = 0, 0

    for run in tqdm(run_ids, disable=not progress_bar):
        run_str = run if isinstance(run, str) else f"{run:06d}"

        if not st.is_stored(run_str, "peak_basics"):
            continue

        chunks = []
        try:
            # Same data_kind → strax auto-merges basics + positions
            for chunk in st.get_iter(
                run_str,
                ["peak_basics", "naive_peak_positions"],
                progress_bar=False,
            ):
                if len(chunk.data) > 0:
                    # chunk.data is the auto-merged array (all columns)
                    chunks.append(chunk.data)
            n_loaded += 1
        except Exception as e:
            n_failed += 1
            if n_failed <= 3:
                print(f"Failed on {run_str}: {e}")
            continue

        if not chunks:
            continue

        arr = np.concatenate(chunks)
        # Build DataFrame: keep basics columns + add r2 from naive positions
        available = [c for c in keep_cols if c in arr.dtype.names]
        run_df = pd.DataFrame(arr[available])
        run_df['x'] = arr['x_naive']
        run_df['y'] = arr['y_naive']
        run_df['r2'] = run_df['x']**2 + run_df['y']**2
        run_df['run_name'] = run_str
        frames.append(run_df)

    print(f"Runs loaded: {n_loaded}, failed: {n_failed}")

    if not frames:
        raise RuntimeError("No data loaded!")

    df = pd.concat(frames, ignore_index=True)

    # Basic alpha cuts
    n_all = len(df)
    df = df[df['type'] == 1]
    df = df[(df['area'] >= 5e3) & (df['area'] <= 150e3)]
    df = df[(df['range_50p_area'] >= 40) & (df['range_50p_area'] <= 200)]
    print(f"S1 alpha candidates: {len(df):,}  (filtered from {n_all:,})")

    return df
