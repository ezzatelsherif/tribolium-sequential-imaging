"""Original Theo handover smoothing, isolated from landmark/WDL logic.

Source: gt_mlpt_analysis/source_notebook.ipynb, zero-based cells4/6/8.
See original_smoothing_config.json for provenance and every spatial setting.
Default normalization matches the original workflow. Set normalize=False
only when input already has the intended full-series/profile normalization.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import savgol_filter

CONFIG = json.loads(Path(__file__).with_name('original_smoothing_config.json').read_text())


def window_length(n, frac=0.10, polyorder=3):
    """Verbatim original window formula, including full-profile sample count."""
    wl = int(max(polyorder + 2, (n * frac) // 2 * 2 + 1))
    return min(wl, n if n % 2 == 1 else n - 1)


def smooth_signal(y, frac=0.10, polyorder=3, axis=-1):
    y = np.asarray(y, dtype=float)
    n = y.shape[axis]
    if n < 5:
        return y.copy()
    return savgol_filter(y, window_length(n, frac, polyorder),
                         polyorder, axis=axis, mode='interp')


def normalize_range(y, axis=-1):
    y = np.asarray(y, dtype=float)
    low = np.min(y, axis=axis, keepdims=True)
    span = np.max(y, axis=axis, keepdims=True) - low
    if np.any(span <= 0):
        raise ValueError('Constant profile has no defined min-max normalization.')
    return (y - low) / span


def temporal(y, axis=-1, normalize=True):
    """Smooth stage means or batches of stage means; returns ordinal x and y.

    Inputs have24 stages along axis (default last). Output replaces this
    dimension with800 values. Bootstrap callers should process batches to
    limit memory, and resample embryos BEFORE computing each input mean.
    This function does not pool cohorts, align genes or resample embryos.
    """
    y = np.asarray(y, dtype=float)
    if y.shape[axis] != CONFIG['temporal']['stage_count']:
        raise ValueError('Original workflow requires24 ordered stage means.')
    if not np.isfinite(y).all():
        raise ValueError('Nonfinite stage mean: handle missing stages explicitly.')
    if normalize:
        y = normalize_range(y, axis)
    xn = np.linspace(0., 1., y.shape[axis])
    xd = np.linspace(0., 1., CONFIG['temporal']['dense_grid_size'])
    yy = PchipInterpolator(xn, y, axis=axis)(xd)
    yy = smooth_signal(yy, frac=0.10, polyorder=3, axis=axis)
    return xd * (len(xn) - 1), yy


def canonical_stage(stage):
    s = str(stage).strip().replace('Stage', '').replace('.', '_')
    if s.endswith('_0'):
        raise ValueError(f'Unrecognized stage{stage!r}.')
    return 'Stage' + s


def spatial_parameters(stage, embryo, gene):
    """Exact original setting. Missing genes error; no silent new settings."""
    g = str(gene).lower()
    if g == 'run':
        g = 'runt'
    s = canonical_stage(stage)
    emb = str(embryo).replace('EMB_', '').replace('Embryo', '').strip()
    try:
        return CONFIG['spatial']['stages'][s]['smooth'][emb][g].copy()
    except KeyError as exc:
        raise KeyError(f'No original WDL smoothing setting for{s}/embryo{emb}/{g}.') from exc


def spatial(x, y, stage, embryo, gene, normalize=True):
    """Return native normalized AP grid and smoothed full native profile.

    Deliberately does not crop or resample: direct Figure7 measurements retain
    full AP coordinates and their pre-existing episode definitions.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y):
        raise ValueError('Spatial input must be equal-length one-dimensional x/y.')
    if not (np.isfinite(x).all() and np.isfinite(y).all()):
        raise ValueError('Nonfinite spatial values require explicit handling.')
    if np.ptp(x) <= 0 or np.any(np.diff(x) <= 0):
        raise ValueError('Spatial x must strictly increase over a nonzero range.')
    x = (x - x.min()) / np.ptp(x)
    if normalize:
        y = normalize_range(y)
    kw = spatial_parameters(stage, embryo, gene)
    return x, smooth_signal(y, **kw)


temporal_smooth = temporal
spatial_smooth = spatial
