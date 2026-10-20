"""Half-height feature definitions preserved from Figure7 ver03.
Only the input curves are smoothed in ver04; episode geometry is unchanged.
"""
import numpy as np
import json
from scipy.signal import find_peaks
LABELS=[f'{c}.{s}' for c in range(1,9) for s in range(1,4)]

def summarize(values):
    observed = float(values[0])
    valid = np.isfinite(values[1:])
    q = np.quantile(values[1:][valid], [.025, .975]) if np.isfinite(observed) and valid.any() else [np.nan, np.nan]
    return dict(estimate=observed, ci95_low=float(q[0]), ci95_high=float(q[1]),
                valid_bootstrap_n=int(valid.sum()), valid_bootstrap_fraction=float(valid.mean()))


def temporal_features(y, peak_lo, peak_hi, bound_lo, bound_hi):
    r = np.arange(len(y))
    baseline = np.nanmin(y, axis=1)
    peak_index = np.argmax(y[:, peak_lo:peak_hi+1], axis=1) + peak_lo
    peak_value = y[r, peak_index]
    level = baseline + .5 * (peak_value-baseline)
    finite = np.isfinite(y).all(axis=1) & (peak_value > baseline)
    above = y >= level[:, None]
    ixs = np.arange(y.shape[1]-1)[None, :]
    out = {'peak': peak_index.astype(float)}
    out['peak'][~finite] = np.nan
    for event in ['rise', 'fall']:
        candidates = ((~above[:, :-1]) & above[:, 1:]) if event == 'rise' else (above[:, :-1] & ~above[:, 1:])
        if event == 'rise':
            candidates &= (ixs >= bound_lo) & (ixs < peak_index[:, None])
            indices = candidates.argmax(axis=1)
            censor = above[:, bound_lo]
        else:
            candidates &= (ixs >= peak_index[:, None]) & (ixs < bound_hi)
            indices = candidates.shape[1]-1-candidates[:, ::-1].argmax(axis=1)
            censor = above[:, bound_hi]
        with np.errstate(divide='ignore', invalid='ignore'):
            val = indices+(level-y[r, indices])/(y[r, indices+1]-y[r, indices])
        val[(~finite) | (~candidates.any(axis=1)) | censor] = np.nan
        out[event] = val
    return out, baseline, peak_value, level


def cut(x, y, lo, hi):
    keep = (x > lo) & (x < hi)
    return np.r_[lo, x[keep], hi], np.r_[np.interp(lo, x, y), y[keep], np.interp(hi, x, y)]


def interpolate_cross(x0, x1, y0, y1, level):
    return float(x0+(level-y0)*(x1-x0)/(y1-y0))


def spatial_feature(x, y, lo, hi):
    baseline = float(np.min(y))
    observed = np.flatnonzero((x >= lo) & (x <= hi))
    k = int(observed[np.argmax(y[observed])])
    peak = float(x[k]); peak_value = float(y[k]); level = baseline+.5*(peak_value-baseline)
    xx, yy = cut(x, y, lo, hi)
    above = yy >= level
    good = bool(peak_value > baseline)
    rise_candidates = np.flatnonzero((~above[:-1]) & above[1:] & (xx[:-1] < peak))
    fall_candidates = np.flatnonzero(above[:-1] & ~above[1:] & (xx[1:] > peak))
    rise_censored = bool(above[0] or not len(rise_candidates) or not good)
    fall_censored = bool(above[-1] or not len(fall_candidates) or not good)
    rise = np.nan if rise_censored else interpolate_cross(xx[rise_candidates[0]], xx[rise_candidates[0]+1], yy[rise_candidates[0]], yy[rise_candidates[0]+1], level)
    fall = np.nan if fall_censored else interpolate_cross(xx[fall_candidates[-1]], xx[fall_candidates[-1]+1], yy[fall_candidates[-1]], yy[fall_candidates[-1]+1], level)
    ids, props = find_peaks(yy, prominence=.15*max(float(np.ptp(yy)), 1e-12))
    prominent = [float(xx[i]) for i in ids]
    multip = any(abs(a-b) >= .04 for i, a in enumerate(prominent) for b in prominent[i+1:])
    component_starts = np.flatnonzero(above & np.r_[True, ~above[:-1]])
    allcross = []
    for event, indices in [('rise', np.flatnonzero(~above[:-1] & above[1:])), ('fall', np.flatnonzero(above[:-1] & ~above[1:]))]:
        for j in indices:
            allcross.append(dict(event=event, position_ap=interpolate_cross(xx[j], xx[j+1], yy[j], yy[j+1], level)))
    info = dict(window_lo=lo, window_hi=hi, baseline_raw=baseline, peak_intensity_raw=peak_value,
        half_height_raw=level, peak=peak, rise=rise, fall=fall,
        rise_censored=rise_censored, fall_censored=fall_censored,
        peak_at_window_edge=bool(k == observed[0] or k == observed[-1]),
        multi_peak_flag=multip, prominent_peak_positions=json.dumps(prominent),
        n_halfheight_components=int(len(component_starts)), disconnected_halfheight_support=bool(len(component_starts) > 1),
        normalized_half_height=.5)
    peaks = []
    for i, p in zip(ids, props['prominences']):
        peaks.append(dict(position_ap=float(xx[i]), intensity_raw=float(yy[i]),
            normalized_intensity=float((yy[i]-baseline)/(peak_value-baseline)), prominence_raw=float(p),
            is_principal_peak=bool(np.isclose(xx[i], peak)),
            separated_from_principal_by_at_least_4pct=bool(abs(xx[i]-peak) >= .04)))
    return info, peaks, allcross


