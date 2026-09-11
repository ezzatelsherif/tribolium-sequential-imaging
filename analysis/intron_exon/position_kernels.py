"""Analytic age kernels for the mapped intron/exon observation model.

All durations use the same stage-coordinate units as the input ages. Signal is
uniform per targeted base, normalized to one when all targets are present.
Mapped coordinates are fixed, and T is the duration from coordinate zero to
``transcript_end_bp``. The target intron disappears following an exponential
processing wait after its 3' site. Exonic RNA survives until completion and
target-intron processing, followed by an exponential decay wait. The model
does not track other introns or individual HCR oligonucleotide occupancies.

``collapsed`` moves all transcription and probe-appearance events to age zero
while keeping the same processing and decay rules. It is a matched biological
comparison with the same processing and degradation rules.

Integrals are analytic. The removable singularity at equal processing and
decay rates uses an analytic series evaluated to machine precision.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.special import gammainc


def _check(T, splice_mean, decay_mean, mode):
    if mode not in ("mapped", "collapsed"):
        raise ValueError("mode must be 'mapped' or 'collapsed'")
    if not all(np.isfinite(x) and x > 0 for x in (splice_mean, decay_mean)):
        raise ValueError("splice_mean and decay_mean must be finite and positive")
    if not np.isfinite(T) or T < 0 or (mode == "mapped" and T == 0):
        raise ValueError("T must be positive in mapped mode or nonnegative in collapsed mode")


def _exp_moment(b, rate, order):
    """Integral of x**order exp(-rate*x) from zero to nonnegative b."""
    if order == 0:
        return -np.expm1(-rate * b) / rate
    return math.factorial(order) * gammainc(order + 1, rate * b) / rate ** (order + 1)


def _sumexp_moments(b, rate1, rate2):
    """First two unnormalized moments of survival of Exp(r1)+Exp(r2)."""
    lo, hi = min(rate1, rate2), max(rate1, rate2)
    diff = hi - lo
    j0 = _exp_moment(b, lo, 0)
    j1 = _exp_moment(b, lo, 1)
    if diff / lo <= 1e-3:
        # H(x)=exp(-lo*x) [1 + lo*(1-exp(-diff*x))/diff].
        # Integrated terms have ratio bounded asymptotically by diff/lo.
        f0, f1 = j0.copy(), j1.copy()
        coefficient = lo
        for n in range(8):
            if n:
                coefficient *= -diff / (n + 1)
            f0 += coefficient * _exp_moment(b, lo, n + 1)
            f1 += coefficient * _exp_moment(b, lo, n + 2)
        return f0, f1
    return (j0 + lo * (j0 - _exp_moment(b, hi, 0)) / diff,
            j1 + lo * (j1 - _exp_moment(b, hi, 1)) / diff)


def _sumexp_survival(b, rate1, rate2):
    lo, hi = min(rate1, rate2), max(rate1, rate2)
    diff = hi - lo
    if diff == 0:
        correction = lo * b
    else:
        correction = lo * (-np.expm1(-diff * b)) / diff
    # Avoid inf*0 in extremely old transcripts.
    z = lo * b
    with np.errstate(over="ignore", invalid="ignore"):
        out = np.exp(-z) * (1.0 + correction)
    return np.where(z > 745, 0.0, out)


def _ramp_moments(a, intervals, time_per_bp):
    """Cumulative moments of a uniform-base transcribed-target fraction."""
    intervals = np.asarray(intervals, dtype=float)
    lengths = intervals[:, 1] - intervals[:, 0]
    if np.any(lengths <= 0) or np.any(intervals[:, 0] < 0):
        raise ValueError("Target intervals require 0 <= start < end")
    weights = lengths / lengths.sum()
    f0 = np.zeros_like(a)
    f1 = np.zeros_like(a)
    for (start_bp, end_bp), weight in zip(intervals, weights):
        start, end = time_per_bp * start_bp, time_per_bp * end_bp
        width = end - start
        u = np.clip(a - start, 0, width)
        tail = np.maximum(a - end, 0)
        f0 += weight * (u * u / (2 * width) + tail)
        f1 += weight * ((start * u * u / 2 + u ** 3 / 3) / width
                        + tail * (a + end) / 2)
    return f0, f1


def _ramp(a, intervals, time_per_bp):
    intervals = np.asarray(intervals, dtype=float)
    lengths = intervals[:, 1] - intervals[:, 0]
    out = np.zeros_like(a)
    for (start, end), weight in zip(intervals, lengths / lengths.sum()):
        out += weight * np.clip((a / time_per_bp - start) / (end - start), 0, 1)
    return out


def cumulative(age, geometry, T, splice_mean, decay_mean, mode="mapped"):
    """Return F0I,F1I,F0E,F1E, where Fn=integral_0^age x**n K(x) dx.

    The return arrays have the input age shape. Ages at or below zero have
    zero cumulative integrals. ``geometry`` is one gene's geometry dictionary.
    """
    _check(T, splice_mean, decay_mean, mode)
    age = np.asarray(age, dtype=float)
    a = np.maximum(age, 0)
    ks, kd = 1.0 / splice_mean, 1.0 / decay_mean
    if mode == "collapsed":
        f0i, f1i = _exp_moment(a, ks, 0), _exp_moment(a, ks, 1)
        f0e, f1e = _sumexp_moments(a, ks, kd)
        return f0i, f1i, f0e, f1e
    tp = T / geometry["transcript_end_bp"]
    ts = tp * geometry["splice_site_bp"]
    if not 0 <= ts <= T:
        raise ValueError("Splice site must precede or equal transcript completion")
    if max(e for _, e in geometry["intron_targets"]) > geometry["splice_site_bp"]:
        raise ValueError("Intron target must end before the target intron's splice site")
    if max(e for _, e in geometry["exon_targets"]) > geometry["transcript_end_bp"]:
        raise ValueError("Exon target must end before transcript completion")
    f0i, f1i = _ramp_moments(np.minimum(a, ts), geometry["intron_targets"], tp)
    b = np.maximum(a - ts, 0)
    j0, j1 = _exp_moment(b, ks, 0), _exp_moment(b, ks, 1)
    f0i += j0
    f1i += ts * j0 + j1
    f0e, f1e = _ramp_moments(np.minimum(a, T), geometry["exon_targets"], tp)
    b = np.maximum(a - T, 0)
    q = np.exp(-(T - ts) * ks)
    h0, h1 = _sumexp_moments(b, ks, kd)
    post0 = (1 - q) * _exp_moment(b, kd, 0) + q * h0
    post1 = (1 - q) * _exp_moment(b, kd, 1) + q * h1
    f0e += post0
    f1e += T * post0 + post1
    return f0i, f1i, f0e, f1e


def age_kernel(age, geometry, T, splice_mean, decay_mean, mode="mapped"):
    """Return intronic and exonic expected signal at transcript age(s)."""
    _check(T, splice_mean, decay_mean, mode)
    age = np.asarray(age, dtype=float)
    a = np.maximum(age, 0)
    ks, kd = 1 / splice_mean, 1 / decay_mean
    if mode == "collapsed":
        return (np.where(age >= 0, np.exp(-ks * a), 0),
                np.where(age >= 0, _sumexp_survival(a, ks, kd), 0))
    tp = T / geometry["transcript_end_bp"]
    ts = tp * geometry["splice_site_bp"]
    ki = _ramp(a, geometry["intron_targets"], tp) * np.exp(-ks * np.maximum(a - ts, 0))
    b = np.maximum(a - T, 0)
    q = np.exp(-(T - ts) * ks)
    survival = (1 - q) * np.exp(-kd * b) + q * _sumexp_survival(b, ks, kd)
    ke = _ramp(a, geometry["exon_targets"], tp) * survival
    return np.where(age >= 0, ki, 0), np.where(age >= 0, ke, 0)


def integrated_ramp(age, geometry, T, splice_mean, decay_mean, mode="mapped"):
    """Return H2(age)=integral_0^age (age-x) K(x) dx for both kernels."""
    a = np.maximum(np.asarray(age, dtype=float), 0)
    f0i, f1i, f0e, f1e = cumulative(a, geometry, T, splice_mean, decay_mean, mode)
    return a * f0i - f1i, a * f0e - f1e


def convolve_hats(times, knots, geometry, T, splice_mean, decay_mean, mode="mapped"):
    """Convolve kernels with piecewise-linear interpolation basis on knots.

    Returns two matrices of shape (len(times),len(knots)). Multiplying either
    by a vector of initiation values at the knots gives its exact convolution
    with the corresponding kernel, assuming initiation is zero outside the
    knots. Endpoint basis functions are the appropriate half hats and may
    jump from/to zero at the endpoints.
    """
    times, knots = np.asarray(times, dtype=float), np.asarray(knots, dtype=float)
    if times.ndim != 1 or knots.ndim != 1 or len(knots) < 2 or np.any(np.diff(knots) <= 0):
        raise ValueError("Require 1-D times and at least two strictly increasing knots")
    age = np.maximum(times[:, None] - knots[None, :], 0)
    f0i, f1i, f0e, f1e = cumulative(age, geometry, T, splice_mean, decay_mean, mode)
    hi, he = age * f0i - f1i, age * f0e - f1e
    widths = np.diff(knots)
    def build(f0, h):
        mat = np.empty_like(h)
        slopes = np.diff(h, axis=1) / widths
        mat[:, 0] = f0[:, 0] + slopes[:, 0]
        mat[:, -1] = -slopes[:, -1] - f0[:, -1]
        mat[:, 1:-1] = slopes[:, 1:] - slopes[:, :-1]
        # Exact entries are nonnegative; cancellation can produce tiny errors
        # in the remote tail when combining cumulative linear asymptotes.
        return np.maximum(mat, 0)
    return build(f0i, hi), build(f0e, he)
