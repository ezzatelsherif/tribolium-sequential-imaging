"""Repeat pair-rule Figure 7 measurements using Theo's original smoothing.

The smoother is applied separately to each of 20,000 jointly resampled
stage-mean curves. No bootstrap resamples spatial pixels. Raw values persist.
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / 'all_gene_direct_analysis'
PREVIOUS = ROOT / 'figure7_ver03_analysis' / 'pairrule'
sys.path.insert(0, str(HERE.parent))
import smoothing

LABELS = [f'{c}.{s}' for c in range(1, 9) for s in range(1, 4)]
CYCLES = [2, 3, 4, 5]


def summary(v):
    obs = float(v[0])
    valid = v[1:][np.isfinite(v[1:])]
    q = np.quantile(valid, [.025, .975]) if len(valid) and np.isfinite(obs) else [np.nan, np.nan]
    return dict(estimate=obs, ci95_low=q[0], ci95_high=q[1],
                valid_bootstrap_n=len(valid), valid_bootstrap_fraction=len(valid) / (len(v) - 1))


def temporal_events(x, y, cycle, half_bin=True):
    """Common window for all genes; continuous-grid peaks and flanking troughs."""
    n, m = y.shape
    rows = np.arange(n)
    lo = 3 * (cycle - 1) - (0.5 if half_bin else 0.)
    hi = 3 * (cycle - 1) + 2 + (0.5 if half_bin else 0.)
    ids = np.flatnonzero((x >= lo) & (x <= hi))
    candidate = ids[np.argmax(y[:, ids], axis=1)]
    local_max = np.zeros_like(y, dtype=bool)
    local_max[:, 1:-1] = ((y[:, 1:-1] >= y[:, :-2]) & (y[:, 1:-1] >= y[:, 2:])
                          & ((y[:, 1:-1] > y[:, :-2]) | (y[:, 1:-1] > y[:, 2:])))
    eligible = local_max[:, ids]
    peak = ids[np.argmax(np.where(eligible, y[:, ids], -np.inf), axis=1)]
    has_maximum = eligible.any(1)
    peak[~has_maximum] = candidate[~has_maximum]
    height = y[rows, peak]
    baseline = np.min(y, axis=1)
    level = baseline + .5 * (height - baseline)
    resolved = has_maximum & (height > baseline)
    local_min = np.zeros_like(y, dtype=bool)
    local_min[:, 1:-1] = ((y[:, 1:-1] <= y[:, :-2]) & (y[:, 1:-1] <= y[:, 2:])
                          & ((y[:, 1:-1] < y[:, :-2]) | (y[:, 1:-1] < y[:, 2:])))
    xx = x[None, :]
    pkx = x[peak]
    left = local_min & (xx < pkx[:, None]) & (xx >= pkx[:, None] - 3)
    right = local_min & (xx > pkx[:, None]) & (xx <= pkx[:, None] + 3)
    li = m - 1 - left[:, ::-1].argmax(1)
    ri = right.argmax(1)
    above = y >= level[:, None]
    ix = np.arange(m - 1)[None, :]
    out = {'peak': pkx.copy()}
    out['peak'][~resolved] = np.nan
    for side in ['rise', 'fall']:
        crosses = (~above[:, :-1]) & above[:, 1:] if side == 'rise' else above[:, :-1] & (~above[:, 1:])
        region = (ix >= li[:, None]) & (ix < peak[:, None]) if side == 'rise' else (ix >= peak[:, None]) & (ix < ri[:, None])
        crosses &= region
        ii = crosses.argmax(1) if side == 'rise' else (m - 2) - crosses[:, ::-1].argmax(1)
        with np.errstate(divide='ignore', invalid='ignore'):
            val = x[ii] + (level - y[rows, ii]) * (x[ii + 1] - x[ii]) / (y[rows, ii + 1] - y[rows, ii])
        has_trough = left.any(1) if side == 'rise' else right.any(1)
        val[~resolved | ~has_trough | ~crosses.any(1)] = np.nan
        out[side] = val
    meta = dict(peak_candidate=float(pkx[0]), baseline=float(baseline[0]),
                peak_height=float(height[0]), halfheight_level=float(level[0]),
                normalized_halfheight_level=float(level[0]), normalized_peak_height=float(height[0]),
                peak_resolved=bool(resolved[0]), peak_clipped=bool(not resolved[0]),
                highest_window_value_at_selected_peak=bool(candidate[0] == peak[0]),
                left_trough_index=float(x[li[0]]) if left.any(1)[0] else np.nan,
                right_trough_index=float(x[ri[0]]) if right.any(1)[0] else np.nan,
                window_lo_index=lo, window_hi_index=hi,
                peak_window_definition='same_three_sample_bins' if half_bin else 'strict_sample_coordinates')
    return out, meta


def cut(x, y, lo, hi):
    inside = (x > lo) & (x < hi)
    return np.r_[lo, x[inside], hi], np.r_[np.interp(lo, x, y), y[inside], np.interp(hi, x, y)]


def crossing(x, y, level, side):
    above = np.flatnonzero(y >= level)
    if not len(above):
        return np.nan
    k = int(above[0] if side == 'rise' else above[-1])
    if (side == 'rise' and k == 0) or (side == 'fall' and k == len(y) - 1):
        return np.nan
    i, j = (k - 1, k) if side == 'rise' else (k, k + 1)
    return float(x[i] + (level - y[i]) * (x[j] - x[i]) / (y[j] - y[i]))


def connected_crossing(x, y, peak, level, side):
    if side == 'rise':
        ids = np.flatnonzero((x[:-1] < peak) & (y[:-1] < level) & (y[1:] >= level))
        i = int(ids[-1]) if len(ids) else None
    else:
        ids = np.flatnonzero((x[1:] > peak) & (y[:-1] >= level) & (y[1:] < level))
        i = int(ids[0]) if len(ids) else None
    if i is None:
        return np.nan
    return float(x[i] + (level - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))


def temporal_analysis():
    z = np.load(SOURCE / 'temporal_bootstrap_curves.npz')
    stage_summary = pd.read_csv(SOURCE / 'temporal_stage_summary.csv')
    profile_rows, raw_rows, peak_rows, feature_rows, strict_rows, draws = [], [], [], [], [], {}
    for gene in ['eve', 'run', 'odd']:
        raw = z[f'paired_pairrule_{gene}']
        x, obs = smoothing.temporal(raw[0])
        smoothed = np.empty((len(raw), len(x)), dtype=np.float64)
        metrics = {(cycle, strict, event): np.empty(len(raw)) for cycle in CYCLES
                   for strict in [False, True] for event in ['peak', 'rise', 'fall']}
        metadata = {}
        for start in range(0, len(raw), 256):
            end = min(len(raw), start + 256)
            _, yy = smoothing.temporal(raw[start:end])
            smoothed[start:end] = yy
            for cycle in CYCLES:
                for strict in [False, True]:
                    events, meta = temporal_events(x, yy, cycle, half_bin=not strict)
                    if start == 0:
                        metadata[(cycle, strict)] = meta
                    for event, vals in events.items():
                        metrics[(cycle, strict, event)][start:end] = vals
        q = np.quantile(smoothed[1:], [.025, .975], axis=0)
        for xi, yi, qlo, qhi in zip(x, obs, q[0], q[1]):
            profile_rows.append(dict(gene=gene, cohort='pairrule', stage='', substage_index=xi,
                                     raw_mean=np.nan, normalized_mean=yi, normalized_ci95_low=qlo,
                                     normalized_ci95_high=qhi, embryo_n=np.nan,
                                     normalization='whole_raw_series_min_max_before_original_smoothing',
                                     is_original_observation=False))
        raw_norm = smoothing.normalize_range(raw[0])
        for i in range(24):
            rr = stage_summary.query('cohort == "pairrule" and gene == @gene and substage_index == @i').iloc[0]
            raw_rows.append(dict(gene=gene, cohort='pairrule', stage=LABELS[i], substage_index=i,
                                 raw_mean=raw[0, i], normalized_raw_mean=raw_norm[i],
                                 raw_sem=rr['sem'], normalized_raw_sem=rr['sem']/np.ptp(raw[0]), embryo_n=int(rr['n'])))
        for cycle in CYCLES:
            for strict in [False, True]:
                meta = metadata[(cycle, strict)]
                for event in ['peak', 'rise', 'fall']:
                    vals = metrics[(cycle, strict, event)]
                    row = dict(gene=gene, cohort='pairrule', cycle=cycle, event=event,
                               normalization='smoothed_full_series_min_to_cycle_max', threshold_fraction=.5,
                               unit='ordinal_substage', smoothing='original_PCHIP800_SG81_poly3', **meta, **summary(vals))
                    (strict_rows if strict else feature_rows).append(row)
                    if not strict:
                        draws[f'{gene}_cycle{cycle}_{event}'] = vals
            meta = metadata[(cycle, False)]
            lo, hi = meta['left_trough_index'], meta['right_trough_index']
            if np.isfinite(lo) and np.isfinite(hi):
                ids = np.flatnonzero((x >= lo) & (x <= hi))
                peaks, _ = find_peaks(obs[ids], prominence=.15*np.ptp(obs[ids]))
                for k in ids[peaks]:
                    peak_rows.append(dict(domain='time', stage='', embryo='', gene=gene, cycle=cycle,
                                          position=x[k], height=obs[k], selected=bool(np.isclose(x[k], meta['peak_candidate'])),
                                          prominence_rule='15_percent_episode_range_after_original_smoothing'))
        del smoothed
    pd.DataFrame(profile_rows).to_csv(HERE / 'temporal_profiles.csv', index=False)
    pd.DataFrame(raw_rows).to_csv(HERE / 'temporal_raw_observations.csv', index=False)
    pd.DataFrame(feature_rows).to_csv(HERE / 'temporal_features.csv', index=False)
    pd.DataFrame(strict_rows).to_csv(HERE / 'temporal_features_strict_window_sensitivity.csv', index=False)
    np.savez_compressed(HERE / 'temporal_feature_bootstraps.npz', **draws)
    contrasts = []
    for cycle in CYCLES:
        for event in ['rise', 'peak', 'fall']:
            for a, b in [('run', 'eve'), ('odd', 'run'), ('odd', 'eve')]:
                d = draws[f'{a}_cycle{cycle}_{event}'] - draws[f'{b}_cycle{cycle}_{event}']
                valid = d[1:][np.isfinite(d[1:])]
                contrasts.append(dict(cycle=cycle, event=event, contrast=f'{a}_minus_{b}', **summary(d),
                                      bootstrap_fraction_positive=float(np.mean(valid > 0)) if len(valid) else np.nan,
                                      bootstrap_fraction_equal=float(np.mean(valid == 0)) if len(valid) else np.nan))
    pd.DataFrame(contrasts).to_csv(HERE / 'temporal_paired_contrasts.csv', index=False)
    return peak_rows


def spatial_analysis(peak_rows):
    raw = pd.read_csv(SOURCE / 'spatial_raw_profiles.csv')
    raw = raw[raw.gene.isin(['eve', 'runt', 'odd'])].copy()
    raw['source_gene'] = raw.gene
    raw['gene'] = raw.gene.replace({'runt': 'run'})
    raw.to_csv(HERE / 'spatial_raw_observations.csv', index=False)
    profiles = []
    lookup = {}
    for key, p in raw.groupby(['stage', 'embryo', 'gene'], sort=False):
        stage, emb, gene = key
        x, y = smoothing.spatial(p.x_ap.to_numpy(), p.raw_intensity.to_numpy(), stage, emb, gene)
        out = p.copy()
        out['raw_normalized_intensity'] = out['normalized_intensity']
        out['normalized_intensity'] = y
        out['smoothed_normalized_intensity'] = y
        kw = smoothing.spatial_parameters(stage, emb, gene)
        out['smoothing_frac'] = kw['frac']
        out['smoothing_polyorder'] = kw['polyorder']
        out['smoothing_window_length'] = smoothing.window_length(len(x), **kw)
        profiles.append(out)
        lookup[key] = out
    pd.concat(profiles).to_csv(HERE / 'spatial_profiles.csv', index=False)
    windows = pd.read_csv(SOURCE / 'spatial_common_window_metrics.csv')
    features = []
    for row in windows.itertuples(index=False):
        gene = 'run' if row.gene == 'runt' else row.gene
        p = lookup[(row.stage, row.embryo, gene)]
        x, y = p.x_ap.to_numpy(), p.normalized_intensity.to_numpy()
        lo, hi = float(row.window_lo), float(row.window_hi)
        xx, yy = cut(x, y, lo, hi)
        ids = np.flatnonzero((x >= lo) & (x <= hi))
        raw_candidate = int(ids[np.argmax(y[ids])])
        all_local_max, _ = find_peaks(y)
        eligible_max = all_local_max[(x[all_local_max] > lo) & (x[all_local_max] < hi)]
        kp = int(eligible_max[np.argmax(y[eligible_max])]) if len(eligible_max) else raw_candidate
        height, baseline = float(y[kp]), float(y.min())
        level = baseline + .5*(height-baseline)
        pp, _ = find_peaks(yy, prominence=.15*np.ptp(yy))
        prom = xx[pp]
        split = bool(len(prom) > 1 and np.ptp(prom) >= .04)
        mapped = bool(row.homologous_event_previously_identified)
        edge = not len(eligible_max)
        developing = bool(row.posterior_developing_cycle)
        resolved = mapped and not edge
        for k in pp:
            peak_rows.append(dict(domain='space', stage=row.stage, embryo=row.embryo, gene=gene, cycle=row.cycle,
                                  position=xx[k], height=yy[k], selected=bool(np.isclose(xx[k], x[kp])),
                                  prominence_rule='15_percent_episode_range_after_original_smoothing'))
        vals = dict(peak=float(x[kp]) if resolved else np.nan,
                    rise=crossing(xx, yy, level, 'rise') if resolved else np.nan,
                    fall=crossing(xx, yy, level, 'fall') if resolved else np.nan)
        flags = []
        if not mapped: flags.append('unmapped')
        if edge: flags.append('peak_clipped')
        if developing: flags.append('posterior_developing_cycle')
        if not row.complete_coarse_window: flags.append('incomplete_window')
        if split: flags.append('multiple_prominent_peaks')
        for side, name in [('rise', 'anterior'), ('fall', 'posterior')]:
            if not np.isfinite(vals[side]): flags.append(f'halfheight_{name}_unresolved')
        meta = dict(stage=row.stage, embryo=int(row.embryo), gene=gene, cycle=int(row.cycle), window_lo=lo, window_hi=hi,
                    normalization='smoothed_full_profile_min_to_cycle_max', threshold_fraction=.5, unit='AP_fraction',
                    baseline=baseline, peak_candidate=float(x[kp]), peak_height=height, halfheight_level=level,
                    highest_window_value_at_selected_peak=bool(raw_candidate == kp),
                    homologous_event_previously_identified=mapped, peak_clipped=edge, multiple_prominent_peaks=split,
                    posterior_developing_cycle=developing, complete_window=bool(row.complete_coarse_window),
                    primary_order_eligible=resolved and not developing and bool(row.complete_coarse_window),
                    smoothing_frac=float(p.iloc[0].smoothing_frac), smoothing_polyorder=int(p.iloc[0].smoothing_polyorder),
                    smoothing_window_length=int(p.iloc[0].smoothing_window_length),
                    flags=';'.join(flags), prominent_peak_positions=json.dumps(prom.tolist()))
        features.extend(dict(**meta, event=event, estimate=val) for event, val in vals.items())
    common = pd.DataFrame(features)
    common.to_csv(HERE / 'spatial_features.csv', index=False)
    pd.DataFrame(peak_rows).to_csv(HERE / 'prominent_peaks.csv', index=False)
    gene_episode_boundaries(common, lookup)


def gene_episode_boundaries(common, lookup):
    own = pd.read_csv(SOURCE / 'spatial_episode_metrics.csv')
    own = own[own.group == 'pair_rule'].copy()
    own['gene'] = own.gene.replace({'runt': 'run'})
    features, profiles, audit = [], [], []
    for row in common[common.event == 'peak'].to_dict('records'):
        stage, emb, gene, cycle = row['stage'], row['embryo'], row['gene'], row['cycle']
        match = own[(own.stage == stage) & (own.embryo == emb) & (own.gene == gene) & (own.cycle == cycle)]
        p = lookup[(stage, emb, gene)]
        x, y = p.x_ap.to_numpy(), p.normalized_intensity.to_numpy()
        q, peak = row['halfheight_level'], row['peak_candidate']
        base = dict(stage=stage, embryo=int(emb), gene=gene, cycle=int(cycle),
                    common_window_lo=row['window_lo'], common_window_hi=row['window_hi'],
                    baseline=row['baseline'], peak_candidate=peak, peak_height=row['peak_height'], halfheight_level=q,
                    common_peak_estimate=row['estimate'], common_peak_clipped=row['peak_clipped'],
                    homologous_event_previously_identified=row['homologous_event_previously_identified'],
                    posterior_developing_cycle=row['posterior_developing_cycle'], multiple_prominent_peaks=row['multiple_prominent_peaks'],
                    normalization='smoothed_full_profile_min_to_common_window_cycle_max', threshold_fraction=.5, unit='AP_fraction')
        if len(match):
            z = match.iloc[0]
            lo, hi = float(z.window_lo), float(z.window_hi)
            xx, yy = cut(x, y, lo, hi)
            inside = lo < peak < hi
            resolved = bool(inside and np.isfinite(row['estimate']))
            for xi, yi in zip(xx, yy):
                profiles.append(dict(stage=stage, embryo=int(emb), gene=gene, cycle=int(cycle), x_ap=xi,
                                     normalized_intensity=yi, raw_normalized_intensity=np.interp(xi, x, p.raw_normalized_intensity),
                                     original_sample=bool(np.any(x == xi)), gene_window_lo=lo, gene_window_hi=hi,
                                     common_window_lo=row['window_lo'], common_window_hi=row['window_hi']))
            ar = dict(**base, gene_window_lo=lo, gene_window_hi=hi, peak_inside_gene_window=inside,
                      gene_window_complete=bool(z.complete_coarse_window))
            for definition in ['outer_episode_envelope', 'connected_selected_peak']:
                vals = {'peak': row['estimate'] if resolved else np.nan}
                for side in ['rise', 'fall']:
                    vals[side] = (crossing(xx, yy, q, side) if definition == 'outer_episode_envelope'
                                  else connected_crossing(xx, yy, peak, q, side)) if resolved else np.nan
                    ar[f'{definition}_{side}'] = vals[side]
                flags = []
                if not inside: flags.append('common_peak_outside_gene_episode')
                if not np.isfinite(row['estimate']): flags.append('common_peak_unresolved')
                if not bool(z.complete_coarse_window): flags.append('incomplete_gene_window')
                if row['posterior_developing_cycle']: flags.append('posterior_developing_cycle')
                if row['multiple_prominent_peaks']: flags.append('multiple_prominent_peaks')
                for side, name in [('rise', 'anterior'), ('fall', 'posterior')]:
                    if not np.isfinite(vals[side]): flags.append(f'{name}_halfheight_unresolved')
                for event, val in vals.items():
                    features.append(dict(**base, gene_window_lo=lo, gene_window_hi=hi,
                                         gene_window_complete=bool(z.complete_coarse_window), peak_inside_gene_window=inside,
                                         boundary_definition=definition, event=event, estimate=val,
                                         valid_estimate=bool(np.isfinite(val)), flags=';'.join(flags)))
            audit.append(ar)
        else:
            for definition in ['outer_episode_envelope', 'connected_selected_peak']:
                for event in ['peak', 'rise', 'fall']:
                    features.append(dict(**base, gene_window_lo=np.nan, gene_window_hi=np.nan,
                                         gene_window_complete=False, peak_inside_gene_window=False,
                                         boundary_definition=definition, event=event, estimate=np.nan,
                                         valid_estimate=False, flags='no_original_gene_episode'))
    pd.DataFrame(features).to_csv(HERE / 'spatial_features_gene_episode.csv', index=False)
    pd.DataFrame(profiles).to_csv(HERE / 'spatial_profiles_gene_episode.csv', index=False)
    pd.DataFrame(audit).to_csv(HERE / 'spatial_boundary_window_audit.csv', index=False)


def annotate_neighbor_tails():
    """Keep the presence of a higher cropped neighboring tail explicit."""
    for name in ['temporal_features.csv', 'temporal_features_strict_window_sensitivity.csv', 'spatial_features.csv']:
        d = pd.read_csv(HERE / name)
        tail = (~d['highest_window_value_at_selected_peak']) & (~d['peak_clipped'])
        d['higher_neighbor_tail_at_window_edge'] = tail
        if 'flags' not in d:
            d['flags'] = ''
        d['flags'] = d['flags'].fillna('')
        for k in d.index[tail]:
            tags = [v for v in d.at[k, 'flags'].split(';') if v]
            if 'higher_neighbor_tail_at_window_edge' not in tags:
                tags.append('higher_neighbor_tail_at_window_edge')
            d.at[k, 'flags'] = ';'.join(tags)
        d.to_csv(HERE / name, index=False)
    common = pd.read_csv(HERE / 'spatial_features.csv')
    keys = ['stage', 'embryo', 'gene', 'cycle']
    flag = common[common.event == 'peak'][keys + ['higher_neighbor_tail_at_window_edge']]
    own = pd.read_csv(HERE / 'spatial_features_gene_episode.csv')
    if 'higher_neighbor_tail_at_window_edge' in own:
        own = own.drop(columns='higher_neighbor_tail_at_window_edge')
    own = own.merge(flag, on=keys, how='left')
    own['flags'] = own['flags'].fillna('')
    for k in own.index[own.higher_neighbor_tail_at_window_edge.fillna(False)]:
        tags = [v for v in own.at[k, 'flags'].split(';') if v]
        if 'higher_neighbor_tail_at_window_edge' not in tags:
            tags.append('higher_neighbor_tail_at_window_edge')
        own.at[k, 'flags'] = ';'.join(tags)
    own.to_csv(HERE / 'spatial_features_gene_episode.csv', index=False)


def compare_and_summarize():
    rows = []
    for domain, filename, keys in [
        ('time', 'temporal_features.csv', ['gene', 'cycle', 'event']),
        ('space', 'spatial_features_gene_episode.csv', ['stage', 'embryo', 'gene', 'cycle', 'event', 'boundary_definition'])]:
        old, new = pd.read_csv(PREVIOUS / filename), pd.read_csv(HERE / filename)
        merge = new.merge(old[keys + ['estimate']].rename(columns={'estimate':'raw_estimate'}), on=keys, how='left')
        merge['smoothed_estimate'] = merge['estimate']
        merge['smoothed_minus_raw'] = merge['estimate'] - merge['raw_estimate']
        merge['domain'] = domain
        rows.append(merge)
    pd.concat(rows, ignore_index=True).to_csv(HERE / 'raw_smoothed_feature_comparison.csv', index=False)
    space = pd.read_csv(HERE / 'spatial_features_gene_episode.csv')
    space = space[space.boundary_definition == 'outer_episode_envelope']
    order_rows = []
    for (stage, emb, cycle, event), p in space.groupby(['stage', 'embryo', 'cycle', 'event']):
        values = p.set_index('gene').estimate.to_dict()
        all_valid = all(np.isfinite(values.get(g, np.nan)) for g in ['eve','run','odd'])
        order_rows.append(dict(stage=stage, embryo=emb, cycle=cycle, event=event,
                               eve=values.get('eve',np.nan), run=values.get('run',np.nan), odd=values.get('odd',np.nan),
                               all_resolved=all_valid, eve_run_odd=bool(all_valid and values['eve'] < values['run'] < values['odd']),
                               flags=';'.join(sorted(set(';'.join(p['flags'].fillna('')).split(';'))))))
    pd.DataFrame(order_rows).to_csv(HERE / 'spatial_order_summary.csv', index=False)
    notes = '''Pair-rule analysis repeated with the original Theo notebook smoother

Temporal: whole raw stage-mean series min-max normalization, PCHIP to 800
equally spaced positions, Savitzky-Golay window81 / degree3, mode interp.
This entire procedure is repeated for each of 20,000 existing paired
within-stage biological bootstrap draws. Confidence intervals are conditional
on the assigned stages and this smoothing model; they do not estimate stage
assignment error or calibrated elapsed time. Dense positions do not add
observed time points. Raw means, SEM and sample sizes remain in a separate file.

Peaks are the highest genuine local maximum inside the candidate interval;
a higher tail from the neighboring cycle at the interval boundary does not
replace an interior maximum. If no local maximum exists, the unresolved
boundary candidate is retained separately and never plotted as a true peak.
The primary peak-search interval represents the same three sampled substage
bins as v03: c.1-0.5 through c.3+0.5. This continuous-grid interpretation is
identical for all genes. A strict c.1 through c.3 sensitivity is also exported.
The latter clips odd peaks shifted slightly beyond c.3 by smoothing, despite
having resolved continuous maxima. No expected gene order is imposed.
Half-heights use the minimum of the whole smoothed series plus half the
distance to the selected smoothed cycle maximum. Nearest local minima within
three ordinal substages bracket the episode. Missing or shallow troughs leave
crossings unresolved. Cycles2-4 are retained and cycle5 is added to match the
spatial companion; complete temporal profiles include all24 sampled stages.

Spatial: each original embryo and channel is min-max normalized then smoothed
on its full native sample grid with the exact stage/embryo/gene parameters in
../original_smoothing_config.json. No stages or embryos are averaged. The
native sample count determines the original SG window. Raw values remain
in spatial_raw_observations.csv and raw_* columns of spatial_profiles.csv.
The common EVE-defined cycle window still selects the peak. Original
gene-specific cycle windows search for the two boundaries, with identical
full-smoothed-profile-minimum-to-cycle-maximum half-height on both sides.
Both the outer envelope and connected-peak sensitivity are provided. Clipped,
unmapped, incomplete and developing episodes remain flagged, not imputed.
SG may overshoot the raw normalized range slightly; no post-smoothing clipping
or renormalization conceals this behavior. Spatial pixels are not bootstrapped.

The original WDL cropping and resampling are not reused because Figure7
measures actual AP coordinates and retains the independently specified cycle
windows. Original smoothing is reused; no WDL score or hand-entered peak order
is used. Results are smoothing-dependent estimates, alongside the raw audit.
'''
    (HERE / 'README.txt').write_text(notes)


if __name__ == '__main__':
    spatial_analysis(temporal_analysis())
    annotate_neighbor_tails()
    compare_and_summarize()
    print(pd.read_csv(HERE / 'temporal_features.csv')[['gene','cycle','event','estimate','ci95_low','ci95_high','valid_bootstrap_fraction']].to_string(index=False))
    print(pd.read_csv(HERE / 'spatial_order_summary.csv').query('stage == "Stage4_3" and cycle == 3').to_string(index=False))
