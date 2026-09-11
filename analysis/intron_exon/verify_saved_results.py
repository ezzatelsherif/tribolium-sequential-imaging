"""Check inputs, stored fits and summaries without running a long optimization.

This is an internal consistency check, not a test of biological validity.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from pulse_fit import ROOT, GENES, PulseFit, load
from sweep import load_geometry, scenario, summarize


def read_fasta(path):
    records = {}
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            records[name] = ''
        elif line.strip():
            records[name] += line.strip()
    return records


def main(skip_manifest=False):
    if not skip_manifest:
        manifest = json.loads((ROOT / 'manifest.json').read_text())
        for name, expected in manifest['sha256'].items():
            actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            if actual != expected:
                raise AssertionError(f'File differs from this release: {name}. '
                                     'After intentional reruns use --skip-manifest.')
        print(f"Verified {len(manifest['sha256'])} release file hashes.")

    targets = read_fasta(ROOT / 'data' / 'target_sequences.fasta')
    spans = read_fasta(ROOT / 'data' / 'reference_spans.fasta')
    geometries = load_geometry()
    for gene, g in geometries.items():
        span = spans[gene]
        if len(span) != g['transcript_end_bp']:
            raise AssertionError(f'{gene} reference span length differs from geometry.')
        for channel, field in (('intron', 'intron_targets'), ('exon', 'exon_targets')):
            sequence = ''.join(span[a:b] for a,b in g[field])
            tail = g['excluded_query_tail'] if channel == 'exon' else ''
            if sequence + tail != targets[f'{gene}_{channel}']:
                raise AssertionError(f'{gene} {channel} target does not match its reference intervals.')
    print('Six target sequences match reference intervals with the specified mlpt terminal-A exclusion.')

    path = ROOT / 'empirical' / 'original_stage_lag_analysis.py'
    spec = importlib.util.spec_from_file_location('empirical_offsets', path)
    empirical = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(empirical)
    estimates = pd.read_csv(ROOT / 'empirical' / 'stage_lag_bootstrap.csv')
    rows = []
    for gene in GENES:
        pairs, records = empirical.load_gene(gene.upper())
        rows.extend(records)
        mean = np.array([a.mean(0) for a in pairs])
        for j, window in enumerate(empirical.WINDOWS[gene.upper()]):
            values = empirical.crossings(mean, window)
            for metric, value in zip(('rising', 'falling', 'falling_minus_rising'), values):
                saved = estimates[(estimates.gene == gene.upper()) &
                    (estimates.pulse == j+1) & (estimates.metric == metric)].estimate.item()
                np.testing.assert_allclose(saved, value, atol=1e-10, rtol=1e-10)
    observed = pd.read_csv(ROOT / 'data' / 'background_corrected_pairs.csv')
    actual = pd.DataFrame(rows)
    actual['stage'] = actual['stage'].astype(float)
    pd.testing.assert_frame_equal(actual, observed, check_dtype=False, atol=1e-12, rtol=1e-12)
    counts = {g: sum(map(len, stages)) for g, stages in load().items()}
    if counts != {'Kr': 116, 'mlpt': 118, 'svb': 116}:
        raise AssertionError(f'Unexpected complete paired-slot counts: {counts}')
    print(f'Workbook corrections match all 360 source slots; complete pairs: {counts}.')

    max_prediction_error = 0.
    pair_arrays = load()
    for gene in GENES:
        means = np.array([p.mean(0) for p in pair_arrays[gene]])
        y = means / means.max(0)
        for shape in ('cosine', 'gaussian'):
            fits = json.loads((ROOT / 'fitresults' / f'fit_{gene}_{shape}.json').read_text())
            for mode, fit in fits.items():
                model = PulseFit(gene, mode=mode, shape=shape, dt=fit['dt'])
                prediction = model.predict(fit)
                np.testing.assert_allclose(prediction, fit['pred'], atol=1e-9, rtol=1e-9)
                np.testing.assert_allclose(np.sum((prediction-y)**2), fit['cost'], atol=1e-9, rtol=1e-9)
                max_prediction_error = max(max_prediction_error,
                    float(np.max(np.abs(prediction-np.asarray(fit['pred'])))))
    print(f'All 12 stored fits reconstructed; maximum prediction difference {max_prediction_error:.3g}.')

    sweep = pd.read_csv(ROOT / 'sweep_results' / 'sweep_results.csv')
    if len(sweep) != 1944 or not sweep['valid'].all():
        raise AssertionError('Expected 1,944 valid mapped sweep scenarios.')
    recalculated = summarize(sweep, ['gene'])
    recorded = pd.read_csv(ROOT / 'sweep_results' / 'summary_by_gene.csv')
    pd.testing.assert_frame_equal(recalculated, recorded, check_dtype=False, atol=1e-10, rtol=1e-10)
    for gene in GENES:
        row = sweep[sweep.gene == gene].iloc[0]
        new = scenario(gene, geometries[gene], row['T'], row.processing_mean,
                       row.decay_mean, row.peak_fraction)
        for key in ('rise_50_offset', 'fall_50_offset'):
            np.testing.assert_allclose(new[key], row[key], atol=1e-9, rtol=1e-9)
    print('Sweep summary and three independently rerun scenarios agree.')

    for gene in GENES:
        result = json.loads((ROOT / 'assessment' / f'assessment_{gene}.json').read_text())
        null = pd.read_csv(ROOT / 'assessment' / f'null_{gene}.csv')
        count = int((null.SSE >= result['observed_SSE']).sum())
        np.testing.assert_allclose((1+count)/(len(null)+1), result['conditional_lack_of_fit_p'])
    print('Conditional lack-of-fit summaries match all recorded null draws.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-manifest', action='store_true')
    main(**vars(parser.parse_args()))
