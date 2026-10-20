"""Assemble reviewed peak-only tables and report content without changing sources."""
from pathlib import Path
import json, shutil
import pandas as pd
R=Path(__file__).resolve().parent
A=R.parent/'figure7_ver05_analysis'
G=A/'gap_peaks'; P=A/'pairrule'

def assemble():
    g=pd.read_csv(G/'all_pairwise_orders.csv'); p=pd.read_csv(P/'pairwise_peak_screen.csv')
    rows=[]
    for z in g.to_dict('records'):
        rows.append(dict(stage=z['stage'],embryo=z['embryo'],gene_a=z['gene_a'],gene_b=z['gene_b'],event_a=z['event_a'],event_b=z['event_b'],cohort_a=z['cohort_a'],cohort_b=z['cohort_b'],matched_temporal_channels=z['matched_temporal_channels'],temporal_difference_b_minus_a=z['temporal_difference_b_minus_a'],temporal_ci95_low=z['temporal_ci95_low'],temporal_ci95_high=z['temporal_ci95_high'],valid_bootstrap_fraction=z['valid_bootstrap_fraction'],position_a=z['position_a'],position_b=z['position_b'],span_a_lo=z['lobe_lo_a'],span_a_hi=z['lobe_hi_a'],span_b_lo=z['lobe_lo_b'],span_b_hi=z['lobe_hi_b'],centre_based_outcome=z['centre_based_outcome'],span_based_outcome=z['robust_outcome'],eligible=z['eligible'],strong_claim_eligible=z.get('strong_claim_eligible',False),availability_a=z['temporal_status_a'],availability_b=z['temporal_status_b'],disjoint_temporal_windows=not z['overlapping_temporal_search_windows'],flags=z['flags'],exclusion_reasons=z['exclusion_reasons'],source_suite='target_related'))
    def nonempty_flags(z):
        return ';'.join(str(z[k]) for k in ['support_role','flags_a','flags_b','qualification_flags','unresolved_reasons'] if k in z and pd.notna(z[k]) and str(z[k]).strip())
    for z in p.to_dict('records'):
        if z['gene_a']==z['gene_b']:continue
        rows.append(dict(stage=z['stage'],embryo=z['embryo'],gene_a=z['gene_a'],gene_b=z['gene_b'],event_a=z['event_a'],event_b=z['event_b'],cohort_a='shared_eve_run_odd',cohort_b='shared_eve_run_odd',matched_temporal_channels=True,temporal_difference_b_minus_a=z['temporal_delta_b_minus_a'],temporal_ci95_low=z['temporal_ci95_low'],temporal_ci95_high=z['temporal_ci95_high'],valid_bootstrap_fraction=z['temporal_valid_bootstrap_fraction'],position_a=z['position_a'],position_b=z['position_b'],span_a_lo=z['lobe_a_lo'],span_a_hi=z['lobe_a_hi'],span_b_lo=z['lobe_b_lo'],span_b_hi=z['lobe_b_hi'],centre_based_outcome=z['centre_based_outcome'],span_based_outcome=z['robust_outcome'],eligible=z.get('eligible',z['centre_based_outcome'] in ['preserved','reversed']),strong_claim_eligible=z.get('strong_claim_eligible',False),availability_a=z['temporal_availability_a'],availability_b=z['temporal_availability_b'],disjoint_temporal_windows=z['window_separated'],flags=nonempty_flags(z),exclusion_reasons=z['exclusion_reasons'],source_suite='shared_eve_run_odd'))
    out=pd.DataFrame(rows).sort_values(['stage','gene_a','gene_b','event_a','event_b','embryo'])
    out.to_csv(R/'Peak_Order_Comparisons_ver05.csv',index=False)
    pairs=sorted(set(tuple(sorted([x.gene_a,x.gene_b])) for x in out.itertuples()))
    assert len(pairs)==21, pairs
    # A compact coverage inventory is not a set of independent hypothesis tests.
    coverage=[]
    for a,b in pairs:
        q=out[((out.gene_a==a)&(out.gene_b==b))|((out.gene_a==b)&(out.gene_b==a))]
        coverage.append(dict(gene_1=a,gene_2=b,audit_rows=len(q),matched_temporal_channels=bool(q.matched_temporal_channels.all()),eligible_rows=int(q.eligible.sum()),note='Relationships are event- and stage-specific; rows share genes and embryos.'))
    pd.DataFrame(coverage).to_csv(R/'Peak_Order_Gene_Pair_Coverage_ver05.csv',index=False)
    shutil.copyfile(P/'posterior_mask.csv',R/'Peak_Order_Posterior_Region_ver05.csv')


if __name__ == '__main__':
    assemble()
