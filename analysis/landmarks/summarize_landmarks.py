#!/usr/bin/env python3
"""Create the inventory and pooled/stratified summaries used in the S1 Data workbook."""
from pathlib import Path
import json, os
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
CATEGORIES = ['preserved_all3', 'reversed_all3', 'spatial_consistent_temporal_unresolved', 'mixed', 'other_unresolved']
CATEGORY_LABELS = {
    'preserved_all3': 'Preserved order in all 3 embryos',
    'reversed_all3': 'Reversed order in all 3 embryos',
    'spatial_consistent_temporal_unresolved': 'Unresolved in time, separated in space',
    'mixed': 'Different spatial orders among embryos',
    'other_unresolved': 'Other unresolved comparisons',
    'excluded': 'Excluded',
}
GROUPS = ['pair-rule', 'gap-gap', 'gap-pair-rule']
EXPECTED = {
    ('all_landmarks','pair-rule'): [32,0,0,1,0],
    ('all_landmarks','gap-gap'): [28,3,5,4,1],
    ('all_landmarks','gap-pair-rule'): [117,13,9,15,5],
    ('peak','pair-rule'): [18,0,0,0,0],
    ('peak','gap-gap'): [12,1,2,2,1],
    ('peak','gap-pair-rule'): [60,5,4,8,3],
    ('half_height','pair-rule'): [14,0,0,1,0],
    ('half_height','gap-gap'): [16,2,3,2,0],
    ('half_height','gap-pair-rule'): [57,8,5,7,2],
}

def truth(x):
    return str(x).lower() in ['true','1','1.0']

def classify(row):
    if not truth(row.eligible):
        return 'excluded'
    if row.outcome in ('preserved_all3','reversed_all3','mixed'):
        return row.outcome
    if row.outcome == 'unresolved_temporal':
        values = np.array([row.ds_emb1,row.ds_emb2,row.ds_emb3],dtype=float)
        if np.isfinite(values).all() and ((values > 1e-12).all() or (values < -1e-12).all()):
            return 'spatial_consistent_temporal_unresolved'
    if row.outcome not in ('unresolved_temporal','unresolved_spatial'):
        raise ValueError(f'Unexpected eligible outcome: {row.outcome}')
    return 'other_unresolved'

def write(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_name(path.name+'.writing')
    with temporary.open('wb') as stream:
        payload=text.encode('utf-8');stream.write(payload);stream.flush();os.fsync(stream.fileno())
    os.replace(temporary,path)
    assert path.stat().st_size==len(payload)

def main():
    data=pd.read_csv(ROOT/'data/landmark_distances.csv')
    assert len(data)==859 and not data.comparison_id.duplicated().any()
    data['summary_panel']=np.where(data.eligible.map(truth),'all_landmarks','excluded')
    data['summary_category']=data.apply(classify,axis=1)
    data['summary_category_label']=data.summary_category.map(CATEGORY_LABELS)
    eligible=data[data.eligible.map(truth)]
    assert len(eligible)==233
    assert int((eligible.landmark=='peak').sum())==116
    assert int(eligible.landmark.isin(['rise','fall']).sum())==117
    assert int((data.summary_category=='excluded').sum())==626
    rows=[]
    for panel,scope in [('all_landmarks','pooled_landmarks'),('peak','peaks'),('half_height','half_heights')]:
        pool=eligible if panel=='all_landmarks' else eligible[eligible.landmark.eq('peak') if panel=='peak' else eligible.landmark.isin(['rise','fall'])]
        for group in GROUPS:
            subset=pool[pool.group==group]
            counts=subset.summary_category.value_counts().reindex(CATEGORIES,fill_value=0)
            assert counts.tolist()==EXPECTED[(panel,group)],(panel,group,counts.tolist())
            for category,count in counts.items():
                rows.append(dict(panel=panel,group=group,category=category,count=int(count),total=len(subset),percentage=100*count/len(subset),summary_scope=scope,figure_panel='B' if panel=='all_landmarks' else '',category_label=CATEGORY_LABELS[category]))
    summary=pd.DataFrame(rows)
    assert summary[summary.panel.eq('all_landmarks')].groupby('category')['count'].sum().reindex(CATEGORIES).tolist()==[177,16,14,20,6]
    write(ROOT/'outputs/S2_Data_Landmark_Comparisons.csv',data.to_csv(index=False))
    write(ROOT/'outputs/S3_Data_Landmark_Summary.csv',summary.to_csv(index=False))
    print(json.dumps(dict(candidates=len(data),eligible=len(eligible),excluded=626,peak_comparisons=116,half_height_comparisons=117,summary_rows=len(summary),plotted_pooled_rows=15,retained_stratified_rows=30,figure_panel='B'),indent=2))

if __name__=='__main__':
    main()
