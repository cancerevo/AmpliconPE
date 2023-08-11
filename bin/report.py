#!/usr/bin/env python3
"""

Questions-Evidence (E) pairs:
------------------

1. Empirical CDF plots

2. Is the library Amplicon?
    E: High unaligned rate -> contamination
    E: Flat Cap -> wrong reference

3. Did any samples have poor PCR?
    E: fwd-rev behind fwd/ref - ref

4. Did any samples have poor index matching?
    E: Index mismatch %
        
6. Did any samples have poor coverage?
    E: Read coverage..above 50%
    
"""
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

threshold = 0.5

info = pd.read_csv("consolidated_info.csv", index_col=["alignment", "sample"])[
    "max_value"
]
index_mismatches = info.pop("Index Mismatches")

gb = info.groupby(level="sample")
max_scores = gb.median().astype(int)
max_final_score = max_scores['core-ref']
if not np.allclose(gb.std(), np.zeros(len(gb))):
    raise ValueError("Samples have different maximum alignment values.")

align_pairs = max_scores.index.values


pileups = pd.read_csv("consolidated_pileups.csv", index_col=['sample', 'barcode']).eval(
    'total alignment = score*read/@max_final_score')

pileups_marginal = pileups.groupby(level=[0,1])['total alignment', 'reads']).sum()
pileups_marginal['mean alignment'] = pileups_marginal.pop('total alignment') / pileups_marginal['reads'] 


scores = pd.read_csv(
    "consolidated_scores.csv", index_col=["sample"] + align_pairs, dtype=np.uint64
)["reads"]
base_scores = pd.DataFrame(
    {name: scores.groupby(level=name).sum() for name in align_pairs}
)
normed_scores = base_scores.div(max_scores, axis=1)
normed_scores['pileups'] = np.histogram(
    pileups_marginal['mean alignment']
    weights=pileups_marginal['reads'],
    bins=len(normed_scores),
    range=(0,1)
    )[0]
print(normed_scores.describe())
normed_scores.columns.name = "alignment"
longform = normed_scores.stack()
longform.name = "reads"


sns.ecdfplot(
    data=longform.reset_index(),
    x="score",
    weights=longform.name,
    hue="alignment",
    stat="count",
    hue_order=normed_scores.columns,
    ax=plt.gca(),
)
plt.savefig("alignments.pdf")

def length_std(bc):
    return pd.Series(bc).str.len().std()

from scipy.stats import entropy

pile_stats = pile_gb_score.agg(dict(barcode=length_stdev, reads=entropy)
pile_stats['fwd-rev alignment'] = scores.groupby('core-ref')['fwd-rev'].mean() / max_scores['fwd-rev']
pile_stats.columns.name = 'stat'
pile_stats /= p.stats.max()/100

line_data = pile_stats.stack()
line_data.name = '% of Max'

ax = plt.gca()
sns.lineplot(
    x='score',
    y=line_data.name,
    hue='stat',
    ax=ax)

ax.vline(threshold, ls=':')
plt.savefig('alignment_profile.pdf')

unaligned = pileups.query('score <= @threshold').groupby(['sample', 'barcode'])['total alignment', 'reads'].sum()
unaligned['mean alignment'] = unaligned.pop('total alignment') / unaligned['reads']
unaligned.sort_values(by='reads', ascending=False).to_csv('unaligned.csv')

aligned = pileups.query('score > @threshold').groupby(['sample', 'barcode']['total alignment', 'reads'].sum()
aligned['mean alignment'] = aligned.pop('total alignment') / aligned['reads']
aligned.to_csv('aligned_pileups.csv')

outcomes = ['Adapter Index Mismatches', 'Unaligned', 'Aligned']
sample_stats = pd.concat([
    index_mismatches, 
    unaligned.groupy(level='sample')['reads'].sum(),
    aligned.groupby(level='sample')['reads'].sum()],
    ,keys=outcomes , names=['outcome'])


sns.barplot(
    data=sample_stats,
    x='sample',
    y='reads',
    hue='outcome',
    hue_order=outcomes,
    ax=plt.gca())

plt.savefig('outcomes.pdf')


