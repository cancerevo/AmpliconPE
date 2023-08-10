#!/usr/bin/env python3
"""

Questions-Evidence (E) pairs:
------------------

1. Empirical CDF plots

2. Is the library Amplicon?
    E: High unaligned rate -> contamination
    E: Flat Cap -> wrong reference

3. Did any samples have poor PCR?
    E: Fwd-Rev behind Fwd/Ref - Ref

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
if not np.allclose(gb.std(), np.zeros(len(gb))):
    raise ValueError("Samples have different maximum alignment values.")

align_pairs = max_scores.index.values

scores = pd.read_csv(
    "consolidated_scores.csv", index_col=["sample"] + align_pairs, dtype=np.uint64
)["reads"]


pileups = pd.read_csv("consolidated_pileups.csv", index_col=['sample', 'barcode'])

pileups_gb = pileups.groupby(level=[0,1])

pileups_marginal = pd.DataFrame({
    'mean alignment' : pileups_gb.agg(lambda df: df.eval('score * reads').sum()),
    'reads': pileups_gb['reads'].sum()})

pileups_marginal['mean_alignment'] /= pileups_marginal['reads'] # Normalize mean

base_scores = pd.DataFrame(
    {name: scores.groupby(level=name).sum() for name in align_pairs}
)
normed_scores = base_scores.div(max_scores, axis=1)


normed_scores['pileups'] = np.histogram(
    pileups_marginal['mean_alignment']
    weights=pileups_marginal['reads'],
    bins=len(normed_scores),
    range=(0,1)
    )[0]

print(normed_scores.describe())
base_scores.columns.name = "alignment"
longform = base_scores.stack()
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

def length(bc):
    return pd.Series(bc).str.len().std()

from scipy.stats import entropy
p_stats = pileups.groupby(level='score').agg(dict(barcode=length_stdev, reads=entropy)
p_stats['Fwd-Rev Alignment'] = scores.groupby('Core-Ref')['Fwd-Rev'].mean() / max_scores['Fwd-Rev']
p_stats.columns.name = 'stat'
p_stats /= p.stats.max()/100

line_data = p_stats.stack()
line_data.name = '% of Max'



sns.lineplot(
    x='score',
    y=line_data.name,
    hue='stat'


[ pileups_gb.groups.keys()]

# What we want:
#   - self-similarity declines with alignment quality
#   - Lengths diversify
#   - barcodes reduce in diversity
#
#

print('f {:}'

hue_order = 
sns.barplot(
    x='hample',
    y'=reads',
    hue='total'




