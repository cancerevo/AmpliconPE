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


base_scores = pd.DataFrame(
    {name: scores.groupby(level=name).sum() for name in align_pairs}
)
normed_scores = base_scores.div(max_scores, axis=1)
print(normed_scores.describe())
base_scores.columns.name = "alignment"
longform = base_scores.stack()
longform.name = "reads"

print(base_scores)


sns.ecdfplot(
    data=longform.reset_index(),
    x="score",
    weights="reads",
    hue="alignment",
    stat="count",
    hue_order=align_pairs,
    ax=plt.gca(),
)
plt.savefig("alignments.pdf")


ref_scores = scores.groupby(level=["Fwd-Ref", "Rev-Ref"]).sum()
composite_scores = ref_scores.groupby(lambda ix: ix[0] + ix[1]).sum()
