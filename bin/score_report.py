#!/usr/bin/env python3
"""

Questions-Evidence (E) pairs:
------------------

1. Histograms of 3-ref pairs

2. Fwd-Rev versus Combined; 

Combined vs. Singleton 
Combined vs. Fwd-Rev

1. Is [forward/reverse] sequencing quality good?
    E: Forward reads match reverse?
       Forward reads match reference?
       Reverse reads match reference?


2. Is the library Amplicon?
    E: High unaligned rate -> contamination
    E: Flat Cap -> wrong reference

3. Did any samples have poor PCR?
    E: Fwd-Rev match, but not ref (for aligned)
        E[N] / E[Max_Score - Score] : aligned

4. Did any samples have poor index matching?
    E: Index mismatch %
        
5. Did any samples have cloning problems?
    E: Length disagreement %
    
6. Did any samples have poor coverage?
    E: Read coverage

7. What is the sgID representation?
    E: sgID coverage (pileups & counts)

LAYOUT:
Histograms  -> 1, 2
Ratios      -> 3, 
%           -> 4, 5
Totals      -> 6
sgID % (2-ways)
"""
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

N_plots = 5

score_pairs = ["Fwd-Ref", "Rev-Ref", "Fwd-Rev"]
scores = pd.read_csv(
    "consolidated_scores.csv", index_col=["sample"] + score_pairs, dtype=np.uint64
)["reads"]

# print(scores, scores.describe(), len(scores)/len(scores.dropna()))
base_scores = pd.DataFrame(
    {name: scores.groupby(level=name).sum() for name in score_pairs}
)
print(base_scores.describe(), len(base_scores) / len(base_scores.dropna()))
# .astype(int)
base_scores.index = pd.Index(np.linspace(0, 1, len(base_scores)), name="score")
base_scores.columns.name = "alignment"
longform = base_scores.stack()
longform.name = "reads"

print(base_scores)

sns.histplot(
    data=longform.reset_index(),
    x="score",
    y="reads",
    discrete=True,
    hue="alignment",
    stat="percent",
    # multiple="dodge",
    hue_order=score_pairs,
    ax=plt.gca(),
)
plt.savefig("Alignments.pdf")


ref_scores = scores.groupby(level=["Fwd-Ref", "Rev-Ref"]).sum()
composite_scores = ref_scores.groupby(lambda ix: ix[0] + ix[1]).sum()
