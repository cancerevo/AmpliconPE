"""

Questions-Evidence (E) pairs:
------------------

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

N_plots = 5


base_scores = pd.DataFrame({
    'Fwd-Ref':scores.sum(axis=2).sum(axis=1),
    'Rev-Ref':scores.sum(axis=2).sum(axis=0),
    'Fwd-Rev':scores.sum(axis=1).sum(axis=0)
}, index=pd.Index(np.linspace(0, 1, len(scores), name='score'))
)
base_scores.columns.name = 'alignment'

composite_scores = np.zeros(2*len(scores) - 1, dtype=np.int64)
for fwd_score, row in enumerate(scores.sum(axis=2)):
    for rev_score, reads in enumerate(row)]):
        composite_score[fwd_score+rev_score] += reads

longform = pd.concat(dict(
    base=base_scores.stack(),
    composite=pd.Series(composite_scores, index=pd.Index(np.linspace(0, 1, len(composite_scores), name='score')))),
    names=['Class']
)
sns.histplot(data=base_scores.stack().reset_index(), x='score', y='reads', discrete=True, hue='alignment', stat='percent', multiple='dodge', 
    hue_order=['Fwd-Ref', 'Rev-Ref', 'Fwd-Rev'])

