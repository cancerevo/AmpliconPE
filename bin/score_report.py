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
Histograms  -> 1,
                2
Ratios      -> 3, 
%           -> 4, 5
Totals      -> 6
sgID % (2-ways)
"""

def histogram(data, ax):

def plot_

