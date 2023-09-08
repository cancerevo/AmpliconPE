import pandas as pd
import numpy as np

alphabet = "ACGTN"


def identify_neighbors(
    Start,
    alpha=1e-10,
    error_rate=1e-3,
    indels=True,
    max_del_length=2,
    add_doubles=False,
):
    """
    For a typical run, the breakdowns of InDels are:
      true   87.83%
        -1    5.94%
        +1    2.30%
        -2    1.81%
        -3    1.18%
        -4    0.73%
        +2    0.20%

    i.e. Deletions are more common (and a much smaller space of barcodes) than insertions. If fact, double deletions are about as likely as single insertions.
    """
    from scipy.stats.distributions import poisson

    S = Start.loc[Start.index.values[0][:2]].sort_values(ascending=False)
    barcodes = S.index
    drop = np.zeros(len(S), dtype=np.bool)
    spawn_ix = 0
    for barcode, abundance in S.iteritems():
        max_tally = int(poisson.ppf(1 - alpha, abundance * error_rate))
        while spawn_ix < len(S) and S.iloc[spawn_ix] > max_tally:
            spawn_ix += 1
        SNVs = {
            barcode[:i] + nuc + barcode[i + 1 :]
            for nuc in alphabet
            for i, bc_nuc in enumerate(barcode)
            if nuc != bc_nuc
        }
        candidate_spawn = SNVs
        if indels:
            deletions = {
                barcode[:i] + barcode[i + del_length :]
                for del_length in range(1, max_del_length + 1)
                for i in range(len(barcode) - del_length + 1)
            }
            candidate_spawn |= deletions | {
                barcode[:i] + nuc + barcode[i:]
                for i in range(len(barcode) + 1)
                for nuc in alphabet
            }
            if add_doubles:
                candidate_spawn |= (
                    # double point mutations
                    {
                        SNV[:i] + nuc + SNV[i + 1 :]
                        for nuc in alphabet
                        for SNV in SNVs
                        for i, snv_nuc in enumerate(SNV)
                        if nuc != snv_nuc
                    }
                    |
                    # Deletions plus point mutations
                    {
                        deletion[:i] + nuc + deletion[i + 1 :]
                        for nuc in alphabet
                        for deletion in deletions
                        for i, del_nuc in enumerate(deletion)
                        if nuc != del_nuc
                    }
                )
        drop[spawn_ix:] |= barcodes[spawn_ix:].isin(candidate_spawn)
    return pd.Series(drop, index=Start.index)
