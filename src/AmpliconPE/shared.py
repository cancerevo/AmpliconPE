import numpy as np
import pandas as pd

DTYPES = dict(
    stat="O",
    sample="O",
    value=np.int64,
    barcode="O",
    reads=np.int64,
    score=np.int64,
)

nucleotides = "ACGTN"


def barcode_content(barcodes):
    return pd.concat(
        {nuc: barcodes.str.count(nuc) for nuc in nucleotides}, names=["content"]
    )
