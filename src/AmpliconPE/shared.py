import numpy as np
import pandas as pd
from itertools import combinations

DTYPES = dict(
    stat="O",
    sample="O",
    value=np.int64,
    barcode="O",
    reads=np.int64,
    score=np.int64,
)

nucleotides = "ACGTN"

mut_categories = {
    pair: "transition"
    if set(pair) == {"A", "G"} or set(pair) == {"C", "T"}
    else "transversion"
    for pair in combinations(nucleotides[:-1], 2)
}
mut_categories.update({(nuc, "N"): "read mismatch" for nuc in nucleotides[:-1]})
mut_categories.update(
    {(pair[1], pair[0]): category for pair, category in mut_categories.items()}
)  # Make Symmetric
mut_categories = pd.Series(mut_categories, name="category")
mut_categories.index.names = ["From", "To"]


def barcode_content(barcodes):
    return pd.concat(
        {nuc: barcodes.str.count(nuc) for nuc in nucleotides}, names=["content"], axis=1
    )


def annotate_mutations(from_barcodes, to_barcodes):
    annotations = (
        barcode_content(to_barcodes) - barcode_content(from_barcodes)
    ).assign(
        From=lambda df: df.idxmin(axis=1),
        To=lambda df: df.idxmax(axis=1),
        Type=lambda df: df.sum(axis=1).map(
            {-1: "Deletion", +1: "Insertion", 0: "Substitution"}
        ),
    )
    substitutions = annotations["Type"] == "Substitution"
    annotations.loc[substitutions] = annotations.loc[substitutions, ["From", "To"]].map(
        mut_categories
    )
    assert not annotations["Type"].isnull().any()
    return annotations
