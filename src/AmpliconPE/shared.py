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


def isTransition(pair_set):
    return pair_set == {"A", "G"} or pair_set == {"C", "T"}


mutation_categories = {
    pair: "transition" if isTransition(set(pair)) else "transversion"
    for pair in combinations(nucleotides[:-1], 2)
}

_extra_categories = {"N": "read mismatch", "Ins": "Insertion", "Del": "Deletion"}

mutation_categories.update(
    {
        (nucleotide, annotation): mut_category
        for annotation, mut_category in _extra_categories.items()
        for nucleotide in nucleotides
    }
)

mutation_categories.update(
    {(pair[1], pair[0]): category for pair, category in mutation_categories.items()}  #
)  # Make Symmetric
mutation_categories = pd.Series(mutation_categories, name="category")
mutation_categories.index.names = ["From", "To"]

mutation_categories = mutation_categories.loc[
    list(nucleotides) + ["Ins"], list(nucleotides) + ["Del"]
]


def barcode_content(barcodes):
    return pd.concat(
        {nuc: barcodes.str.count(nuc) for nuc in nucleotides}, names=["content"], axis=1
    )


def delta_content_constructor(row):
    change = {nuc: 0 for nuc in nucleotides}
    if row["From"] != "Ins":
        change[row["From"]] = -1
    if row["To"] != "Del":
        change[row["To"]] = +1
    return pd.Series(change)


delta_content_map = mutation_categories.reset_index().apply(
    delta_content_constructor, axis=1
)
delta_content_map.index = mutation_categories.index


def annotate_mutations(barcodes):
    from_barcodes = barcodes["From"]
    to_barcodes = barcodes["To"]
    d_content = barcode_content(to_barcodes) - barcode_content(from_barcodes)
    annotations = pd.DataFrame(
        dict(From=d_content.idxmin(axis=1), To=d_content.idxmax(axis=1))
    )
    d_nucleotides = d_content.sum(axis=1)
    insertions = d_nucleotides == +1
    deletions = d_nucleotides == -1
    annotations.loc[insertions, "From"] = "Ins"
    annotations.loc[deletions, "To"] = "Del"
    return annotations
