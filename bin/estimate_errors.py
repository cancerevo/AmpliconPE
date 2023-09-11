#!/bin/env python3
import pandas as pd
from AmpliconPE.shared import nucleotides, barcode_content
from AmpliconPE import BarcodeSet
from pathlib import Path
import numpy as np
from jsonargparse import CLI


def estimate_errors(aligned_filename: Path, model_filename: Path = "error_model.csv"):

    model = pd.read_csv(model_filename)

    # Converting barcode pairs into a From->To mutation pair is performance limiting.
    # Hashing into an ndarray helps...

    def delta_content(row):
        return pd.Series(
            [
                -1 if row.From == nucleotide else (+1 if row.To == nucleotide else 0)
                for nucleotide in nucleotides
            ],
            index=list(nucleotides),
        )

    d_content_model = model.apply(delta_content, axis=1)

    HASHING_VECTOR = 3 ** np.arange(len(nucleotides))
    model_hashes = (d_content_model.values).dot(HASHING_VECTOR)
    hashed_model = np.zeros(2 * model_hashes.max() + 1)
    hashed_model[model_hashes] = model.Probability.values
    assert hashed_model.sum() == model.Probability.sum()

    aligned = pd.read_csv(aligned_filename)
    pileups = aligned.set_index("barcode")["reads"]

    mismatcher = BarcodeSet(InDels=True, fixed_length=False).mismatcher
    possible_neighbors = frozenset(aligned.barcode)

    content = barcode_content(aligned.barcode)
    content_map = content.set_index(aligned.barcode)

    def sequencing_errors(barcode):
        neighbors = pileups.loc[
            (
                mismatch
                for mismatch in mismatcher(barcode)
                if mismatch != barcode and mismatch in possible_neighbors
            )
        ]
        d_content = content_map.loc[barcode] - content_map.loc[neighbors.index]
        return hashed_model[d_content.values.dot(HASHING_VECTOR)].dot(neighbors)

    aligned["sequencing errors"] = aligned["barcode"].apply(sequencing_errors)
    error_loss = model.groupby("From")["Probability"].sum()
    Ins = error_loss.pop("Ins")
    error_loss += Ins

    aligned["lost reads"] = (
        content.mul(error_loss, axis=1).sum(axis=1) * aligned["reads"]
    )

    from scipy.stats.distributions import poisson

    aligned["P_sequencing_error"] = 1 - poisson.cdf(
        aligned.reads, aligned["sequencing errors"]
    )
    aligned.to_csv(aligned_filename, index=False)


if __name__ == "__main__":
    CLI()
