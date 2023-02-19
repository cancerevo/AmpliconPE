#!/usr/bin/env python3
"""
All functions in filters return a boolean array of barcode pileups to ~keep~
"""

import pandas as pd
import numpy as np

filters = [
    "Insufficient_Reads",
    "Singletons",
    "Unexpected_sgID",
    "Aberrant_Length",
    "Unknown_sgID",
    "Recurrent_Sequencing_Error",
    "Neighbor_Contamination",
    "iHopping",
]


def Singletons(S):
    return S == 1


def Insufficient_Reads(S):
    gb = S.groupby(level="Sample")
    depth = gb.size()
    bad_samples = depth.loc[depth < args.min_read_depth].index
    if len(bad_samples) > 0:
        info = stats.loc[bad_samples]
        Log("The following Samples had insufficient read depth:", True)
        Log(
            info.div(info.sum(axis=1))
            .assign(**{"Passed (#)": info["Passed"].astype(int)})
            .to_string(float_format="{:.0%}".format),
            True,
        )
    return gb.transform(lambda s: s.name in bad_samples)


unexpected_sgIDs = list(set(other_sgIDs.values()) - set(sgID_map.values()))


def Unexpected_sgID(S):
    return S.index.get_level_values("target").isin(unexpected_sgIDs)


def Aberrant_Length(S):
    return (
        np.abs(
            S.index.get_level_values("barcode").str.len()
            - (barcode_length - sgID_length)
        )
        > args.indel_tolerated
    )


def Unknown_sgID(S):
    return ~S.index.get_level_values("target").isin(set(all_sgIDs.values()))


def Recurrent_Sequencing_Error(S):
    return groupby_apply(
        S.groupby(level=["Sample", "target"], group_keys=False), identify_neighbors
    )


def Neighbor_Contamination(S):
    M = S.unstack(level="Sample")
    order = list(M.columns)
    N_mice = len(M.columns)
    contaminants, _report_str = contamination(
        M, alpha=0.05 * (N_mice * 100 - 1) / 2, map=map, graph_contaminations=False
    )
    Ixs = S.isnull()
    for ix, row in contaminants.iterrows():
        contaminant, contaminator = ix
        if abs(order.index(contaminant) - order.index(contaminator)) == 1:
            print("Spillover between:", *ix)
            toDrop = (
                M[contaminant]
                .loc[lambda S: S < M[contaminator] * 0.75]
                .nsmallest(
                    int(row["Overlapping Barcodes"] - row["Expected Barcode Overlap"])
                )
                .sort_index()
            )
            Ixs.loc[contaminant].loc[toDrop.index] = True
    return Ixs


def iHopping(S):
    from scipy.stats.distributions import poisson
    from .estimates import iHop_rate

    barcode_matrix = S.unstack("Sample")
    if len(barcode_matrix.columns) < 2:
        Log("Cannot investigate i-Hopping with only 1 sample.")
        return S < -1
    error_rate = iHop_rate(barcode_matrix)
    combined_errors["Index-Hopping"] = error_rate
    if type(error_rate) is str:
        return S > -1
    largest = barcode_matrix.max(axis=1)
    min_tally = pd.Series(
        poisson.ppf(1 - args.pHop, largest * error_rate), index=largest.index
    ).fillna(0)
    if "Spike" in min_tally.index.get_level_values("target").unique():
        min_tally.loc["Spike"] = 0
    return (
        barcode_matrix.lt(min_tally, axis=0)
        .stack()
        .reorder_levels(["Sample", "target", "barcode"])
    )
