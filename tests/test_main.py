from AmpliconPE import (
    MasterRead,
    SimplexMasterRead,
    BarcodeSet,
    pairedFASTQiter,
    reverse_compliment,
    get_PE_FASTQs,
)
import pandas as pd
from pathlib import Path

base_dir = Path.cwd()
data_dir = base_dir / "data"


def test_pairedFASTQiter():
    directory = data_dir / "TuBa-seq_FASTQs" / "77"

    fwd_file, rev_file = get_PE_FASTQs(directory)
    Iter = pairedFASTQiter(fwd_file, rev_file)
    next(Iter)
    next(Iter)
    real_fwd_read = b"ACTAGGTGCGCACGTCTGCCGCGCTGTTCTCCTCTTCCTCATCTCCGGGACCCGGGATGCCCAAGAAGAAGAGGAAGGTGTCCAATTTACTGACCGTACACCAAAATTTGCCTGCATTACCGGTCGATGCAACGAGTGATGAGGTTCGCAA"
    real_rev_read = b"CTGGTCNTCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNGGCATCCCGGGTCCCGGAGATGAGGAAGAGGAGAACAGCGCGGCAG"
    fwd_read, rev_read = next(Iter)
    assert fwd_read == real_fwd_read
    assert rev_read == real_rev_read


read = 2 * "ACGT" + 2 * "NNNNAA" + "AACCGGTT"
alignment_params = dict(match=6, mismatch=-2, gap_open=6, gap_extend=1)

fwd_read = read.encode("ascii")
rev_read = reverse_compliment(fwd_read)

perfect_seq_pair = fwd_read.replace(b"N", b"T"), rev_read.replace(b"N", b"A")
perfect_barcode = "TTTTAATTTT"


imperfect_seq_pair = (fwd_read[:4] + fwd_read[5:]).replace(b"N", b"T"), (
    rev_read[:-3] + b"A" + rev_read[-2:]
).replace(b"N", b"T")

imperfect_barcode_error = "NNNAANNNNA"
perfect_imperfect_barcode = "NNNNAANNNN"


def test_MasterRead():
    perfect_score = (len(read) - read.count("N")) * alignment_params["match"] * 2
    imperfect_score = (
        perfect_score
        - alignment_params["match"]
        + alignment_params["mismatch"]
        - 2 * alignment_params["gap_open"]
    )

    master_read = MasterRead(read, **alignment_params)
    assert master_read.max_score == perfect_score
    perfect_alignment = master_read.align(*perfect_seq_pair)
    assert perfect_alignment.final_score == master_read.max_score
    assert perfect_alignment.extract_barcode() == perfect_barcode

    imperfect_alignment = master_read.align(*imperfect_seq_pair)
    assert imperfect_alignment.final_score == imperfect_score
    assert imperfect_alignment.extract_barcode() == imperfect_barcode_error


def test_simplexMasterRead():
    perfect_scores = (120, 120, 168, 42)
    imperfect_score = 42
    master_read = SimplexMasterRead(read, **alignment_params)
    perfect_alignment = master_read.align(*perfect_seq_pair)

    assert perfect_alignment.fwd_read == perfect_seq_pair[0]
    assert perfect_alignment.rev_read == perfect_seq_pair[1]
    assert perfect_alignment.core_consensus == perfect_seq_pair[0]

    assert perfect_alignment.get_scores() == perfect_scores
    assert perfect_alignment.extract_barcode() == perfect_barcode

    imperfect_alignment = master_read.align(*imperfect_seq_pair)
    assert imperfect_alignment.final_score == imperfect_score
    assert imperfect_alignment.extract_barcode() == perfect_imperfect_barcode


def test_BarcodeSet():
    barcodes = {10 * "A": "As", 10 * "C": "Cs", 10 * "G": "Gs", 10 * "T": "Ts"}
    barcode_map = BarcodeSet(barcodes, n_mismatches=3)
    assert barcode_map[5 * "A" + 3 * "G" + 2 * "A"] == "As"


"""
def test_TuBa_pileups():
    pileups = (data_dir / "TuBa-seq_FASTQs").glob("**/pileups.csv")
    new = pd.concat(
        {
            filename.parts[-2]: pd.read_csv(filename, index_col=["target", "barcode"])
            for filename in pileups
        },
        names=["Sample"],
    )["reads"]
    reference = pd.read_csv(
        data_dir / "reference_pileups.csv.gz",
        index_col=["sample", "target", "barcode"],
        dtype=dict(sample=object),
    )
    joined = reference.join(new, on=["sample", "target", "barcode"], rsuffix="new")
    R = joined.corr().iloc[0, 1]
    assert 1 - R < 1e-6, "New Pileups are appreciably different than old pileups"

"""
