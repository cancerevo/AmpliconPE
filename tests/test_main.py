from AmpliconPE import (
    MasterRead,
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


def test_MasterRead():
    read = 2 * "ACGT" + 2 * "NNNNAA" + 2 * "ACGT"
    master_read = MasterRead(read)
    assert master_read.max_score == 120
    rc = reverse_compliment(read.encode("ascii"))

    perfect_alignment = master_read.align(
        read.replace("N", "T").encode("ascii"), rc.replace(b"N", b"A")
    )
    assert perfect_alignment.score == 120
    assert perfect_alignment.extract_barcode() == "TTTTAATTTT"

    fwd_read = (read[:4] + read[5:]).replace("N", "T")
    rev_read = (rc[:-3] + b"A" + rc[-2:]).replace(b"N", b"T")
    imperfect_alignment = master_read.align(fwd_read.encode("ascii"), rev_read)
    assert imperfect_alignment.score == 107
    assert imperfect_alignment.extract_barcode() == "NNNAANNNNA"


def test_BarcodeSet():
    barcodes = {10 * "A": "As", 10 * "C": "Cs", 10 * "G": "Gs", 10 * "T": "Ts"}
    barcode_map = BarcodeSet(barcodes, n_mismatches=3)
    assert barcode_map[5 * "A" + 3 * "G" + 2 * "A"] == "As"


def test_TuBa_pileups():
    new = pd.read_hdf(data_dir / "consolidated_outputs.h5", "pileups")
    reference = pd.read_csv(
        data_dir / "reference_pileups.csv.gz",
        index_col=["sample", "target", "barcode"],
        dtype=dict(sample=object),
    )
    joined = reference.join(new, on=["sample", "target", "barcode"], rsuffix="new")
    R = joined.corr().iloc[0, 1]
    assert 1 - R < 1e-6, "New Pileups are appreciably different than old pileups"
