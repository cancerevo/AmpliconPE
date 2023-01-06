from AmpliconPE import MasterRead, BarcodeSet, pairedFASTQiter, reverse_compliment

from pathlib import Path

base_dir = Path.cwd()


def test_pairedFASTQiter():
    examples = base_dir / "examples"
    fastq = "77.fastq.gz"
    Iter = pairedFASTQiter(
        examples / "forward_reads" / fastq, examples / "reverse_reads" / fastq
    )
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
