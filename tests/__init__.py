from AmpliconPE import MasterRead, BarcodeSet, pairedFASTQiter, reverse_compliment

from os import Path

base_dir = Path.getcwd()


def test_pairedFASTQiter():
    examples = base_dir / 'examples'
    fastq = '77.fastq.gz'
    Iter = pairedFASTQiter(examples / 'forward_reads' / fastq, examples / 'reverse_reads' / fastq)
    next(Iter)
    next(Iter)
    fwd_read, rev_read = next(Iter)
    assert fwd_read == 'AAA'
    assert rev_read == 'AAA'


def test_MasterRead():
    read = 2*'ACGT' + 2*"NNNNAA" + 2*"ACGT"
    rc = reverse_compliment(read)
    master_read = MasterRead(read)
    assert master_read.max_score == 90
    fwd_read = (read[:4] + read[5:]).replace("N", "T")
    rev_read = (read[:-3] + 'A' + read[-2:]).replace('N', 'T')

    Alignment = master_read.align(bytes(fwd_read), bytes(rev_read))
    assert Alignment.score == 100
    assert Alignment.extract_barcode() == 'TTTTAATTTT'

def test_BarcodeSet():
    barcodes = {10*'A':'As', 10*'C':'Cs', 10*'G':'Gs', 10*'T':'Ts'}
    barcode_map = BarcodeSet(barcodes, n_mismatches=3)
    assert barcode_map[5*'A'+3*'G'+2*'A'] == 'As'


