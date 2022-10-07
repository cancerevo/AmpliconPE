# AmpliconPE
Extracts DNA barcodes from Paired-End Amplicon Sequencing using Smith Waterman alignments

## Why use AmpliconPE for your barcoding project?

* Paried-end (redundant) seqencing of amplicons is cheaper and more accurate
* Homology-based alignment to a reference read is the most robust way to correct all read anomalies 
* AmpliconPE extracts both known barcodes (barcode sets) and random barcodes intelligently
* Its fast & simple

## Usage

AmpliconPE is a highly-extensible python library packaged with convient command-line scripts to 
process common amplicon constructs (e.g. TuBa-seq, CLONtracer, Brunello/Brie CRISPRko libraries). 

Usage is built around a `MasterRead` class, which aligns Paired-End (PE) FASTQ reads to a reference
sequence

```python


from AmpliconPE import MasterRead
tuba_seq_reference_seq = 'GACCCGGA'            +    # 5' flanking sequence of double-barcode (8 nts is good)
                         '********'            +    # Known barcode - specified by '*'
                         'AA'                  +    # spacer
                         'NNNNNTTNNNNNAANNNNN' +    # Random Barcode w/ spacers - specified by 'N'
                         'ATGCCCAA'
                       
master_read = MasterRead(tuba_seq_reference_seq)
print(vars(master_read))
# {'sequence': '', 
# 'max_score': XX, 
# 'rc_sequence': 
# 'SW_kwargs': {}, 
# 'known_barcode':True, 
# 'random_barcode':True, 
# 'known_slice':slice(8, 16, None),
# 'random_slice':slice(18,37))}
```
Known barcodes are assigned labels using the `BarcodeSet` --- a sub-class of `dict` that tolerates 
mismatches when assigning labels. Use [BARCOSEL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2262-7)
to generate barcodes that can tolerate sequencing errors. 

```python
from AmpliconPE import BarcodeSet
import pandas as pd

sgIDs = pd.read_csv('sgID_info.csv').set_index("target")['ID']
known_barcodes = BarcodeSet(sgIDs, n_mismatches=1, indel=1)
```
Reads are then processed using the `IterPairedFASTQ` iterator:

```python
from AmpliconPE import IterPairedFASTQ

Iter = IterPairedFASTQ('forward_file.fastq.gz', 'reverse_file.fastq.gz', check_indecies=True)

FWD_read, REV_read = next(Iter)
```
Index-Hopping dramatically undermines the power of barcode sequencing. Be sure to (1) sequence only primer-free libraries, 
(2) use Dual-Unique Indecies, and (3) filter Forward/Reverse index pairs that do not match. 

Barcodes are then extracted using an internal-loop that generally looks something like this: 

```python

  score = master_read.align(FWD_read, REV_read)
  if score < 0.8 * master_read.max_score:
    continue  # poor alignment

  sgID = known_barcodes.get(master_read.extract_known_barcode(), 'Unknown sgID')
  random_barcode = master_read.extract_random_barcode()
  if 'N' in random_barcode or random_barcode == 'Length Mismatch':
    continue
```

In short, `MasterRead` internally-stores the reference alignment to both the forward and reverse read. Its methods process 
both reads in tandem. Reads are kept if they align to the reference read well and if the barcodes between the forward and revese reads match. 
Single nucleotide differences between reads are replaced with an 'N', when an InDel difference generally return a 'Length Mismatch' string. 

Barcode pileups are then tallied and processed using downstream software. You may want to use a barcode clusterer, e.g. Shepard, to de-noise
random barcodes; however, PE sequencing generally resolves most reccurrent-read errors. The function `AmpliconPE.identify_neighbors` cleans-up
simple recurrent read errors. 
