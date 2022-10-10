# AmpliconPE
Extracts DNA barcodes from Paired-End (PE) amplicon sequencing using Smith Waterman alignments

## Why use AmpliconPE for your barcoding project?

* PE (redundant) seqencing of amplicons is cheaper and more accurate
* Homology-based alignment to a reference read is the most robust way to correct all read anomalies 
* AmpliconPE extracts both known barcodes (barcode sets) and random barcodes intelligently
* Its fast & simple

## Why _not_ use a RegEx expression?

Fuzzy RegEx libraries use sequence aligners, so you _could_ mimic our functionality; however,
this package is built around years of expreience addressing biases and contaminations in barcoding project. 
Also, please consider merging your forward and reverse reads using [PEAR][1] before applying a RegEx expression. 

## Usage

AmpliconPE is an extensible python library. We provide command-line scripts to process common 
amplicon constructs (e.g. TuBa-seq, CLONtracer, Brunello/Brie CRISPRko libraries). 

Usage is built around a `MasterRead` class, which aligns PE reads to a reference sequence:

```python


from AmpliconPE import MasterRead
tuba_seq_reference_seq = 'GACCCGGA'            +    # 5' flanking sequence of double-barcode (8 nts is good)
                         '********'            +    # Known barcode - specified by '*'
                         'AA'                  +    # spacer
                         'NNNNNTTNNNNNAANNNNN' +    # Random Barcode w/ spacers - specified by 'N'
                         'ATGCCCAA'
                       
master_read = MasterRead(tuba_seq_reference_seq)
```
MasterRead automatically extracts barcode locations...
```python
print(master_read.known_barcode)
> slice(8, 16, None)
print(master_read.random_barcode)
> slice(18, 37 None)
```
You can map known barcodes to labels using `BarcodeSet`, a sub-class of `dict` that tolerates 
mismatches when assigning labels. Use [BARCOSEL][2] to generate barcode sets that are robust to sequencing errors. 

```python
from AmpliconPE import BarcodeSet
import pandas as pd

sgIDs = pd.read_csv('sgID_info.csv').set_index("Targeted Gene")['barcode']
known_barcodes = BarcodeSet(sgIDs, n_mismatches=1, indel=1)
```
Reads are then processed using the `IterPairedFASTQ` iterator:

```python
from AmpliconPE import IterPairedFASTQ

Iter = IterPairedFASTQ('forward_file.fastq.gz', 'reverse_file.fastq.gz', check_indecies=True)

FWD_read, REV_read = next(Iter)
```
Index-Hopping dramatically undermines the power of barcode sequencing. Be sure to (1) sequence only primer-free libraries, 
(2) use Dual-Unique Indecies, and (3) filter Forward/Reverse index pairs that do not match using the _check_indecies_ keyword argument. 

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

In short, `MasterRead` internally-stores the reference alignment to both the forward and reverse read, as sequence alignment is the performance-
limiting step (implemented as embeded C code). . Its methods process 
both reads in tandem. Reads are kept if they align to the reference read well and if the barcodes between the forward and revese reads match. 
Single nucleotide differences between reads are replaced with an 'N', when an InDel difference generally return a 'Length Mismatch' string. 

Barcode pileups are then tallied and processed using downstream software. You may want to use a barcode clusterer, e.g. Shepard, to de-noise
random barcodes; however, PE sequencing generally resolves most reccurrent-read errors. The function `AmpliconPE.identify_neighbors` cleans-up
simple recurrent read errors. 

# References

[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933873/ "Zhang J, Kobert K, Flouri T, Stamatakis A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics. 2014 Mar 1;30(5):614-20. doi: 10.1093/bioinformatics/btt593. Epub 2013 Oct 18. PMID: 24142950; PMCID: PMC3933873."

[2]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2262-7 "Somervuo, P., Koskinen, P., Mei, P. et al. BARCOSEL: a tool for selecting an optimal barcode set for high-throughput sequencing. BMC Bioinformatics 19, 257 (2018). https://doi.org/10.1186/s12859-018-2262-7"
