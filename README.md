# AmpliconPE
Extracts DNA barcodes from Paired-End (PE) amplicon sequencing using Smith Waterman alignments

## Why use AmpliconPE for your barcoding project?

From [Johnson et al 2022][0]:
> Extracting barcodes from the sequencing reads may appear as a trivial problem at first glance, given that the structure of the read is known by design. 
> However, the challenge is that not all reads may have identical structure...

* PE (redundant) seqencing is cheaper (!) than Single-End sequencing and eliminates 99% of errors
* Homology-based alignment to a reference read is the most robust way to correct read anomalies 
* Its _SIMD_-fast, simple & _ultra_-extensible library.

## Why _not_ use a RegEx expression?

Fuzzy RegEx libraries use sequence aligners, which _could_ mimic our functionality; however,
this package is built around years of expreience addressing biases and contaminations in barcoding project. 
Also, please consider merging your forward and reverse reads using [PEAR][1] before applying a RegEx expression. 

## Usage

We provide command-line scripts to process common amplicon constructs (e.g. TuBa-seq, CLONtracer, Brunello/Brie CRISPRko libraries). 

Usage is built around a `MasterRead` class, which aligns PE reads to a reference sequence:

```python


from AmpliconPE import MasterRead
tuba_seq_reference_seq = 'GACCCGGA'            +    # 5' flanking sequence of double-barcode (8 nts is good)
                         '********'            +    # Known barcode - specified by '*'
                         'AA'                  +    # spacer
                         'NNNNNTTNNNNNAANNNNN' +    # Random Barcode w/ spacers - specified by 'N'
                         'ATGCCCAA'                 # 3' flank
                       
master_read = MasterRead(tuba_seq_reference_seq)
```
MasterRead automatically extracts barcode start and end locations...
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

taget_genes = pd.read_csv('sgRNA_info.csv').set_index("Targeted Gene")['barcode']
known_barcodes = BarcodeSet(target_genes, n_mismatches=1, indel=1)
```
Reads are then processed using the `IterPairedFASTQ` iterator:

```python
from AmpliconPE import IterPairedFASTQ

Iter = IterPairedFASTQ('forward_file.fastq.gz', 'reverse_file.fastq.gz', check_indecies=True)

FWD_read, REV_read = next(Iter)
```
Index-Hopping dramatically undermines the power of barcode sequencing. For this reason, we recommend 
(1) sequence only primer-free libraries, (2) use Dual-Unique Indecies, and (3) filter non-matching 
Forward/Reverse index pairs using the _check_indecies_ keyword argument. 

Barcodes are then extracted using an internal-loop that generally looks something like this: 

```python

  score = master_read.align(FWD_read, REV_read)
  if score < 0.8 * master_read.max_score:
    continue  # poor alignment

  target = known_barcodes.get(master_read.extract_known_barcode(), 'Unknown Target')
  random_barcode = master_read.extract_random_barcode()
  if 'N' in random_barcode or random_barcode == 'Length Mismatch':
    continue
```

In short, `MasterRead` internally-stores the reference alignment to both the forward and reverse read when you call the `align` method--- alignment 
is the performance-limiting step and implemented using [SSW Library][3]. Single nucleotide differences between the forward and reverse barcodes 
are replaced with `'N'`, and an InDel difference between barcodes returns `'Length Mismatch'`. 

Generally, reads are kept if they match the reference read well and if forward and reverse barcodes match. Barcode pileups are then tallied and 
processed using downstream software. You may want to use a barcode clusterer, e.g. Shepard, to de-noise random barcodes; however, PE sequencing 
generally resolves most reccurrent-read errors. The function `AmpliconPE.identify_neighbors` cleans-up simple recurrent read errors. 

[0]: https://ecoevorxiv.org/t58xw/ "Johnson, M. S., Venkataram, S., & Kryazhimskiy, S. (2022, September 28). Best practices in designing, sequencing and identifying random DNA barcodes."

[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3933873/ "Zhang J, Kobert K, Flouri T, Stamatakis A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics. 2014 Mar 1;30(5):614-20. doi: 10.1093/bioinformatics/btt593. Epub 2013 Oct 18. PMID: 24142950; PMCID: PMC3933873."

[2]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2262-7 "Somervuo, P., Koskinen, P., Mei, P. et al. BARCOSEL: a tool for selecting an optimal barcode set for high-throughput sequencing. BMC Bioinformatics 19, 257 (2018)."

[3]: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138 "Zhao, Mengyao, et al. SSW library: an SIMD Smith-Waterman C/C++ library for use in genomic applications. PloS one 8.12 (2013): e82138."
