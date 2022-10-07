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
                         '^^^^^^^^'            +    # Known barcode - specified by '^'
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
