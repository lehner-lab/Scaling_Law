# Pipeline

This document describes the pipeline used to process the raw sequencing files and count the number of times each genotyope shows up in our HEK293 phylogeny library. To run this pipeline, you'll need the following software installed:

* FastQC ([download](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
* Sabre ([download](https://github.com/najoshi/sabre))
* PEAR ([download](https://sco.h-its.org/exelixis/web/software/pear/))
* Fastq-Tools ([download](https://github.com/dcjones/fastq-tools))
* FASTX-Toolkit ([download](http://hannonlab.cshl.edu/fastx_toolkit/))
* Seqtk ([download](https://github.com/lh3/seqtk))

All the code in this document is written in bash unless stated otherwise.

## 1. Quality control

We first checked that the file looked OK using FastQC. You can run the following code for this:

```bash
bash 001_quality_control.sh
```

## 2. Demultiplex

The nine biological replicates ("output" from Figure 1C) and the three technical replicates ("input" from Figure 1C) were all sequenced in the same lane. To distinguish all these replicates from each other, [barcodes](./002_barcodes.txt) were added to the sequenced reads. We used Sabre to demultiplex the replicates:

```bash
bash 002_demultiplex_replicates.sh
```

The three technical replicates and the first three biological replicates have the same barcode and so will be demultiplexed together. In the merging step (below) we will separate them.

**Note that you can skip this step if you downloaded the already demultiplexed files from GEO.**


## 3. Merge forward & reverse reads

Up to now we have been dealing with paired-end reads. In this step we use PEAR to combine forward and reverse reads:

```bash
bash 003_merge_forward_and_reverse_seqs
```

**Note that you can skip this step if you downloaded the already demultiplexed files from GEO.**

## 4. Count genotypes

Count the genotypes in the library:

```bash
bash 004_count_genotypes.sh
```

This step produces a text file with a table like this, showing the number of times each genotype was found:

| Genotype                                                        | Occurrences | MutationCount |
|-----------------------------------------------------------------|-------------|---------------|
| GATCCAGATCTAACTTGGGGTGGCTTTGTCTTCTTCTTTTGCCAATTCCACTAATTGTTTGGG | 27336       | 0             |
| GATCCAGATCTAACTTGGTGTGGTTGTGTCTTCTGCTTTTGCCGATTCTAGTAATTGTTTGGG | 9559        | 7             |
| GATCCAGATCTAACTTGGGGTGGTTTTGTCTTCTGCTTTTGCCGATTCCAGTAATTGTTTGGG | 9514        | 4             |
| GATCCAGATCTAACTTGGGGTGGTTGTGTCTTCTGCTTTTGCCGATTCCAGTAATTGTTTGGG | 9377        | 5             |
| GATCCAGATCTAACTTGGGGTGGCTGTGTCTTCTGCTTTTGCCGATTCTAGTAATTGTTTGGG | 8960        | 5             |


## 5. Remove sequencing errors

Genotypes containing a mutation not included in the original library design (which allows for a maximum of 3072 different genotypes) was considered a sequencing error. We filtered these sequencing errors with the code in `005_remove_sequencing_errors`:

```bash
bash 005_remove_sequencing_errors
``` 
