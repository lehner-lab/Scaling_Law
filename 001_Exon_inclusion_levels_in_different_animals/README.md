# Pipeline

This document describes the method used to estimate the PSI of exons in different animal genomes. As an example, I'll describe the process involved in extracting this information in chimp heart, but the process should be the same for other any other species/tissue. To run this pipeline, you will need:

* **Genome annotations** in .gtf format (the links to download all the annotations used in this analysis are found in [gtf_downloads.txt](./gtf_downloads.txt))
* **Genome sequence** in .fasta format (the links to download all the genomes used in this analysis are found in [genome_downloads.txt](./genome_downloads.txt))
* **Sample files** (the links to download all sample files are found in [sample_downloads.txt](./sample_downloads.txt))
* **STAR v2.5.2a** installed (downloadable from [here](https://github.com/alexdobin/STAR/archive/2.5.2a.tar.gz); the latest version of the software can be downloaded from the [STAR Github repo](https://github.com/alexdobin/STAR/releases) and should also work).
* **SAMtools v1.3.1** installed (downloadable from [here](https://sourceforge.net/projects/samtools/files/samtools/1.3.1/), later versions should also work).

Unless stated otherwise, all code in this markdown file is written in bash.


## 1. Prep work

Before calculating exon PSI values, we need to build the species genome indices and filter the genome annotations. This preparatory work only needs to be done once per species (e.g. in the case of Chimp, the files generated in this step will work for the heart, liver, lung, kidney and lymph node samples). All the code from this section can also be found in [example\_prep\_work.sh](./example_prep_work.sh).

### 1.1. Build STAR indices

We first download the genome file

```bash
# downloade chimp genome sequence
wget ftp://ftp.ensembl.org/pub/release-84/fasta/pan_troglodytes/dna/Pan_troglodytes.CHIMP2.1.4.dna.toplevel.fa.gz

# extract file
gunzip Pan_troglodytes.CHIMP2.1.4.dna.toplevel.fa.gz
```

and the annotations:

```bash
# download chimp genome annotations
wget ftp://ftp.ensembl.org/pub/release-84/gtf/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.84.gtf.gz

# extract the file
gunzip Pan_troglodytes.CHIMP2.1.4.84.gtf.gz
```
And build the STAR indices. This can take quite a bit of time and will use up a lot of RAM so feel free to change the `runThreadN` and `limitGenomeGenerateRAM` parameters to whatever your computer/cluster can handle:

```bash
# build folder to dump STAR indices
mkdir genomeDir

# build indices
STAR --runThreadN 5 --limitGenomeGenerateRAM 80000000000 --runMode genomeGenerate --genomeDir ./genomeDir --genomeFastaFiles Pan_troglodytes.CHIMP2.1.4.dna.toplevel.fa --sjdbGTFfile Pan_troglodytes.CHIMP2.1.4.84.gtf
```


### 1.2. Filter GTF file

We will subset the .gtf annotations so that we're only left with information relative to exons. For this, we can run the `dexseq_prepare_annotation.py` script from the [DEXseq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) package:

```bash
# subset annotations file and save it as "Exonic_Parts.gtf"
python dexseq_prepare_annotation.py Pan_troglodytes.CHIMP2.1.4.84.gtf Exonic_Parts.gtf
```



## 2. Calculate exon PSI levels

All the code from this section can also be found in [example\_calculate\_PSI.sh](./example_calculate_PSI.sh).


### 2.1. Download sample files

Download and extract the files from the [NHPRTR server](http://www.nhprtr.org/data/2014_NHP_tissuespecific/10NHP_Illumina_Baylor/):

```bash
# download chim heart RNA-seq files
wget http://www.nhprtr.org/data/2014_NHP_tissuespecific/10NHP_Illumina_Baylor/Chimpanzee/Chimpanzee_Heart_359_D2EKWACXX-8-ID22_1_sequence.txt.bz2
wget http://www.nhprtr.org/data/2014_NHP_tissuespecific/10NHP_Illumina_Baylor/Chimpanzee/Chimpanzee_Heart_359_D2EKWACXX-8-ID22_2_sequence.txt.bz2

# extract the files
bunzip2 Chimpanzee_Heart_359_D2EKWACXX-8-ID22_1_sequence.txt.bz2 Chimpanzee_Heart_359_D2EKWACXX-8-ID22_2_sequence.txt.bz2
```

### 2.2. Align RNA-seq reads to the reference genome

Align the sample reads to the reference genome. This step can also take a long time so feel free to change the `runThreadN` and parameter to whatever your computer/cluster can handle:

```bash
# build a directory to store STAR output
mkdir STAR_Output

# run STAR
STAR --runThreadN 10  --genomeDir ./genomeDir --readFilesIn Chimpanzee_Heart_359_D2EKWACXX-8-ID22_1_sequence.txt Chimpanzee_Heart_359_D2EKWACXX-8-ID22_2_sequence.txt --outFileNamePrefix STAR_Output/Pan_troglodytes.CHIMP2.1.4.84_Heart_
```

### 2.3. Change format of STAR output

The script used to calculate PSI values requires a TopHat-like input format, so we'll output the STAR output so it looks more like the output from TopHat:

```bash
# convert the output from STAR into a format
# that is similar to the output from TopHat
# oneliner code taken from this protocol:
# http://onlinelibrary.wiley.com/doi/10.1002/0471142905.hg1116s87/abstract
cd STAR_Output
awk 'BEGIN{OFS="\t"}{print $1, $2-20-1, $3+20, "JUNCBJ"NR, $7, ($4 == 1)? "+":"-",$2-20-1, $3+20, "255,0,0", 2, "20,20", "0,300" }' *SJ.out.tab > junctions.bed
```
Next, make a .bam file from the .sam file produced by STAR:

```bash
# convert STAR .sam output into .bam
samtools view -hb *Aligned.out.sam > accepted_hits.bam
```
Move files to a folder where we'll store all the files needed to calculate PSI

```bash
# go back up
mv ..

# make a directory
mkdir PSI_Calculation

# move some files...
mv STAR_Output/accepted_hits.bam PSI_Calculation/
mv STAR_Output/junctions.bed PSI_Calculation/

```


### 2.4. Calculate exon PSI values

Calculate exon PSI values using `PSI.sh`, the "compact pipeline" from [Schafer et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/26439713):

```bash
# move to relevant folder
cd PSI_Calculation

# start counting
bash ../PSI.sh StartPSI ../Exonic_Parts.gtf 101 accepted_hits.bam junctions.bed Chimp_Heart
```


### 2.5. FAS exon 6

The PSI of FAS exon 6 can be accessed by using the FAS gene identifier for each species and looking for the 63-nucleotide long exon (all identifiers in [fas_identifiers.txt](./fas_identifiers.txt))

