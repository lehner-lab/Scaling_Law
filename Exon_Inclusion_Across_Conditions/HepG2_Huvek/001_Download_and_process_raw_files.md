# HepG2 vs Huvek


To compare exon inclusion in two different cell lines, we chose two cell lines representing different lineages (HepG2 and Huvek) from the [ENCODE Project Common Cell Types](https://www.genome.gov/26524238/encode-project-common-cell-types/). Transcriptomic RNA-seq data from these cell lines is stored in GEO. To run the code written in this document, you'll also need the 
`Samples_ID.txt` file found in this folder.

## Download files

The easiest way to download GEO data is to use the [SRA toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/). First, we will build a table linking the GSM IDs we have to the SRA IDs required by the SRA toolkit:

```bash
# Make SRR-to-GSM table
wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab
grep ^SRR SRA_Accessions.tab | grep GSM | awk 'BEGIN {FS="\t"; print "GSM" FS "SRR"}; {print $10 FS $1}' > GSM_SRR.txt
rm SRA_Accessions.tab
```
Once we have this table, we can use it to extract the sample SRA IDs and prefetch the files:

```bash
# prefetch SRR files from GEO
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=$(grep ${myArray[0]} GSM_SRR.txt | cut -f 2)
	prefetch -v $SRR_ID
done < Sample_IDs.txt
```
And now you need to extract the fastq files. I extracted them to a folder called `fastq_files` in the current working directory, but feel free to replace this with something else. By default, the files prefetched in the previous step are stored in `~/ncbi/public/sra/`, but if you saved them elsewhere, then this directory should also be changed in the code below:

```bash
# extract fastq files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=$(grep ${myArray[0]} GSM_SRR.txt | cut -f 2)
	fastq-dump --outdir fastq_files/ --split-files ~/ncbi/public/sra/${SRR_ID}.sra
done < Sample_IDs.txt
```
Finally, you can change the name of the fastq files to something more meaningful such as the sample name:

```bash
# rename files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=$(grep ${myArray[0]} GSM_SRR.txt | cut -f 2)
	mv fastq_files/${SRR_ID}_1.fastq fastq_files/${myArray[1]}_1.fastq
	mv fastq_files/${SRR_ID}_2.fastq fastq_files/${myArray[1]}_2.fastq
done < Sample_IDs.txt
```

## Process files with VAST-TOOLS

To process the raw sequencing files and extract genomewide exon inclusion levels (as well as information about other alternative splicing events), I used [VAST-TOOLS](https://github.com/vastgroup/vast-tools). The first step in the VAST-TOOLS pipeline is the alignment. This can take quite a bit of time, which is why I recommend using multiple cores if possible. The second step is `combine` and should be fairly quick.

```bash
# run vast-tools
while IFS=$'\t' read -r -a myArray
do
	# make a directory for each sample and go inside
	mkdir ${myArray[1]}
	cd ${myArray[1]}
	
	# store the path to the fastq files we want to process
	FILE1=../fastq_files/${myArray[1]}_1.fastq
	FILE2=../fastq_files/${myArray[1]}_2.fastq
	
	# run vast-tools align
	vast-tools align $FILE1 $FILE2 --sp Hsa -c 10   # using 10 cores
	# run vast-tools combine
	vast-tools combine -o vast_out -sp Hsa
	
	# move back up before creating a directory for the next sample
	cd ..
done < Sample_IDs.txt
```



## Build table to use in downstream analysis

The output of `vast-tools combine` has [many columns](https://github.com/vastgroup/vast-tools#combine-output-format). I removed some of them which I was not interested in, and I also changed the quality scores column. If the quality of an event was `LOW` or better, I set the event as a `Pass`. Otherwise, I called it a `Fail`.

```bash
# for each sample, remove columns I don't want AND say whether quality is overall good ('Pass') or bad ('Fail')
for i in `seq 1 4`; do
	# save sample name as a variable
	SAMPLE=$(head -n $i Sample_IDs.txt | tail -n 1)
	
	# move into this sample's directory
	cd $SAMPLE
	
	# select columns I want to keep AND change the value of the quality score column to either 'Pass' or 'Fail'
	cut -f 1,2,4,7,8 vast_out/INCLUSION_LEVELS_FULL-Hsa1-hg19.tab | head -n 1 | awk '{print $0}' > INCLUSION_LEVELS_TRIMMED.tab
	tail -n +2 vast_out/INCLUSION_LEVELS_FULL-Hsa1-hg19.tab | awk 'BEGIN{OFS="\t"}{
		if ($8 ~ /^(SOK|OK|LOW),(SOK|OK|LOW),(SOK|OK|LOW|NA),/ )
			print $1 OFS $2 OFS $4 OFS $7 OFS "Pass";
		else
			print $1 OFS $2 OFS $4 OFS $7 OFS "Fail";
		}' >> INCLUSION_LEVELS_TRIMMED.tab
	
	# move back up before working on the next sample
	cd ..
done
```
After trimming (removing unwanted columns) the files and modifying the quality score column, I merged all the tables into one:

```bash
# take the location of all the trimmed files
for i in `seq 1 4`; do
	SAMPLE=$(head -n $i Sample_IDs.txt | tail -n 1)
	SAMPLE_PATH=${SAMPLE}/INCLUSION_LEVELS_TRIMMED.tab
	SAMPLES_ARRAY[i]=${SAMPLE_PATH}
done

# and paste them together
paste -d "\t" ${SAMPLES_ARRAY[@]} | cut -f 1,2,3,4,5,9,10,14,15,19,20 > Huvec_HepG2_TABLE.txt
```
Finally, since I was only interested in alternative exon events, I removed all other events for downstream analyses:


```bash
awk 'NR==1 {print}; $2 ~ /^HsaEX/ {print $0}' Huvec_HepG2_TABLE.txt >> Huvec_HepG2_TABLE_EXONS.txt
```