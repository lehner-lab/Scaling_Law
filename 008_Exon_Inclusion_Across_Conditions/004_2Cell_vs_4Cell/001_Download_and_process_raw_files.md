# Raw sequencing files

In this document, I explain the code found in [001a\_download\_raw\_files.sh](001a_download_raw_files.sh) and [001b\_process\_raw\_files.sh](001b_process_raw_files.sh). All the code in this document is written in bash. To run the code written in this document, you'll also need the 
`Samples_ID.txt` file found in the `Data` folder.


## 1. Source of data

To compare exon levels in the two- and four-cell states of the human embryo, we looked at the data generated by [Yan et al, 2013](https://www.nature.com/articles/nsmb.2660).


## 2. Download files

The raw RNA-seq datasets are stored in GEO. The easiest way to download GEO data is to use the [SRA toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/). For this datasets I already had the relevant SRA IDs, so I used this information to prefetch the files:

```bash
# prefetch SRR files from GEO
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=${myArray[0]}
	prefetch -v $SRR_ID
done < Data/Sample_IDs.txt
```
The next step involves extracting the fastq files. I extracted them to a folder called `fastq_files` in the current working directory, but feel free to replace this with something else. By default, the files prefetched in the previous step are stored in `~/ncbi/public/sra/`, but if you saved them elsewhere, then this directory should also be changed in the code below:

```bash
# make a directory to store the fastq files
mkdir fastq_files

# extract fastq files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=${myArray[0]}
	fastq-dump --outdir fastq_files/ --split-files ~/ncbi/public/sra/${SRR_ID}.sra
done < Data/Sample_IDs.txt
```
Finally, I changed the name of the fastq files to the sample name since that is more meaningful:

```bash
# rename files
while IFS=$'\t' read -r -a myArray
do
	SRR_ID=${myArray[0]}
	mv fastq_files/${SRR_ID}_1.fastq fastq_files/${myArray[1]}.fastq
done < Data/Sample_IDs.txt
```

## 3. Process files with VAST-TOOLS

To process the raw sequencing files and extract genomewide exon inclusion levels (as well as information about other alternative splicing events), I used [VAST-TOOLS](https://github.com/vastgroup/vast-tools). The first step in the VAST-TOOLS pipeline is the alignment. This can take quite a bit of time, which is why I recommend using multiple cores if possible. 


```bash
# run vast-tools
while IFS=$'\t' read -r -a myArray
do
	# make a directory for each sample and go inside
	mkdir ${myArray[1]}
	cd ${myArray[1]}
	
	# store the path to the fastq files we want to process
	FILE1=../fastq_files/${myArray[1]}.fastq
	
	# run vast-tools align
	vast-tools align $FILE1 --sp Hsa -c 10   # using 10 cores
		
	# move back up before creating a directory for the next sample
	cd ..
done < Data/Sample_IDs.txt
```
These samples come from single-cell RNA-seq experiments which were [not sequenced deep enough](https://github.com/vastgroup/vast-tools#merging-outputs) for VAST-TOOLS to work properly. The workaround proposed by the VAST-TOOLS team is to use `vast-tools merge` and merge the information from all the different single-cell samples into one file before running `vast-tools combine`:

```bash
EMBRYOS=(2C_Embryo1 2C_Embryo2 2C_Embryo3 4C_Embryo1 4C_Embryo2 4C_Embryo3)

for i in ${EMBRYOS[@]}; do
	GROUP_FILE=Data/${i}_Samples_to_Merge.txt
	
	# create a directory where we'll put all the samples to be merged
	mkdir -p ${i}/vast_out/to_combine
	
	# move relevant files to this new directory
	while IFS=$'\t' read -r -a myArray
	do
		ORIGINAL_SAMPLE=${myArray[0]}
		mv ${ORIGINAL_SAMPLE}/vast_out/to_combine/* ${i}/vast_out/to_combine/
	done < $GROUP_FILE
	
	# move into the merged sample directory and run vast-tools merge
	cd ${i}
	vast-tools merge --groups $GROUP_FILE --outDir vast_out
	
	# run vast-tools combine
	vast-tools combine -o vast_out -sp Hsa
	
	# move back up before creating a directory for the next embryo
	cd ..
done
```


## Build table to use in downstream analysis

The output of `vast-tools combine` has [many columns](https://github.com/vastgroup/vast-tools#combine-output-format). I removed some of them which I was not interested in, and I also changed the quality scores column. If the quality of an event was `LOW` or better, I set the event as a `Pass`. Otherwise, I called it a `Fail`.

```bash
EMBRYOS=(2C_Embryo1 2C_Embryo2 2C_Embryo3 4C_Embryo1 4C_Embryo2 4C_Embryo3)

# for each embryo, remove columns I don't want AND say whether quality is overall good ('Pass') or bad ('Fail')
for i in ${EMBRYOS[@]}; do
	# move into this embryo's directory
	cd ${i}
	
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
EMBRYOS=(2C_Embryo1 2C_Embryo2 2C_Embryo3 4C_Embryo1 4C_Embryo2 4C_Embryo3)

# take the location of all the trimmed files
for i in `seq 1 6`; do
	SAMPLE_NAME=${EMBRYOS[${i}]}
	SAMPLE_PATH=${SAMPLE_NAME}/INCLUSION_LEVELS_TRIMMED.tab
	SAMPLES_ARRAY[${i}]=${SAMPLE_PATH}
done

# and paste them together
paste -d "\t" ${SAMPLES_ARRAY[@]} | cut -f 1,2,3,4,5,9,10,14,15,19,20,24,25,29,30 > 2C_4C_TABLE.txt
```
Finally, since I was only interested in alternative exon events, I removed all other events for downstream analyses:

```bash
awk 'NR==1 {print}; $2 ~ /^HsaEX/ {print $0}' 2C_4C_TABLE.txt >> Data/2C_4C_TABLE_EXONS.txt
```
The `2C_4C_TABLE_EXONS.txt` file produced is provided in the `Data` folder.