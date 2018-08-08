#!/bin/bash

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
done < Data/Sample_IDs.txt

# for each sample, remove columns I don't want AND say whether quality is overall good ('Pass') or bad ('Fail')
for i in `seq 1 12`; do
	# save sample name as a variable
	SAMPLE=$(head -n $i Data/Sample_IDs.txt | tail -n 1 | cut -f 2)
	
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

# take the location of all the trimmed files
for i in `seq 1 12`; do
	SAMPLE=$(head -n $i Sample_IDs.txt | tail -n 1 | cut -f 2)
	SAMPLE_PATH=${SAMPLE}/INCLUSION_LEVELS_TRIMMED.tab
	SAMPLES_ARRAY[i]=${SAMPLE_PATH}
done

# and paste them together
paste -d "\t" ${SAMPLES_ARRAY[@]} | cut -f 1,2,3,4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60 > SF3B1wt_SF3B1mut_TABLE.txt

# take only alternative exons
awk 'NR==1 {print}; $2 ~ /^HsaEX/ {print $0}' SF3B1wt_SF3B1mut_TABLE.txt >> Data/SF3B1wt_SF3B1mut_TABLE_EXONS.txt