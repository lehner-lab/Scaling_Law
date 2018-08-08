#!/bin/bash

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

# vast tools merge and combine
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

# take the location of all the trimmed files
for i in `seq 1 6`; do
	SAMPLE_NAME=${EMBRYOS[${i}]}
	SAMPLE_PATH=${SAMPLE_NAME}/INCLUSION_LEVELS_TRIMMED.tab
	SAMPLES_ARRAY[${i}]=${SAMPLE_PATH}
done

# and paste them together
paste -d "\t" ${SAMPLES_ARRAY[@]} | cut -f 1,2,3,4,5,9,10,14,15,19,20,24,25,29,30 > 2C_4C_TABLE.txt

# take only alternative exons
awk 'NR==1 {print}; $2 ~ /^HsaEX/ {print $0}' 2C_4C_TABLE.txt >> Data/2C_4C_TABLE_EXONS.txt