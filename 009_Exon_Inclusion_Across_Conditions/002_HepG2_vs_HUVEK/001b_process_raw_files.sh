#!/bin/bash

# for each sample, remove columns I don't want AND say whether quality is overall good ('Pass') or bad ('Fail')
for i in `seq 1 4`; do
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
for i in `seq 1 4`; do
	SAMPLE=$(head -n $i Sample_IDs.txt | tail -n 1 | cut -f 2)
	SAMPLE_PATH=${SAMPLE}/INCLUSION_LEVELS_TRIMMED.tab
	SAMPLES_ARRAY[i]=${SAMPLE_PATH}
done

# and paste them together
paste -d "\t" ${SAMPLES_ARRAY[@]} | cut -f 1,2,3,4,5,9,10,14,15,19,20 > Huvec_HepG2_TABLE.txt

# take only alternative exons
awk 'NR==1 {print}; $2 ~ /^HsaEX/ {print $0}' Huvec_HepG2_TABLE.txt >> Data/Huvec_HepG2_TABLE_EXONS.txt