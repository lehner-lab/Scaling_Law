#!/bin/bash

# store the absolute location of two files we'll use in this script
GENE_REGIONS_BED=absolute_path_to_Gene_Regions.bed_file_generated_earlier
PRUNED_GENOTYPES_FILE=absolute_path_to_genotypes.MAF01.vcf.gz_file_generated_earlier

# move to directory where I'll keep all the information
mkdir sQTL_tests
cd sQTL_tests

# iterate through the 35153 exon skipping events
for this_number in `seq 1 35153`;
do
    # make a new directory and move to it
    mkdir $this_number
    cd $this_number
    
    # subset Gene_Regions.bed file to only include the gene we're interested in right now
    tail -n +$(( this_number + 1 )) $GENE_REGIONS_BED  | head -n 1 > ./Gene_Regions.bed
    
    # now subset the vcf file and extract it
    bcftools view --regions-file ./Gene_Regions.bed $PRUNED_GENOTYPES_FILE -O z -o ./genotypes.MAF01.GeneRegions.vcf.gz
    gunzip ./genotypes.MAF01.GeneRegions.vcf.gz
    
    # process the vcf file
    awk -f ../../006b_process_vcf.awk ./genotypes.MAF01.GeneRegions.vcf > ./genotypes.MAF01.GeneRegions.Parsed.txt
    
    # calculate p values for all possible qtls
    Rscript ../../006c_sQTL_test.R
    
    # move back up
    cd ..

done



# Code to send to an SGE queue
# (the sQTL_tests folder must already be created before submitting the job)
############################################################################################################################
##                                                                                                                        ##
## #!/bin/bash                                                                                                            ##
## #$ -t 1-35153                                                                                                          ##
## #$ -o location_of_queue_output_files                                                                                   ##
## #$ -e location_of_queue_error_files                                                                                    ##
##                                                                                                                        ##
## # save SGE_TASK_ID variable with more user-friendly name                                                               ##
## this_number=${SGE_TASK_ID}                                                                                             ##
##                                                                                                                        ##
## # store the absolute location of two files we'll use in this script                                                    ##
## GENE_REGIONS_BED=absolute_path_to_Gene_Regions.bed_file_generated_earlier                                              ##
## PRUNED_GENOTYPES_FILE=absolute_path_to_0.2.Pruned.genotypes.MAF01.vcf.gz_file_generated_earlier                        ##
##                                                                                                                        ##
## # move to sQTL_tests  directory where I'll keep all the information                                                    ##
## # (must have been created before submitting this job to the queue)                                                     ##
## cd absolute_path_to_sQTL_tests                                                                                         ##
##                                                                                                                        ##
## # make a new directory and move to it                                                                                  ##
## mkdir $this_number                                                                                                     ##
## cd $this_number                                                                                                        ##
##                                                                                                                        ##
## # subset Gene_Regions.bed file to only include the gene we're interested in right now                                  ##
## tail -n +$(( this_number + 1 )) $GENE_REGIONS_BED  | head -n 1 > ./Gene_Regions.bed                                    ##
##                                                                                                                        ##
## # now subset the vcf file and extract it                                                                               ##
## bcftools view --regions-file ./Gene_Regions.bed $PRUNED_GENOTYPES_FILE -O z -o ./genotypes.MAF01.GeneRegions.vcf.gz    ##
## gunzip ./genotypes.MAF01.GeneRegions.vcf.gz                                                                            ##
##                                                                                                                        ##
## # process the vcf file                                                                                                 ##
## awk -f ../../006b_process_vcf.awk ./genotypes.MAF01.GeneRegions.vcf > ./genotypes.MAF01.GeneRegions.Parsed.txt         ##
##                                                                                                                        ##
## # calculate p values for all possible qtls                                                                             ##
## Rscript ../../006c_sQTL_test.R                                                                                         ##
##                                                                                                                        ##
############################################################################################################################

