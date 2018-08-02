#!/bin/bash

# move to directory where I'm keeping all the information
cd sQTL_tests

# iterate through the 35153 exon skipping events
for this_number in `seq 1 35153`;
do
    # move to directory corresponding to this splicing event
    cd ${this_number}
    
    # extract the splicing event ID from Gene_Regions.bed
    this_splicing_event_id=$(cat Gene_Regions.bed | cut -f 4)
    
    # get the line numbers where this event appears inside
    # 007_Significant_Splicing_Events_sQTLs_Combos.txt
    line_numbers=$(grep -n $this_splicing_event_id ../../Data/007_Significant_Splicing_Events_sQTLs_Combos.txt | cut -f 1 -d :)
    
    # count how many lines contained this splicing event
    number_of_lines=$(echo $line_numbers | awk '{print NF}')
    
    # make sure there are some significant sQTL's here...
	if (( $number_of_lines > 0 )); then
	
		# iterate through each of the sQTL's
		for i in $line_numbers; do
		
			# get sQTL ID
			sQTL_id=$(head -n $i ../../Data/007_Significant_Splicing_Events_sQTLs_Combos.txt | tail -n 1 | cut -f 2)
		
			# run the R script
			Rscript ../../008b_sQTL_effects_per_tissue.R $sQTL_id
		
		done
	fi
    
    # move back up
    cd ..

done




# Code to send to an SGE queue
##################################################################################################################################
##                                                                                                                              ##
## #!/bin/bash                                                                                                                  ##
## #$ -t 1-35153                                                                                                                ##
## #$ -o location_of_queue_output_files                                                                                         ##
## #$ -e location_of_queue_error_files                                                                                          ##
##                                                                                                                              ##
## # save SGE_TASK_ID variable with more user-friendly name                                                                     ##
## this_number=${SGE_TASK_ID}                                                                                                   ##
##                                                                                                                              ##
## # move to directory corresponding to this splicing event                                                                     ##
## cd ${this_number}                                                                                                            ##
##                                                                                                                              ##
## # extract the splicing event ID from Gene_Regions.bed                                                                        ##
## this_splicing_event_id=$(cat Gene_Regions.bed | cut -f 4)                                                                    ##
##                                                                                                                              ##
## # get the line numbers where this event appears inside                                                                       ##
## # 007_Significant_Splicing_Events_sQTLs_Combos.txt                                                                           ##
## line_numbers=$(grep -n $this_splicing_event_id ../../Data/007_Significant_Splicing_Events_sQTLs_Combos.txt | cut -f 1 -d :)  ##
##                                                                                                                              ##
## # count how many lines contained this splicing event                                                                         ##
## number_of_lines=$(echo $line_numbers | awk '{print NF}')                                                                     ##
##                                                                                                                              ##
## # make sure there are some significant sQTL's here...                                                                        ##
## if (( $number_of_lines > 0 )); then                                                                                          ##
## 	                                                                                                                            ##
## 	# iterate through each of the sQTL's                                                                                        ##
## 	for i in $line_numbers; do                                                                                                  ##
## 		                                                                                                                        ##
## 		# get sQTL ID                                                                                                           ##
## 		sQTL_id=$(head -n $i ../../Data/007_Significant_Splicing_Events_sQTLs_Combos.txt | tail -n 1 | cut -f 2)                ##
## 		                                                                                                                        ##
## 		# run the R script                                                                                                      ##
## 		Rscript ../../008b_sQTL_effects_per_tissue.R $sQTL_id                                                                   ##
##                                                                                                                              ##
## 	done                                                                                                                        ##
## fi                                                                                                                           ##
##                                                                                                                              ##
##################################################################################################################################