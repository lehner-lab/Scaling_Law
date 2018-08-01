# sQTL effects per tissue

In this document, I explain the code found in [008a\_all\_sQTL\_effects\_per\_tissue.sh](008a_all_sQTL_effects_per_tissue.sh), which goes through each of the folders inside the `sQTL_tests/` directory and calculates the effect of the given sQTL in all the different tissues for which we have data. All the code in this document is written in bash.


## 1. Slow code

First, move to the `sQTL_tests/` directory 

```bash
# move to directory where I'm keeping all the information
cd sQTL_tests
```
There are 35153 exon skipping events in Gene_Regions.bed that we are studying, so we'll run a for loop where the code will iterate through each of these alternative splicing events, and work with each one of them, one by one:

```bash
# iterate through the 35153 exon skipping events
for this_number in `seq 1 35153`;
do
    # code here
done
```
Inside the for loop, we first move into the folder corresponding to the given alternative exon:

```bash
# move to directory corresponding to this splicing event
cd ${this_number}
```
Next, we exract the ID of the exon skipping event:

```bash
# extract the splicing event ID from Gene_Regions.bed
this_splicing_event_id=$(cat Gene_Regions.bed | cut -f 4)
```
And use this ID to get the line numbers where this splicing event shows up inside `007_Significant_Splicing_Events_sQTLs_Combos.txt` (generated earlier):

```bash
# get the line numbers where this event appears inside
# 007_Significant_Splicing_Events_sQTLs_Combos.txt
line_numbers=$(grep -n $this_splicing_event_id ../../Data/007_Significant_Splicing_Events_sQTLs_Combos.txt | cut -f 1 -d :)
```
Use the information about the line numbers to find out in how many lines this alternative splicing event is found:

```bash
# count how many lines contained this splicing event
number_of_lines=$(echo $line_numbers | awk '{print NF}')
```
This is useful because if the splicing event was not found in `007_Significant_Splicing_Events_sQTLs_Combos.txt`, that means that no sQTL was found associated with it, and we will not be considering this splicing event any longer. Therefore, the code now has an 'if' statement that checks whether any sQTLs were found for this splicing event:

```bash
# make sure there are some significant sQTL's here...
if (( $number_of_lines > 0 )); then
	
	# iterate through each of the sQTL's
	for i in $line_numbers; do
		
		# get sQTL ID
		sQTL_id=$(head -n $i ../Data/011_ALT_Significant_Splicing_Events_sQTLs_Combos.txt | tail -n 1 | cut -f 2)
		
		# run the R script
		Rscript ../../008b_sQTL_effects_per_tissue.R $sQTL_id
		
	done
fi
```
As you can see, inside the 'if' statement, we go through each of the sQTLs significantly associated with our given splicing event, and run [008b\_sQTL\_effects\_per\_tissue.R](008b_sQTL_effects_per_tissue.R) to calculate what this effect is. Finally, we move back up before moving to the next iteration of the for loop:

```bash
cd ..
```


## 2. Faster code

Alternatively, if you have access to a cluster with an SGE batch system, you can simply submit the following code to the queue. A separate job will be sent for each of the 35153 exon skipping events, which means the whole thing will take far less time:

```bash
#!/bin/bash
#$ -t 1-35153
#$ -o location_of_queue_output_files
#$ -e location_of_queue_error_files

# save SGE_TASK_ID variable with more user-friendly name
this_number=${SGE_TASK_ID}

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
		sQTL_id=$(head -n $i ../Data/011_ALT_Significant_Splicing_Events_sQTLs_Combos.txt | tail -n 1 | cut -f 2)
		
		# run the R script
		Rscript ../../008b_sQTL_effects_per_tissue.R $sQTL_id
		
	done
fi
```