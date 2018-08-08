#!/bin/bash

cut -f 1,2,46,81,305 GSM1911086_A5SS_spliced_reads.txt | awk -v min=100 '$2 + $3 + $4 + $5 >= min {print $0}' > A5SS_dataset.txt

mv GSM1911085_A5SS_seq.txt A5SS_sequences.txt