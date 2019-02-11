#!/bin/bash

# alternative 5' splice sites
awk 'NR==1 {print}; $2 ~ /^HsaALTD/ {print $0}' SF3B1wt_SF3B1mut_TABLE.txt >> Data/SF3B1wt_SF3B1mut_TABLE_ALTD.txt

# alternative 3' splice sites
awk 'NR==1 {print}; $2 ~ /^HsaALTA/ {print $0}' SF3B1wt_SF3B1mut_TABLE.txt >> Data/SF3B1wt_SF3B1mut_TABLE_ALTA.txt
