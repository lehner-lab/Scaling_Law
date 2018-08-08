#!/bin/bash

# alternative 5' splice sites
awk 'NR==1 {print}; $2 ~ /^HsaALTD / {print $0}' 2C_4C_TABLE.txt >> Data/2C_4C_TABLE_ALTD.txt

# alternative 3' splice sites
awk 'NR==1 {print}; $2 ~ /^HsaALTA / {print $0}' 2C_4C_TABLE.txt >> Data/2C_4C_TABLE_ALTA.txt
