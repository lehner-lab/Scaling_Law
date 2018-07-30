#!/bin/bash

cut -f1,6 GTEx_v7_Annotations_SampleAttributesDS.txt | sed s/-/\./g > IDs_Tissues.txt

gunzip GTEx_Analysis_2016-01-15_v7_STARv2.4.2a_junctions.gct.gz