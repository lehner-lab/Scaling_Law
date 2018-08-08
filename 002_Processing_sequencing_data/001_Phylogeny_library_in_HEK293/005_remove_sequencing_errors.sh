#!/bin/bash

for i in `seq 1 9`;
do
	grep GATCCAGATCTAACTTG[CGT][GT]GTGG[CT]T[GT]TGTCT[CT]CT[GT]CTT[CT]T[GC]CC[AG]ATTC[CT]A[GC]TAATTGTTTGGG Phylogeny_HEK293_BR_Rep_${i}.counts > Phylogeny_HEK293_BR_Rep_${i}.3072.counts
done

for i in `seq 1 3`;
do
	grep GATCCAGATCTAACTTG[CGT][GT]GTGG[CT]T[GT]TGTCT[CT]CT[GT]CTT[CT]T[GC]CC[AG]ATTC[CT]A[GC]TAATTGTTTGGG Phylogeny_HEK293_TR_Rep_${i}.counts > Phylogeny_HEK293_TR_Rep_${i}.3072.counts
done