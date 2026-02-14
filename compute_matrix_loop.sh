#!/bin/bash

for i in '0.1' '0.2' '0.5' '1' '2' '5'
do
	for j in GSM5652176 GSM5652296 GSM5652297
	do
		python compute_breakpoints.py -i roadmap/$j''_WGBS_Seg_L$i''.bigwig -c 2 -o matrix/$j''_WGBS_Seg_L$i''.csv
	done
done

for i in '0.1' '0.2' '0.5' '1' '2' '5'
do
	for j in GSM5652176 GSM5652296 GSM5652297
	do
		Rscript cluster_breakpoints.R matrix/$j''_WGBS_Seg_L$i''.csv matrix/$j''_WGBS_Seg_L$i
	done
done
