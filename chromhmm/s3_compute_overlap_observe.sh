#!/bin/bash

for gran in 40 20 10 5 2.5 1.25 0.625; do

    mkdir -p overlap/observed/gran_${gran}

    for group in short medium long; do
        bedtools intersect \
            -a output_${gran}'_'${group}.bed \
            -b E038_15_coreMarks_hg38lift_mnemonics.bed.gz \
            -wao \
            > overlap/observed/gran_${gran}/${group}_E038.tsv
    done
done