#!/bin/bash

for gran in 40 20 10 5 2.5 1.25 0.625; do
    for group in short medium long; do

        mkdir -p overlap/random/gran_${gran}/${group}

        for i in $(seq 1 100); do
            bedtools intersect \
                -a random/gran_${gran}/${group}/random_${i}.bed \
                -b E038_15_coreMarks_hg38lift_mnemonics.bed.gz \
                -wao \
                > overlap/random/gran_${gran}/${group}/random_${i}_E038.tsv
        done
    done
done