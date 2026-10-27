#!/bin/bash

for gran in 40 20 10 5 2.5 1.25 0.625; do

    for group in short medium long; do
        mkdir -p random/gran_${gran}/${group}

        for i in $(seq 1 100); do
            bedtools shuffle \
                -i output_${gran}'_'${group}.bed \
                -g hg38_chr2.genome \
                -chrom \
                -seed $i \
                > random/gran_${gran}/${group}/random_${i}.bed
        done
    done
done