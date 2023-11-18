#!/usr/bin/env bash

tename=$1
rg "\t${tename}\t" ../../data/edta-out/fish11t2t.fasta.mod.out.bed > "${tename}_1.bed"
bedtools slop \
    -i "${tename}_1.bed" \
    -g ../../data/genome/fish11t2t.genome \
    -b 50 > "${tename}_2.bed"
bedtools getfasta \
    -fi ../../data/genome/fish11t2t.fasta \
    -bed "${tename}_2.bed" \
    -s > "../../data/seqs/${tename}.fa"
mafft --thread 16 "../../data/seqs/${tename}.fa" > "../../data/seqs/${tename}.mafft.fa"
rm ${tename}_{1,2}.bed

