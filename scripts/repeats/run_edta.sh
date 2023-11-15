#!/usr/bin/env bash

singularity exec -C --bind $PWD --pwd $PWD /programs/edta-2.0.0/EDTA.sif EDTA.pl \
    --genome fish11t2t.fasta \
    --curatedlib zebrep.ref \
    --cds cds_nodup.fa \
    --species others \
    --step all \
    --overwrite 1 \
    --evaluate 1 \
    --anno 1 \
    --sensitive 0 \
    --threads 16
