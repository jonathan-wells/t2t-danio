#!/usr/bin/env bash

    
genomefile="../../data/genome/fish11t2t.fasta"

############################################################################
## Run RepeatModeller
############################################################################

echo "RepeatModeler - processing T2T Danio rerio genome"

/programs/RepeatModeler-2.0.1/BuildDatabase \
    -name "../../data/repeatmodeler-out/fish11t2t_rmodeler_database" \
    $genomefile

/programs/RepeatModeler-2.0.1/RepeatModeler \
    -database "../../data/repeatmodeler-out/fish11t2t_rmodeler_database" \
    -LTRStruct \
    -pa 5

echo "RepeatModeler - Finished T2T Danio rerio"

