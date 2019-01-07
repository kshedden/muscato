#!/bin/bash

TARGET=/var/tmp/muscato_test
mkdir -p $TARGET

muscato_gendat -NumRead=10000000 -NumGene=10000000 -Dir=${TARGET}

muscato_prep_targets ${TARGET}/genes.txt.gz

muscato -GeneIdFileName=${TARGET}/genes_ids.txt.sz -GeneFileName=${TARGET}/genes.txt.sz \
        -ReadFileName=${TARGET}/reads.fastq -WindowWidth=20 -Windows=10,30,50,70 -MaxReadLength=200 \
        -TempDir=${TARGET}/muscato_tmp -NoCleanTemp -ResultsFileName=${TARGET}/results.txt
