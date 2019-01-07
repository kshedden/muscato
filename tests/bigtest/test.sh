#!/bin/bash

TARGET=/var/tmp/muscato_test
mkdir -p $TARGET

muscato_gendat -NumRead=100000 -NumGene=100000 -Dir=${TARGET}

echo "Starting prep_targets"
muscato_prep_targets ${TARGET}/genes.txt.sz

muscato -GeneIdFileName=${TARGET}/musc_ids_genes.txt.sz -GeneFileName=${TARGET}/musc_genes.txt.sz \
        -ReadFileName=${TARGET}/reads.fastq -WindowWidth=20 -Windows=10,30,50,70 -MaxReadLength=200 \
        -TempDir=${TARGET}/muscato_tmp -NoCleanTemp -ResultsFileName=${TARGET}/results.txt

