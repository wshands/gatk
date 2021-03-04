#!/bin/bash

# make sure you have the most updated version
#cd ~/gatk
#./gradlew installDist
#cd -

~/gatk/gatk --java-options "-Xmx4g" \
    ExtractFeatures \
        --ref-version 38  \
        -R ~/gatk/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz  \
        -O acmg_feature_extract_debug.vcf.gz \
        --local-sort-max-records-in-ram 1000000 \
        --sample-table spec-ops-aou.kc_acmg_tieout_v6.metadata \
        --alt-allele-table  spec-ops-aou.kc_acmg_tieout_v6.alt_allele \
        --min-location 1000000000000 --max-location 1000035055462 \
        --project-id spec-ops-aou
