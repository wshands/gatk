#!/bin/bash

cd ~/gatk
./gradlew installDist
cd -

#
#WORKFLOW_ID=68026724-7fdc-4a0e-bcfa-6a5e4a86cc0a
#gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-ImportGVCFs/shard-0/attempt-2/genomicsdb.tar .
#
#tar -xf genomicsdb.tar
WORKSPACE=genomicsdb

reference="/Users/marymorg/gatk/src/test/resources/large/Homo_sapiens_assembly38.fasta.gz"
#reference="/Users/kcibul/projects/references/hg38/v0/Homo_sapiens_assembly38.fasta"

~/gatk/gatk --java-options "-Xms8g -Xdebug -Xrunjdwp:transport=dt_socket,address=5005,server=y,suspend=n" \
  GnarlyGenotyper \
  -R $reference \
  -O debug_warp_10492.vcf \
  -V gendb://$WORKSPACE \
  --only-output-calls-starting-in-intervals \
  -stand-call-conf 10 \
  -L chr1:10492
