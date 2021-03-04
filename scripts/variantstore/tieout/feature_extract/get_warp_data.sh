#!/bin/bash

# Get a shard of pre-hard filtered extracted data from the WARP run:

WORKFLOW_ID=68026724-7fdc-4a0e-bcfa-6a5e4a86cc0a
gsutil cp gs://broad-dsp-spec-ops-cromwell-execution/JointGenotyping/${WORKFLOW_ID}/call-TotallyRadicalGatherVcfs/shard-0/*.gnarly.vcf.gz bq_validation_v5_35.0.gnarly.vcf.gz
#gunzip bq_validation_v5_35.0.gnarly.vcf.gz
