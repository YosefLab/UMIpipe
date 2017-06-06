#!/bin/bash

#PBS  -N david1
#PBS -q yosef2
#PBS -l mem=2g
cd /data/yosef2/users/chenling/YosefCode/packages/UMI/Nir_David

python ../UMIpipe.py --fq1 FASTQs/Run_R2D2_462/S0003/Nir_d5_minus_3_S3_R1_001.fastq.gz --fq2 FASTQs/Run_R2D2_462/S0003/Nir_d5_minus_3_S3_R2_001.fastq.gz --ref hg38 --samplename S3 --ncells 500 --cend 8 --mstart 9 --mend 18

