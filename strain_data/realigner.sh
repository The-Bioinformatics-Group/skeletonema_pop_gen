#!/bin/bash

#Identify regions in need of realignment:

java -Xmx2g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
  -T RealignerTargetCreator \
  -R ./long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta \
  -o merged_output.intervals \
  -I merged.bam 
  
#  defaults for optional parameters:
#  --minReadsAtLocus N [the minimum coverage at a locus for the entropy calculation to be enabled; default=4]
#  --windowSize N [any two SNP calls and/or high entropy positions are considered clustered when they occur no more than N basepairs apart; default=10]
#  --mismatchFraction f [fraction of total sum of base qualities at a position that need to mismatch for the position to be considered to have high entropy; default=0.15; to disable, set to <= 0 or > 1]
#  Note that this fraction should be adjusted based on your particular data set. For deep coverage and/or when looking for indels with low allele frequency, this number should be smaller.
#  --maxIntervalSize [max size in bp of intervals that we'll pass to the realigner; default=500]

#-------------------------------------------------------------------------

#  Run realigner over intervals:

java -Xmx4g -jar ../scripts/GenomeAnalysisTK_2_5.jar \
  -I merged.bam \
  -R ./long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta \
  -T IndelRealigner \
  -targetIntervals merged_output.intervals \
  -o merged_realigned.bam 
