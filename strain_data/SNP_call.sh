#!/bin/bash

java -Xmx16g -jar ~/scripts/GenomeAnalysisTK_3_6.jar \
   -T HaplotypeCaller \
   -R ./long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta \
   -I merged_realigned.bam \
   -stand_call_conf 30.0 \
   -stand_emit_conf 10.0 \
   -L exons.bed \
   -o raw_snps_indels_Q30_exons.vcf

java -Xmx16g -jar ~/scripts/GenomeAnalysisTK_3_6.jar \
   -T HaplotypeCaller \
   -R ./long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta \
   -I merged_realigned.bam \
   -stand_call_conf 30.0 \
   -stand_emit_conf 10.0 \
   -L introns.bed \
   -o raw_snps_indels_Q30_introns.vcf

java -Xmx16g -jar ~/scripts/GenomeAnalysisTK_3_6.jar \
   -T HaplotypeCaller \
   -R ./long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta \
   -I merged_realigned.bam \
   -stand_call_conf 30.0 \
   -stand_emit_conf 10.0 \
   -L intergenic.bed \
   -o raw_snps_indels_Q30_intergenic.vcf
