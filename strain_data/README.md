### On this page, I describe the pipeline I have used in order to analyse the Deeply sequenced strains on individual basis, for use in the genotype-phenotype mapping performed by Staffan Nilsson and Marina Axelson-Fisk (Chalmers)

### First, I needed to map the new data against the reference fasta, which was the "80 long contigs" used also for the population data:

#Reference fasta:
/nobackup/data6/skeletonema_resequencing/gene_search/20160323_SNPs_HFMF/long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta

I did this using bowtie2, with the command (e.g. for Sample 1_160225_AC86M1ANXX_P4005_1003):

> bowtie2-2.2.7/bowtie2-align-s --wrapper basic-0 -x long -p 8 --passthrough -1 1_160225_AC86M1ANXX_P4005_1003.paired.1.fastq -2 1_160225_AC86M1ANXX_P4005_1003.paired.2.fastq -U 1_160225_AC86M1ANXX_P4005_1003.unpaired.1.fastq,1_160225_AC86M1ANXX_P4005_1003.unpaired.2.fastq

Then sorting and converting to bam using convert_to_bam.sh 

All sorted bam files are located in (on Albiorix):

/nobackup/data8/pierre_skeletonema_project/skeletonema_DeepSeq_strains

### Then, merging the bam files for joint SNP calling, while retaining read group information:

samtools merge -h rg.txt merged.bam *.sorted.bam

rg.txt is here: (https://github.com/The-Bioinformatics-Group/skeletonema_pop_gen/blob/master/strain_data/rg.txt)

### InDel realigning the merged bam file using GATK InDel realigner, using the script realigner.sh

> realigner.sh

(https://github.com/The-Bioinformatics-Group/skeletonema_pop_gen/blob/master/strain_data/realigner.sh)

## Then, calling SNPs and InDels using GATK Unified Genotyper, with a quality cutoff of Q 30. 
This is done separately for the three different regions (exon, intron and intergenic, specified in the bed files).

> SNP_call.sh

(https://github.com/The-Bioinformatics-Group/skeletonema_pop_gen/blob/master/strain_data/SNP_call.sh)

