### In this repository, I list all the commands and analyses conducted using the Skeletonema historical core data.

### First, for the population level analysis, I acquired bam files from here:

/nobackup/data6/skeletonema_resequencing/gene_search/20160504_GenStrains/RE_RUN_MF02.5_478-494/*_1_cleaned/bam/*rg.bam
(n=2)

/nobackup/data6/skeletonema_resequencing/gene_search/20160504_GenStrains/RE_RUN/*_1_cleaned/bam/*rg.bam.gz
(n=14)

/nobackup/data6/skeletonema_resequencing/gene_search/20160504_GenStrains/LONG/*_1_cleaned/bam/*rg.bam.gz
(n=125, but MF02.5 478 and 494 are corrupt)

And the reference fasta from here:
/nobackup/data6/skeletonema_resequencing/gene_search/20160323_SNPs_HFMF/long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta

### Then, using the interval information listed in /nobackup/data6/skeletonema_resequencing/gene_search/Pierre_SNPs/ProtCodReg_sorted_Introns-Exons_All.xslx, I created three bed files.

exons.bed
introns.bed
intergenic.bed

### I merged the bam files per layer using samtools merge.


### Then, using the Popoolation2 pipeline, I first used samtools mpileup to call SNPs, for each of the three regions:

>samtools mpileup -B *.bam -l introns.bed --reference ~/long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta > mpileup_introns

>java -ea -Xmx7g -jar ~/popoolation2_1201/mpileup2sync.jar --input mpileup_introns --output intron.sync --fastq-type sanger --min-qual 10 --threads 16 

>samtools mpileup -B *.bam -l exons.bed --reference ~/long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta > mpileup_exons

>java -ea -Xmx7g -jar ~/popoolation2_1201/mpileup2sync.jar --input mpileup_exons --output exon.sync --fastq-type sanger --min-qual 10 --threads 16 

>samtools mpileup -B *.bam -l intergenic.bed --reference ~/long_fwrev_fixedheaders_c44fix_3discarded_7_mar.fasta > mpileup_intergenic

>java -ea -Xmx7g -jar ~/popoolation2_1201/mpileup2sync.jar --input mpileup_intergenic --output intergenic.sync --fastq-type sanger --min-qual 10 --threads 16

# Then, calcuating allele frequencies and allele frequency differences at all variant sites:

>perl ~/popoolation2_1201/snp-frequency-diff.pl --input intron.sync --output-prefix intron --min-count 1 --min-coverage 1 --max-coverage 1000

>perl ~/popoolation2_1201/snp-frequency-diff.pl --input exon.sync --output-prefix exon --min-count 1 --min-coverage 1 --max-coverage 1000

>perl ~/popoolation2_1201/snp-frequency-diff.pl --input intergenic.sync --output-prefix intergenic --min-count 1 --min-coverage 1 --max-coverage 1000

# Doing a Fisher's exact test to determine statistical significance of allele frequcney changes at all SNP sites.
# NOTE: Worked only after PROPER install of perl lib "Text::NSP::Measures::2D::Fisher::twotailed" 
# Also, need to export env variable first: export PERL5LIB=~/MyPerlLib

>perl ~/popoolation2_1201/fisher-test.pl --input intron.sync --output intron.fet --min-count 2 --min-coverage 20 --max-coverage 1000 --suppress-noninformative

>perl ~/popoolation2_1201/fisher-test.pl --input exon.sync --output exon.fet --min-count 2 --min-coverage 20 --max-coverage 1000 --suppress-noninformative

>perl ~/popoolation2_1201/fisher-test.pl --input intergenic.sync --output intergenic.fet --min-count 2 --min-coverage 20 --max-coverage 1000 --suppress-noninformative


# Calculating pairwise Fst between all pops:

>perl ~/popoolation2_1201/fst-sliding.pl --input MF_intron.sync --output intron.fst --suppress-noninformative --min-count 2 --min-coverage 20 --max-coverage 1000 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 18:9:14:17:21:35:25

>perl ~/popoolation2_1201/fst-sliding.pl --input exon.sync --output exon.fst --suppress-noninformative --min-count 2 --min-coverage 20 --max-coverage 1000 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 18:9:14:17:21:35:25

>perl ~/popoolation2_1201/fst-sliding.pl --input intergenic.sync --output intergenic.fst --suppress-noninformative --min-count 2 --min-coverage 20 --max-coverage 1000 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 18:9:14:17:21:35:25
# Pop sizes : MF: 21:35:25   HF: 18:9:14:17