### In this repository, I list all the commands and analyses run in early 2018 to analyze the low-coverage Skeletonema data using ANGSD.

### The purpose was to generate Fst values of sliding windows, while taking uncertainties at all levels into account.

### First, for the population level analysis, I acquired new bam files, mapped by Matt Pinder, and a new shortened reference fasta, as the original one had some problems.

bam files were found on Albiorix, at:
/proj/data17/skeletonema_resequencing/06_remapping

I made lists of the bam files from each layer (MF02_5_bams.txt, MF09_5_bams.txt and MF19_5_bams.txt, located in this repository)

### STEP 1: CALCULATING THE SITE ALLELE FREQUENCY SPECTRUM ### 

# First, need to load in all necessary modules:

>cd /proj/data17/pierre/lowCov_angsd_Dec2017/genic
>mkdir Results_SNPfilt

>module load msms/v3.2rc-b163
>module load samtools/v1.3.1
>module load angsd/v0.918
>module load ngsTools/vX.X

Then, calculating the folded site allele frequency (saf) probability matrix for each layer with ANGSD, while filtering out obviously monomorphic sites:

>angsd -bam MF02_5_bams.txt -anc long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta -ref long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta -fai long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta.fai -nInd 21 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF02_5
>angsd -bam MF09_5_bams.txt -anc long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta -ref long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta -fai long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta.fai -nInd 35 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF09_5
>angsd -bam MF19_5_bams.txt -anc long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta -ref long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta -fai long_fwrev_fixedheaders_c44fix_3+3discarded_03_Nov_2017.fasta.fai -nInd 25 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF19_5



### STEP 2: 2D-SITE FREQUENCY SPECTRA ###

# Now, use the saf's to calculate 2-dimensional site-frequency spectra for each pairwise comparison:
# (We will use these as priors for the Fst calculation below)

>cd Results_SNPfilt
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF02_5.saf.idx MF09_5.saf.idx > MF02_5_MF09_5.2dsfs
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF02_5.saf.idx MF19_5.saf.idx > MF02_5_MF19_5.2dsfs
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF09_5.saf.idx MF19_5.saf.idx > MF09_5_MF19_5.2dsfs

# Testing to plot the first one, just to see what it looks like:

>Rscript /usr/local/packages/ngsTools/Scripts/plot2DSFS.R MF02_5_MF09_5.2dsfs MF02_5-MF09_5 21-35



### STEP 3: FST CALCULATION USING THE 2D-SFS AS PRIOR ###

## First, let's calculate the Fst index:

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF09_5.saf.idx -sfs MF02_5_MF09_5.2dsfs -fstout MF02_5_MF09_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF19_5.saf.idx -sfs MF02_5_MF19_5.2dsfs -fstout MF02_5_MF19_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF09_5.saf.idx MF19_5.saf.idx -sfs MF09_5_MF19_5.2dsfs -fstout MF09_5_MF19_5 -whichFst 1

## Now, let's do it with a sliding window approach.
## Should try to keep exons and introns separate. For this, need to generate a bash scripts with one line for each exon/intron position (given by exons.bed and introns.bed in the parent folder) for each of the three pairwise layer comparisons
 
## Each exon/intron region should thus have three lines, and contain this with CHROM:START-STOP modified.

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF02_5_MF09_5.fst.idx -r CHROM:START-STOP -win 200 -step 50 >> MF02_5_MF09_5_INTRON/EXON.fst.wins
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF02_5_MF19_5.fst.idx -r CHROM:START-STOP -win 200 -step 50 >> MF02_5_MF19_5_INTRON/EXON.fst.wins
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF09_5_MF19_5.fst.idx -r CHROM:START-STOP -win 200 -step 50 >> MF09_5_MF19_5_INTRON/EXON.fst.wins

# I put this into the scripts called exon_fst.sh and intron_fst.sh
# Executing the scripts:

>./exon_fst.sh
>./intron_fst.sh

------------------------------------------------

### INTERGENIC REGIONS ###

## Let's do the same for the intergenic regions as for the genic regions, now using the nc_regions reference fasta and bams mapped to that region:

>cd /proj/data17/pierre/lowCov_angsd_Dec2017/intergenic
>mkdir Results_SNPfilt

>angsd -bam MF02_5_bams.txt -anc nc_regions_to_use.fasta -ref nc_regions_to_use.fasta -fai nc_regions_to_use.fasta.fai -nInd 21 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF02_5
>angsd -bam MF09_5_bams.txt -anc nc_regions_to_use.fasta -ref nc_regions_to_use.fasta -fai nc_regions_to_use.fasta.fai -nInd 35 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF09_5
>angsd -bam MF19_5_bams.txt -anc nc_regions_to_use.fasta -ref nc_regions_to_use.fasta -fai nc_regions_to_use.fasta.fai -nInd 25 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF19_5

cd Results_SNPfilt

# Now, use the saf's to calculate 2-dimensional site-frequency spectra for each pairwise comparison:

>/usr/local/packages/ngsTools/angsd/misc/realSFS MF02_5.saf.idx MF09_5.saf.idx > MF02_5_MF09_5.2dsfs
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF02_5.saf.idx MF19_5.saf.idx > MF02_5_MF19_5.2dsfs
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF09_5.saf.idx MF19_5.saf.idx > MF09_5_MF19_5.2dsfs

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF09_5.saf.idx -sfs MF02_5_MF09_5.2dsfs -fstout MF02_5_MF09_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF19_5.saf.idx -sfs MF02_5_MF19_5.2dsfs -fstout MF02_5_MF19_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF09_5.saf.idx MF19_5.saf.idx -sfs MF09_5_MF19_5.2dsfs -fstout MF09_5_MF19_5 -whichFst 1

# Then, using the bash script with the intergenic regions to calculate Fst over the known intergenic regions with a sliding window. 

>./intergenic_fst.sh

---------------------------------------------------

### contig 2985 revcomp missing ###

# Found out that contig_2985_revcomp should have been included in the reference but it was not, so needed to redo the Fst analysis over that contig only:
# For this, Matt Pinder remapped the raw data to the contig, bam files were found here: /proj/data17/skeletonema_resequencing/06_remapping with filenames MF##.#_2985rc_merged.bam

>cd /proj/data17/pierre/lowCov_angsd_Dec2017/genic
>mkdir Results_SNPfilt
>angsd -bam MF02_5_bams.txt -anc contig_2985_revcomp.fasta -ref contig_2985_revcomp.fasta -fai contig_2985_revcomp.fasta.fai -nInd 21 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF02_5
>angsd -bam MF09_5_bams.txt -anc contig_2985_revcomp.fasta -ref contig_2985_revcomp.fasta -fai contig_2985_revcomp.fasta.fai -nInd 35 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF09_5
>angsd -bam MF19_5_bams.txt -anc contig_2985_revcomp.fasta -ref contig_2985_revcomp.fasta -fai contig_2985_revcomp.fasta.fai -nInd 25 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF19_5

>cd Results_SNPfilt

## Fst (using the 2d-sfs data from the larger reference used above as prior):

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF09_5.saf.idx -sfs MF02_5_MF09_5.2dsfs -fstout MF02_5_MF09_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF19_5.saf.idx -sfs MF02_5_MF19_5.2dsfs -fstout MF02_5_MF19_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF09_5.saf.idx MF19_5.saf.idx -sfs MF09_5_MF19_5.2dsfs -fstout MF09_5_MF19_5 -whichFst 1

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF02_5_MF09_5.fst.idx -r contig_2985_revcomp:5085-5966 -win 200 -step 50 >> MF02_5_MF09_5_contig_2985_revcomp_EXON.fst.wins
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF02_5_MF19_5.fst.idx -r contig_2985_revcomp:5085-5966 -win 200 -step 50 >> MF02_5_MF19_5_contig_2985_revcomp_EXON.fst.wins
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF09_5_MF19_5.fst.idx -r contig_2985_revcomp:5085-5966 -win 200 -step 50 >> MF09_5_MF19_5_contig_2985_revcomp_EXON.fst.wins


------------------------------------------------------

### MICROSATELLITE REGIONS ###

## Let's do the same for the microsatellite regions as for the others:

# Also, here, I got a reference sequence and bam files provided.

>cd /proj/data17/pierre/lowCov_angsd_Dec2017/microsats
>mkdir Results_SNPfilt

>angsd -bam MF02_5_bams.txt -anc Long_BetweenPrimersResults.fasta -ref Long_BetweenPrimersResults.fasta -fai Long_BetweenPrimersResults.fasta.fai -nInd 21 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF02_5
>angsd -bam MF09_5_bams.txt -anc Long_BetweenPrimersResults.fasta -ref Long_BetweenPrimersResults.fasta -fai Long_BetweenPrimersResults.fasta.fai -nInd 35 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF09_5
>angsd -bam MF19_5_bams.txt -anc Long_BetweenPrimersResults.fasta -ref Long_BetweenPrimersResults.fasta -fai Long_BetweenPrimersResults.fasta.fai -nInd 25 -GL 2 -doSaf 1 -fold 1 -doMajorMinor 4 -doMaf 1 -snp_pval 1e-2 -out Results_SNPfilt/MF19_5

>cd Results_SNPfilt

# Now, use the saf's to calculate 2-dimensional site-frequency spectra for each pairwise comparison:

>/usr/local/packages/ngsTools/angsd/misc/realSFS MF02_5.saf.idx MF09_5.saf.idx > MF02_5_MF09_5.2dsfs
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF02_5.saf.idx MF19_5.saf.idx > MF02_5_MF19_5.2dsfs
>/usr/local/packages/ngsTools/angsd/misc/realSFS MF09_5.saf.idx MF19_5.saf.idx > MF09_5_MF19_5.2dsfs

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF09_5.saf.idx -sfs MF02_5_MF09_5.2dsfs -fstout MF02_5_MF09_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF02_5.saf.idx MF19_5.saf.idx -sfs MF02_5_MF19_5.2dsfs -fstout MF02_5_MF19_5 -whichFst 1
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst index MF09_5.saf.idx MF19_5.saf.idx -sfs MF09_5_MF19_5.2dsfs -fstout MF09_5_MF19_5 -whichFst 1

# Then, calculate Fst over the microsat regions with a sliding window. 

>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF02_5_MF09_5.fst.idx -win 200 -step 50 >> MF02_5_MF09_5_microsat.fst.wins
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF02_5_MF19_5.fst.idx -win 200 -step 50 >> MF02_5_MF19_5_microsat.fst.wins
>/usr/local/packages/ngsTools/angsd/misc/realSFS fst stats2 MF09_5_MF19_5.fst.idx -win 200 -step 50 >> MF09_5_MF19_5_microsat.fst.wins





