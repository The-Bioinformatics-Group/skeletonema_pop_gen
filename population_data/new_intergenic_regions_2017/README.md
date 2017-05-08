### As the "intergenic regions" found in the reference file previously provided (the 80 "long contigs") were indicated to perhaps not be so intergenic, we decided to find new regions.

I was provided with a new reference file, now on Albiorix in: /nobackup/data8/pierre_skeletonema_project/intergenic_Feb2017/nc_regions_to_use.fasta

### The first task was to map the raw data against this new reference file in the same way as had been done for the other regions, my mapping script is here:

https://github.com/The-Bioinformatics-Group/skeletonema_pop_gen/blob/master/population_data/bwa_map_skeletonema.sge

Briefly, the scripts is mapping all MF strains against the reference using bwa mem (default settings), sorting and converting to bam, then removing duplicated with Picard MarkDuplicates.

All bam files can be found on Albiorix, here: /nobackup/data8/pierre_skeletonema_project/intergenic_Feb2017/MF_data

### I merged the bam files per layer using samtools merge (three files: MF02_5.bam, MF09_5.bam and MF19_5.bam)

### Then, using the Popoolation2 pipeline, I did the same thing as for the other regions described in the previous page. All commands are found here:

https://github.com/The-Bioinformatics-Group/skeletonema_pop_gen/blob/master/population_data/skeletonema_mpileup_script_2017.sge

