#!/bin/bash

#PBS -N GENOTYPE_merge
#PBS -l ncpus=8
#PBS -l mem=72gb
#PBS -l walltime=48:00:00
#PBS -M christine.d.roper@student.uts.edu.au
#PBS -m abe

export MODULEPATH=/shared/c3/apps/modules:$MODULEPATH

echo "#####################START#####################"

module load bio/fastqc-current
module load bio/trimmomatic-current
module load bio/samtools-current
module load bio/bwa-current
module load bio/picard-current
module load bio/gatk-current
module load bio/star-current
module load bio/stringtie-current
module load bio/HTSeq-current

cd /shared/c3/projects/Future_Reefs/coral_deoxygenation.christine/

ref="/shared/c3/projects/Future_Reefs/coral_deoxygenation.christine/Pocillopora_acuta_ref/Pocillopora_acuta_genome_v1.fasta"
GFF="/shared/c3/projects/Future_Reefs/coral_deoxygenation.christine/Pocillopora_acuta_ref/Structural_annotation_experimental.gff"
GTF="/shared/c3/projects/Future_Reefs/coral_deoxygenation.christine/Pocillopora_acuta_ref/Structural_annotation_experimental.gtf"

#step 1 merge fastqz
#1.1 R1
#cat \
#./raw_data/GENOTYPE_L001_R1.fastq.gz \
#./raw_data/GENOTYPE_L002_R1.fastq.gz \
#./raw_data/GENOTYPE_L003_R1.fastq.gz \
#./raw_data/GENOTYPE_L004_R1.fastq.gz \
#./raw_data/GENOTYPE_L005_R1.fastq.gz \
#> ./merged_fastqz/GENOTYPE_merged_R1.fastq.gz

#1.2 R2
#cat \
#./raw_data/GENOTYPE_L001_R2.fastq.gz \
#./raw_data/GENOTYPE_L002_R2.fastq.gz \
#./raw_data/GENOTYPE_L003_R2.fastq.gz \
#./raw_data/GENOTYPE_L004_R2.fastq.gz \
#./raw_data/GENOTYPE_L005_R2.fastq.gz \
#> ./merged_fastqz/GENOTYPE_merged_R2.fastq.gz

#step 2 FASTQC
#fastqc \
#-o ./fastqc/ \
#-t 4 \
#./merged_fastqz/GENOTYPE_merged_R1.fastq.gz

#fastqc \
#-o ./fastqc/ \
#-t 4 \
#./merged_fastqz/GENOTYPE_merged_R2.fastq.gz

#step 3 trim adapter
#java -jar /shared/c3/apps/bio/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar \
#PE \
#-threads 8 \
#./merged_fastqz/GENOTYPE_merged_R1.fastq.gz  \
#./merged_fastqz/GENOTYPE_merged_R2.fastq.gz \
#./trim/GENOTYPE_1.paired.fastq.gz \
#./trim/GENOTYPE_1.unpaired.fastq.gz \
#./trim/GENOTYPE_2.paired.fastq.gz \
#./trim/GENOTYPE_2.unpaired.fastq.gz \
#ILLUMINACLIP:/shared/c3/apps/bio/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 \
#LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
#-phred64
#-phred33

#step 3 align to reference genome
STAR \
--runThreadN 8 \
--genomeDir /shared/c3/projects/Future_Reefs/coral_deoxygenation.christine/Pocillopora_acuta_ref/ \
--readFilesIn \
./trim/GENOTYPE_1.paired.fastq.gz \
./trim/GENOTYPE_2.paired.fastq.gz \
--outFileNamePrefix ./star/bam_star/GENOTYPE_ \
--readFilesCommand zcat \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate

#step 4 gene count
#step 4.1 stringtie
#stringtie \
#./star/bam_star/GENOTYPE_Aligned.sortedByCoord.out.bam \
#-p 8 \
#-G $GTF \
#-e \
#-o ./stringtie/GENOTYPE.gtf \
#-A ./stringtie/GENOTYPE.gene.abundances.tsv -B

#step 4.2 also HT-seq
htseq-count \
-f bam \
-t gene \
-i ID \
./star/bam_star/GENOTYPE_Aligned.sortedByCoord.out.bam \
$GFF \
--stranded=no \
>./htseq/GENOTYPE.htseq.out

#step 5 on local

echo "#####################END2#####################"
