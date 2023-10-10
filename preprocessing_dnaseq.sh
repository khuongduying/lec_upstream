# PATH
p_work='/DATA/bio06/duydao/lec_preprocessing'
raw_bam
raw_fq


##~~~~~~~~------------
# SAMtools
# sort paired read alignment .bam file (sort by name -n)
##Syntax: samtools sort -n <file.bam> -o <file_sorted.bam>

for file in $(ls *chr22.bam); do samtools sort -n $file -o ${file//\.bam/}_sorted.bam;done

# Example of string subtitution using pure bash: for file in $(ls *chr22.bam); do echo ${file//\.bam/};done
# Why sorted .bam has more line than the original?

# save fastq reads in separate R1 and R2 files

samtools fastq -@ 4 $p_work/0raw_bam/ERR5862175_chr22_sorted.bam \
-1 $p_work/1raw_fq/ERR5862175_chr22_R1.fastq.gz \
-2 $p_work/1raw_fq/ERR5862175_chr22_R2.fastq.gz \
-0 /dev/null -s /dev/null -n
"""[M::bam2fq_mainloop] discarded 43095 singletons
[M::bam2fq_mainloop] processed 3747402 reads"""

samtools fastq -@ 4 $p_work/0raw_bam/SRR9216862_chr22_sorted.bam \
-1 $p_work/1raw_fq/SRR9216862_chr22_R1.fastq.gz \
-2 $p_work/1raw_fq/SRR9216862_chr22_R2.fastq.gz \
-0 /dev/null -s /dev/null -n

"""[M::bam2fq_mainloop] discarded 32909 singletons
[M::bam2fq_mainloop] processed 14213663 reads"""

# 

samtools bam2fq -@ 4 $p_work/0raw_bam/ERR5862175_chr22_sorted.bam \
-1 $p_work/2raw_fq/ERR5862175_chr22_R1.fastq.gz \
-2 $p_work/2raw_fq/ERR5862175_chr22_R2.fastq.gz \
-0 /dev/null -s /dev/null -n

#~~~~~~~~~

# FastQC

for i in $(ls *.fastq.gz); do fastqc $i -o qc_checked; done

'''[bio06]
real	2m50,526s
user	2m56,306s
sys	0m2,295s
'''

# trim
fastp

fastp \
-i in.R1.fq.gz \
-I in.R2.fq.gz \
-o out.R1.fq.gz \
-O out.R2.fq.gz

fastp \
-i LowQuality_Reads.fastq.gz \
-o LowQuality_Reads.fastq_trimmed.gz 

fastp \
-i 1_raw/sample_1/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz \
-I 1_raw/sample_1/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz \
-o 2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_trimmed.fastq.gz \
-O 2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_trimmed.fastq.gz

'''[bio06]
time used: 21 seconds
'''

# Trimmomatic

# using -phred64
java -jar /opt/biotools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred64 \
-trimlog SLGFSK-N_231335.log \
1_raw/sample_1/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz \
1_raw/sample_1/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_unpaired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_unpaired.fastq.gz \
ILLUMINACLIP:./adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

'''[bio06]
real	3m34,599s
user	1m46,331s
sys	1m51,236s
'''


# using -phred33
java -jar /opt/biotools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 \
-trimlog SLGFSK-N_231335.log \
1_raw/sample_1/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz \
1_raw/sample_1/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_unpaired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_unpaired.fastq.gz \
ILLUMINACLIP:./adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

'''[bio06]
real	7m18,986s
user	5m27,095s
sys	1m54,638s
'''

# using -threads 4
java -jar /opt/biotools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 \
-threads 4 \
-trimlog SLGFSK-N_231335.log \
1_raw/sample_1/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz \
1_raw/sample_1/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_unpaired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_unpaired.fastq.gz \
ILLUMINACLIP:./adapters/NexteraPE-PE.fa:2:30:10 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:36

'''[bio06]
real	3m11,778s
user	5m57,531s
sys	1m58,202s
'''
##~~~~

java -jar /opt/biotools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-phred33 \
-threads 4 \
-trimlog SLGFSK-T_231335.log \
1_raw/sample_1/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz \
1_raw/sample_1/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz \
2_trimmed/sample_1/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-T_231336_r1_chr5_12_17_unpaired.fastq.gz \
2_trimmed/sample_1/SLGFSK-T_231336_r2_chr5_12_17_paired.fastq.gz \
2_trimmed/sample_1/SLGFSK-T_231336_r2_chr5_12_17_unpaired.fastq.gz \
ILLUMINACLIP:/opt/biotools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:3:true \
HEADCROP:3 \
LEADING:3 \
TRAILING:10 \
SLIDINGWINDOW:4:15 \
MINLEN:25

'''[bio06]
real	4m52,345s
user	8m58,877s
sys	3m0,733s
'''

# Single end 
java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Q30
java -jar /opt/biotools/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
-phred33 \
-threads 4 \
-trimlog LowQuality_Reads.log \
$p_raw/sample1/LowQuality_Reads.fastq.gz \
$p_trim/sample1/LowQuality_Reads_trimmed.fastq.gz \
ILLUMINACLIP:/opt/biotools/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10:2:3 \
HEADCROP:3
CROP:3 \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:30 \
MINLEN:36

# Alignment
## Indexing the reference genome
bwa index -a bwtsw hg19.chr5_12_17.fa.gz

'''
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -a bwtsw hg19.chr5_12_17.fa.gz
[main] Real time: 221.161 sec; CPU: 220.344 sec
'''

## Alignment
bwa mem -t 4 \
/DATA/bio06/duydao/Exome/ref_hg/ref_hg19/hg19.chr5_12_17.fa \
/DATA/bio06/duydao/Exome/2_trimmed/sample_1/SLGFSK-N_231335_r1_chr5_12_17_paired.fastq.gz \
/DATA/bio06/duydao/Exome/2_trimmed/sample_1/SLGFSK-N_231335_r2_chr5_12_17_paired.fastq.gz > 3_aligned/SLGFSK-N_231335_chr5_12_17_aln.sam
'''
[main] Real time: 272.486 sec; CPU: 1108.553 sec
'''

@ST-K00265:137:HT33CBBXX:3:2225:24271:46750

cd /DATA/bio06/duydao/Exome/3_aligned/sample_1

## Convert sam to bam using samtools (install with http://www.htslib.org/download/)

samtools view -Sb SLGFSK-N_231335_chr5_12_17_aln.sam > SLGFSK-N_231335_chr5_12_17_aln.bam

# POST PROCCESSING
## Validating the SAM/BAM file (No need to run)
java -jar $picard_path/picard.jar ValidateSamFile \
--INPUT 3_aligned/sample_1/SLGFSK-N_231335_chr5_12_17_aln.bam \
--MODE SUMMARY

## Sorting and marking duplicates, and indexing the BAM file
## Sort BAM with picardtools
### Install picard
### Requires Java-17

picard_path='/opt/biotools/picard/build/libs'

java -jar $picard_path/picard.jar SortSam \
--INPUT 3_aligned/sample_1/SLGFSK-N_231335_chr5_12_17_aln.bam \
--OUTPUT 4_aln_sorted/sample_1/SLGFSK-N_231335_chr5_12_17_aln_sorted.bam \
--SORT_ORDER coordinate

'''1.68 minutes'''

## Make duplications
java -jar $picard_path/picard.jar MarkDuplicates \
--INPUT 4_aln_sorted/sample_1/SLGFSK-N_231335_chr5_12_17_aln_sorted.bam \
--OUTPUT 5_aln_deduped/sample_1/SLGFSK-N_231335_chr5_12_17_aln_dedup.bam \
--METRICS_FILE SLGFSK-N_231335.metrics


## Check
samtools view -c -f 0x400 5_aln_deduped/sample_1/SLGFSK-N_231335_chr5_12_17_aln_dedup.bam 
#3767423

samtools view -c -f 0x400 4_aln_sorted/sample_1/SLGFSK-N_231335_chr5_12_17_aln_sorted.bam
#0



## Re-alignment

