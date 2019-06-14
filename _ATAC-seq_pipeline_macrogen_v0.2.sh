#!/bin/bash

###################################################################################################################
###################################################################################################################
######################################### ATAC-seq Analysis Pipeline ##############################################
###################################################################################################################
############################################################### written by Min-hwan R. Sohn  ######################
####################################################################### on May 24th 2019     ######################
###################################################################################################################
############################################################## modified by Min-hwan R. Sohn  ######################
####################################################################### on June 7th 2019     ######################
###################################################################################################################
###################################################################################################################

### Edit history

#(2019-05-27) Put in extra '--bdg' option for peak calling step with MACS2 ==> generating optional bedGraph file
#(2019-05-27) Remove '--nomodel' option for peak calling step with MACS2 ==> not necessary when --f BEDPE is turned on 
#(2019-05-27) Put in extra Quality Control measure: Nonredundant Fraction (NRF) 
#             NRF = #unique start positions of uniquely mappable reads / #uniquely mappable reads (>0.8 recommended)
#(2019-06-05) TSS Enrichment, FRiP score calculation Added (for hg19, mm10) for version 0.2
#             Added features require specific bedfiles in working directory. 
#             _TSS-enrichment-score_calculation.py and several bedfiles needed! see below for details
#(2019-06-07) Primary contig selection for mm10 is added ==> works by manual code modification
#(2019-06-07) code error in Blacklisted region filtering ==> amended
#(2019-06-13) MT proportion counts code error ==> amended
#(2019-06-13) NRF code error ==> amended

### Analysis pipeline Overview ###

## 1.Adapter and low quality bases Trimming and Reporting QC ##
## 2.Sequence Alignment ##
## 3-1.Mitochondrial Reads ##
## 3-2.Read Duplicates Removal ##
## 4.Peak Calling ##
## 5.ATAC-seq quality assessment via checking on Fragment Length Distribution, FRiP score, TSS enrichment...etc. ##
##
## (optional) Concordance or Correlation between replicates ##
## (optional) Differential Peaks Identification ##
## (optional) Motif Enrichment ##

## Required softwares, packages and annotation files

## TrimGalore! (and cutadapt)
## Bowtie2
## Samtools
## Picard
## MACS2 (MUST be the latest one!)
## !!! NOTE: If your samples of interest are Mouse ==> --gsize mm
##                                           C.Elegans ==> --gsize ce
##                                           D.melanogaster ==> --gsize dm

## Reference Genome indexed using bowtie2-build (ends with *.bt2) ==> refdir variable
## Reference Genome chromSize file
## !!! NOTE: File format Example
## chr1    249250621
## chr2    243199373
##         .
##         .
##         .
## !!! NOTE: Please check if Mitochodrial genome header is chrMT or chrM
##           If it's chrM, then modify script below to substitute chrM for chrMT
## Gene annotation file (GTF file)
## Transcription Start Site bedfile from GTF file for TSS enrichment Calculation
## NOTE: 
## ENCODE Blacklisted genomic regions 
## !!! NOTE: Try download here http://mitra.stanford.edu/kundaje/akundaje/release/blacklists
##           For hg19: http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg19-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz
## (optional) HOMER version 4.9
## (optional) Several R packages e.g. DESeq2
## _TSS-enrichment-score_calculation.py code
## Refseq TSS+1kb-1kb file and TSS+1kb-1kb_flanking100bp file for TSS enrichment score calculation

###################################################################################################################

## Softwares, packages and annotation files directory

bowtie2=''
java=''
picard=''
samtools=''
bedtools=''
python=''
trim=''
macs2=''
refdir='' # Directory where you put your bowtie2 indexed reference file and basename ex)/clinix1/Analysis/mongol/phenomata/Reference/hg19/hg19_primary ==> where hg19_primary.1.bt2 is
refsize='' # ChromSize file for bedtools slop. See upper description for details on this.
blacklist='' # ENCODE blacklisted region bed file
tss_updown1kb='' # hg19_refseq_TSS+1kb-1kb.sorted.bed, mm10_refseq_TSS+1kb-1kb.sorted.bed
tss_updown1kb_flanking100bp='' # hg19_refseq_TSS+1kb-1kb_flanking100bp.sorted.bed, mm10_refseq_TSS+1kb-1kb_flanking100bp.sorted.bed

###################################################################################################################

### Usage ###

## nohup sh _ATAC-seq_pipeline.sh sample(a) sample(a+k) samplelist.txt n > NAME.out &
## where, 'a' represents the integer number from which the analysis begins and 'a+k' refers to where you want to end the analyiss.
## 'n' is the user-selected number of CPU core usage (for Bowtie2 and samtools)

## check 


## Input file (samplelist.txt) format ==> SampleName read1.fastq.gz read2.fastq.gz
## Input file example
## AK1-rep1 AK1-rep1_R1.fastq.gz AK1-rep1_R2.fastq.gz
## AK1-rep2 AK1-rep2_R1.fastq.gz AK1-rep2_R2.fastq.gz
## AK1-rep3 AK1-rep3_R1.fastq.gz AK1-rep3_R2.fastq.gz

###################################################################################################################

for((i=$1;i<=$2;i++))
do
	sed -n ${i}p $3 > tmp
	sample=$(awk '{print $1}' tmp)
	read1=$(awk '{print $2}' tmp)
	read2=$(awk '{print $3}' tmp)

echo $sample $read1 $read2

echo 'MESSAGE: ATAC-seq Analysis Pipeline Activated'

echo 'MESSAGE: Now Generating sub-directories'
mkdir -p 01.Alignment
mkdir -p 02.Peak_Calling
mkdir -p 03.Analysis

###################################################################################################################

## 1.Adapter and low quality bases Trimming and Reporting QC ##
echo 'MESSAGE: Now Performing Step1.Adapter and low quality bases Trimming and Reporting QC using TrimGalore!'

# !!! Please make sure to specify a path to the Cutadapt by --path_to_cutadapt </path/to/cutadapt> or Else it is assumed that Cutadapt is in the PATH !!!
# Assuming TruSeq adapters...
${trim} \
--nextera \
--path_to_cutadapt /root/.local/bin/cutadapt \
--quality 15 \
--phred33 \
--fastqc \
--trim-n \
--gzip \
--length 36 \
--paired \
${read1} \
${read2} 

echo 'MESSAGE: Step1 has been completed'

## 2.Sequence Alignment ##
echo 'MESSAGE: Now Performing Step2.Sequence Alignment using Bowtie2'

mkdir -p ./01.Alignment/Unaligned
mkdir -p ./01.Alignment/Metrics
mkdir -p ./01.Alignment/Aligned

read1_1=`echo ${read1} | sed 's/.fastq.gz//'`
read2_1=`echo ${read2} | sed 's/.fastq.gz//'`

${bowtie2} \
-q \
--phred33 \
--maxins 2000 \
--un-gz ./01.Alignment/Unaligned \
--time \
--met-file ./01.Alignment/Metrics/${sample}_metrics \
--no-unal \
--very-sensitive \
--no-discordant \
--no-mixed \
--mm \
--threads ${4} \
-x ${refdir} \
-1 ${read1_1}_val_1.fq.gz \
-2 ${read2_1}_val_2.fq.gz \
-S ./01.Alignment/Aligned/${sample}.sam 

${samtools} \
view \
-@ ${4} \
-bS ./01.Alignment/Aligned/${sample}.sam \
> ./01.Alignment/Aligned/${sample}.bam

/bin/rm -rf ./01.Alignment/Aligned/${sample}.sam

${samtools} \
sort \
-@ ${4} \
./01.Alignment/Aligned/${sample}.bam \
-o ./01.Alignment/Aligned/${sample}.sorted.bam  

/bin/rm -rf ./01.Alignment/Aligned/${sample}.bam

${samtools} \
index \
./01.Alignment/Aligned/${sample}.sorted.bam

${samtools} \
flagstat \
./01.Alignment/Aligned/${sample}.sorted.bam \
> ./01.Alignment/Metrics/${sample}.sorted.flagstat.txt

echo 'MESSAGE: Step2 has been completed'

## 3-1.Mitochondrial Reads ##
echo 'MESSAGE: Now performing Step3-1.Mitochondrial Reads'

# if it's mm10 your working on please specify the chromosomes by: chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
${samtools} \
view \
-@ ${4} \
-b ./01.Alignment/Aligned/${sample}.sorted.bam \
chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
> ./01.Alignment/Aligned/${sample}.sorted.noMT.bam

# QC with the proportion of mitochondrial reads 
mtcount=`${samtools} view -c -@ ${4} ./01.Alignment/Aligned/${sample}.sorted.bam chrMT | awk '{print $1}'`
allcount=`${samtools} view -c -@ ${4} ./01.Alignment/Aligned/${sample}.sorted.bam | awk '{print $1}'`
mtproportion=`awk "BEGIN { print ${mtcount} * 100 / ${allcount} }"`
echo 'Mitochondrial reads Proportion:' ${mtproportion} '%'

echo 'MESSAGE: Step3-1 has been completed'

## 3-2.Read Duplicates Removal ##
echo 'MESSAGE: Now performing Step3-2.Read Duplicates Removal'

${java} -Xmx5g -jar ${picard} \
MarkDuplicates \
INPUT=./01.Alignment/Aligned/${sample}.sorted.noMT.bam \
OUTPUT=./01.Alignment/Aligned/${sample}.sorted.noMT.dp.bam \
METRICS_FILE=./01.Alignment/Metrics/${sample}.sorted.noMT.dp.metrix \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=true

${samtools} \
index \
./01.Alignment/Aligned/${sample}.sorted.noMT.dp.bam 

# QC with Nonredundant fraction (NRF)
uspumr=`${samtools} view -c -@ ${4} ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.bam | awk '{print $1}'`
umr=`${samtools} view -c -@ ${4} ./01.Alignment/Aligned/${sample}.sorted.noMT.bam | awk '{print $1}'`
nrf=`awk "BEGIN { print ${uspumr} / ${umr} }"`
echo 'Nonredundant fraction (NRF):' ${nrf}

${samtools} \
flagstat \
./01.Alignment/Aligned/${sample}.sorted.noMT.dp.bam \
> ./01.Alignment/Metrics/${sample}.sorted.noMT.dp.flagstat.txt

# Filtering for high Quality, properly paired reads###
${samtools} \
view \
-@ ${4} \
-bq 30 \
-f 0x2 \
./01.Alignment/Aligned/${sample}.sorted.noMT.dp.bam \
> ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bam

${samtools} \
index \
./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bam

${samtools} \
flagstat \
./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bam \
> ./01.Alignment/Metrics/${sample}.sorted.noMT.dp.proper.flagstat.txt

echo 'MESSAGE: Step3-2 has been completed'

## 4.Peak Calling ##
echo 'MESSAGE: Now performing Step4.Peak Calling using MACS2'

# Peak calling using MACS2
${macs2} \
callpeak \
--treatment ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bam \
--format BAMPE \
--gsize hs \
--bdg \
--nolambda \
--keep-dup all \
--broad \
--outdir ./02.Peak_Calling \
-n ${sample}

# Filter blacklisted region
${bedtools} subtract \
-a ./02.Peak_Calling/${sample}_peaks.broadPeak \
-b $blacklist \
> ./02.Peak_Calling/${sample}_peaks.blfiltered.broadPeak

# Generating BigWig file for visualization


echo 'MESSAGE: Step4 has been completed'

## 5.ATAC-seq quality assessment via checking on Fragment Length Distribution, FRiP score, TSS enrichment...etc. ##
echo 'MESSAGE: Now performing Step5.ATAC-seq quality assessment via checking on Fragment Length Distribution, FRiP score, TSS enrichment...etc.'

# TSS enrichment
${bedtools} bamtobed \
-i ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bam \
> ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bed

sort -k1,1 -k2,2n ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bed > ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.sorted.bed

$python \
$bedtools \
_TSS-enrichment-score_calculation.py \
$tss_updown1kb \
$tss_updown1kb_flanking100bp \
./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.sorted.bed 

tes=`awk '{print $1}' tmp3.bed`
echo 'Transcription Start Site Enrichment Score (TSSES): ' ${tes}

/bin/rm -rf tmp3.bed

# FRiP score
frip=`${bedtools} coverage \
-counts \
-sorted \
-a ./02.Peak_Calling/${sample}_peaks.blfiltered.broadPeak \
-b ./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.sorted.bed \
| awk 'BEGIN { total=0 } { total = total + $NF } END { printf total }'`
echo 'Fraction of reads in Peaks (FRiP) score: ' ${frip}

# Fragment Length Distribution
${java} -Xmx5g -jar ${picard} \
CollectInsertSizeMetrics \
I=./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.bam \
O=./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.txt \
H=./01.Alignment/Aligned/${sample}.sorted.noMT.dp.proper.pdf


echo 'MESSAGE: Step5 has been completed'
echo 'MESSAGE: Sequence Procesing Complete !!'
done
###################################################################################################################
