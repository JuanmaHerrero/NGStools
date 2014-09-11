#!/bin/bash
#SBATCH --time=10000
#SBATCH --mem=16000
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --output=outputGG-WUR-002_%j.txt
#SBATCH --error=error_outputGG-WUR-002_%j.txt
#SBATCH --job-name=GG-WUR-002
#SBATCH --partition=ABGC_Std
#SBATCH --constraint=normalmem
mkdir tmpGG-WUR-002
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'starting time: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
# archive number 1: ABGSAGG0002
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/backup/WUR/ABGC/shared/ABGC_datastore/SequenceData/Chicken/ABGSA/ABGSAGG0002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.gz | pigz >tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/backup/WUR/ABGC/shared/ABGC_datastore/SequenceData/Chicken/ABGSA/ABGSAGG0002/C40J5ANXX_102445-02_AGTTCC_L006_R2.fastq.gz | pigz >tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.gz -r tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R2.fastq.gz -o tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.tr -p tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R2.fastq.tr -s tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.tr
pigz tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'starting bwa-mem mapping of GG-WUR-002 archive 1: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSAGG0002_1\tSM:GG-WUR-002\tPL:ILLUMINA' /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R1.fastq.tr.gz tmpGG-WUR-002/C40J5ANXX_102445-02_AGTTCC_L006_R2.fastq.tr.gz >tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Shb tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sam > tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.bam
rm tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sam
echo 'start sorting'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.bam tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sorted
rm tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'finished, produced BAM file tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sorted.bam archive 1: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sorted.bam`; echo "size of file tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sorted.bam is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
# archive number 2: ABGSAGG0002
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/backup/WUR/ABGC/shared/ABGC_datastore/SequenceData/Chicken/ABGSA/ABGSAGG0002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.gz | pigz >tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.gz
python2 /cm/shared/apps/WUR/ABGC/abgsascripts/fix_fq_names.py /lustre/backup/WUR/ABGC/shared/ABGC_datastore/SequenceData/Chicken/ABGSA/ABGSAGG0002/C46M7ANXX_102445-02_AGTTCC_L008_R2.fastq.gz | pigz >tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R2.fastq.gz
# quality trimming of reads by sickle
sickle pe -f tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.gz -r tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R2.fastq.gz -o tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.tr -p tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R2.fastq.tr -s tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.singles.tr -l 50 -t sanger
pigz tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.tr
pigz tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R2.fastq.tr
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'starting bwa-mem mapping of GG-WUR-002 archive 2: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
# maping using the bwa-mem algorithm, including sorting of bam
echo 'start mapping using BWA-mem algorithm'
/cm/shared/apps/WUR/ABGC/bwa/bwa-0.7.5a/bwa mem -t 4 -M -R '@RG\tID:ABGSAGG0002_2\tSM:GG-WUR-002\tPL:ILLUMINA' /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R1.fastq.tr.gz tmpGG-WUR-002/C46M7ANXX_102445-02_AGTTCC_L008_R2.fastq.tr.gz >tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -Shb tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sam > tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.bam
rm tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sam
echo 'start sorting'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools sort tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.bam tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sorted
rm tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'finished, produced BAM file tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sorted.bam archive 2: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sorted.bam`; echo "size of file tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sorted.bam is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
#number of bams: 2
#multiple bam files --> do merge
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sorted.bam | sed 's/SM:unknown/SM:GG-WUR-002/' | sed 's/PL:sanger/PL:ILLUMINA/' >tmpGG-WUR-002/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools view -H tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sorted.bam | sed 's/SM:unknown/SM:GG-WUR-002/' | sed 's/PL:sanger/PL:ILLUMINA/' | grep @RG >>tmpGG-WUR-002/newheader.txt
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools merge tmpGG-WUR-002/tmpmergedGG-WUR-002.bam tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-1-pe.sorted.bam tmpGG-WUR-002/aln-ABGSAGG0002-GG-WUR-002-2-pe.sorted.bam 
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools reheader tmpGG-WUR-002/newheader.txt tmpGG-WUR-002/tmpmergedGG-WUR-002.bam >tmpGG-WUR-002/GG-WUR-002_rh.bam
rm tmpGG-WUR-002/tmpmergedGG-WUR-002.bam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools index tmpGG-WUR-002/GG-WUR-002_rh.bam
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced BAM file tmpGG-WUR-002/GG-WUR-002_rh.bam: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.bam`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.bam is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5BAM=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.bam is "$MD5BAM >>tmpGG-WUR-002/GG-WUR-002.log
# dedup using samtools
echo 'dedupping using samtools'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools rmdup tmpGG-WUR-002/GG-WUR-002_rh.bam tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools index tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam
cp tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam.bai tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced BAM file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5BAM=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam is "$MD5BAM >>tmpGG-WUR-002/GG-WUR-002.log
# re-alignment using GATK-RealignmentTargetCreator+IndelRealigner
java -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -T RealignerTargetCreator -R /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa -I tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam -o tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.intervals
java -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -T IndelRealigner -R /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa -I tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.bam -targetIntervals tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.intervals -o tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam
cp tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bai tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam.bai
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced BAM file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5BAM=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam is "$MD5BAM >>tmpGG-WUR-002/GG-WUR-002.log
# Calculate coverage statistics
java -Xmx8G -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -T DepthOfCoverage -R /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa -I tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam --omitDepthOutputAtEachBase --logging_level ERROR --summaryCoverageThreshold 10 --summaryCoverageThreshold 20 --summaryCoverageThreshold 30 --summaryCoverageThreshold 40 --summaryCoverageThreshold 50 --summaryCoverageThreshold 80 --summaryCoverageThreshold 90 --summaryCoverageThreshold 100 --summaryCoverageThreshold 150 --minBaseQuality 15 --minMappingQuality 30 --start 1 --stop 1000 --nBins 999 -dt NONE -o tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam.coverage
# old-school variant calling using the pileup algorithm
echo 'old-school variant calling using the pileup algorithm'
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools view -u tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools pileup -vcf /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa - >tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA_vars-raw.txt
VAR=`cat tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA_vars-raw.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let MAX=2*VAR
echo "max depth is $MAX"
MIN=$(( $VAR / 3 )) ; if [ $MIN -lt 5 ]; then MIN=4; fi
echo "min_depth is $MIN"
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/misc/samtools.pl varFilter -D$MAX -d$MIN tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA_vars-raw.txt | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' >tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced variant file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5VAR=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt is "$MD5VAR >>tmpGG-WUR-002/GG-WUR-002.log
# variant calling using the mpileup function of samtools
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/samtools mpileup -C50 -ugf /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/bcftools view -bvcg -| /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/bcftools view - | perl /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.19/bcftools/vcfutils.pl varFilter -D $MAX -d $MIN >tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced variant file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5VAR=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf is "$MD5VAR >>tmpGG-WUR-002/GG-WUR-002.log
# Variant calling using GATK UnifiedGenotyper - parameters need tweaking
java -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK2.6/GenomeAnalysisTK.jar -nt 4 -R /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa -T UnifiedGenotyper -I tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam --genotype_likelihoods_model BOTH -o tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200
bgzip tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf
tabix -p vcf tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced variant file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf.gz: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf.gz`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf.gz is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5VAR=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf.gz | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.UG.raw.vcf.gz is "$MD5VAR >>tmpGG-WUR-002/GG-WUR-002.log
# Create gVCF file using modified GATK UnifiedGenotyper - parameters need tweaking
/cm/shared/apps/WUR/ABGC/gvcftools/gvcftools-0.16/bin/getBamAvgChromDepth.pl tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam >tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.avgdepth.txt
java -Xmx8g -jar /cm/shared/apps/WUR/ABGC/GATK/GATK_gVCFmod/GenomeAnalysisTK.jar -nt 4 -R /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa -T UnifiedGenotyper -I tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam -glm BOTH -l OFF -stand_call_conf 20.0 -stand_emit_conf 10.0 -dcov 200 -out_mode EMIT_ALL_SITES | /cm/shared/apps/WUR/ABGC/gvcftools/gvcftools-0.16/bin/gatk_to_gvcf --chrom-depth-file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.avgdepth.txt | bgzip -c >tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz
tabix -p vcf tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'Produced variant file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
FSIZE=`stat --printf="%s" tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz`; echo "size of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz is "$FSIZE >>tmpGG-WUR-002/GG-WUR-002.log
MD5VAR=`md5sum tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz | sed 's/ \+/	/' | cut -f1`; echo "md5sum of file tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.gvcf.gz is "$MD5VAR >>tmpGG-WUR-002/GG-WUR-002.log
DATE=`date`; echo "++++++++++++++++++++++++++++" >>tmpGG-WUR-002/GG-WUR-002.log; echo 'finished variant calling: '$DATE >>tmpGG-WUR-002/GG-WUR-002.log
# predicting function using VEP
perl /cm/shared/apps/WUR/ABGC/variant_effect_predictor/VEP75/variant_effect_predictor.pl -i tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf --dir /cm/shared/apps/WUR/ABGC/variant_effect_predictor/VEP75//cache --species gallus_gallus -o tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vep.txt --fork 4 --canonical --symbol --sift b --no_intergenic --offline --force_overwrite
# course nucleotide diversity stat generator
VAR=`cat tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.vars-flt_final.txt | cut -f8 | head -1000000 | sort | uniq -c | sed 's/^ \+//' | sed 's/ \+/\t/' | sort -k1 -nr | head -1 | cut -f2`
let VAR=2*VAR
echo "max depth is $VAR"
/cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools view -u tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.bam | /cm/shared/apps/WUR/ABGC/samtools/samtools-0.1.12a/samtools pileup -f /lustre/nobackup/WUR/ABGC/shared/public_data_store/genomes/chicken/Ensembl74/Gallus_gallus.Galgal4.74.dna.toplevel.fa -c - | awk '$8>4' | awk -v VAR=$VAR '$8<VAR' | perl /cm/shared/apps/WUR/ABGC/abgsascripts/extract_stats-pileup-bins_allchroms.pl -f tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA
# Annovar: Convert into annovar format and annotate genetic variants
/cm/shared/apps/WUR/ABGC/Annovar/annovar/convert2annovar.pl -format vcf4 tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.var.mpileup.flt.vcf > tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.annov.txt
/cm/shared/apps/WUR/ABGC/Annovar/annovar/annotate_variation.pl -out tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.annov --buildver galGal4 -dbtype ensgene tmpGG-WUR-002/GG-WUR-002_rh.dedup_st.reA.annov.txt /cm/shared/apps/WUR/ABGC/Annovar/annovar/chickendb
