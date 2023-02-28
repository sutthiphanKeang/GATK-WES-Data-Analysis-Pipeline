# ref https://www.youtube.com/watch?v=iHkiQvxyr5c
# ref https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/

# ******** how to use ********
# copy command after arrow to run that command

# using > echo "Hello world" means print("Hello world") to check the continuity of the program. 
# example
# echo "start variant calling"
# > variant calling command
# echo "complete variant calling"

# [PATH] means The address of the folder.

# run this file using ./test.gatk.sh



# ====================================required (to be generated only once) ==================================== #

# 1. create directory using
mkdir aligned_reads reads scripts results data reference

# check directory using
ls *

# keep this file in scripts folder.


# note: At this step you can substitute your inputs for [DNA_FORWARD] and [DNA_REVERSE] in fastq.gz format.
# 2. download data file using ***you can skip this step if you already have inputs data
wget -P [PATH OF FOLDER]/reads [DNA_FORWARD].fastq.gz
wget -P [PATH OF FOLDER]/reads [DNA_REVERSE].fastq.gz 

# check correctness using
ls *
# sample files have to be in reads folder. 



# 3. download reference file using
wget -P [PATH OF FOLDER]/reference https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# you will get a zipped file. To unzip file using
gunzip [PATH OF FOLDER]/reference/hg38.fa.gz



# 4. download reference file for VariantRecalibrator
wget -P [PATH OF FOLDER]/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz
wget -P [PATH OF FOLDER]/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
wget -P [PATH OF FOLDER]/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -P [PATH OF FOLDER]/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz

# creact index for reference file
gatk IndexFeatureFile -I [PATH OF FOLDER]/reference/hapmap_3.3.hg38.vcf.gz 
gatk IndexFeatureFile -I [PATH OF FOLDER]/reference/1000G_omni2.5.hg38.vcf.gz
gatk IndexFeatureFile -I [PATH OF FOLDER]/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz
gatk IndexFeatureFile -I [PATH OF FOLDER]/reference/dbsnp_138.hg38.vcf.gz



# 5. index refernce file -> .fai file
samtools faidx [PATH OF FOLDER]/reference/reference/hg38.fa



# 6. -> .dict file
gatk CreateSequenceDictionary R=[PATH OF FOLDER]/reference/hg38.fa O=[PATH OF FOLDER]/reference/hg38.dict



# 7. download known sites files for BQSR from GATK resource bundle (use in recalibation)
wget -P [PATH OF FOLDER]/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P [PATH OF FOLDER]/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# ====================================set PATH==================================== #

ref="[PATH OF FOLDER]/reference/hg38.fa"
known_sites="[PATH OF FOLDER]/reference/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="[PATH OF FOLDER]/aligned_reads"
reads="[PATH OF FOLDER]/reads"
results="[PATH OF FOLDER]/results"
data="[PATH OF FOLDER]/data"

# ====================================alignment==================================== #

# ------1. Run fast qc------ check quility
# output: fastqc.html, fastqc.zip in reads folder.

fastqc ${reads}/[DNA_FORWARD].fastq.gz -o ${reads}/
fastqc ${reads}/[DNA_REVERSE].fastq.gz -o ${reads}/



# ------2. BWA-MEM------ map to reference
# output: .paired.sam file in aligned_reads folder.
# *if use others sample data, you have to change tID, tSM, and .paired.sam name.
# [DNA_ID] is the id of the dna which you can set yourself  But most of them already have an id with the input.

bwa index ${ref}
bwa mem -t 4 -R "@RG\tID:[DNA_ID]\tPL:ILLUMINA\tSM:[DNA_ID]" ${ref} ${reads}/[DNA_FORWARD].fastq.gz ${reads}/[DNA_REVERSE].fastq.gz > ${aligned_reads}/[DNA_ID].paired.sam



# ------3. GATK------ mark duplicates and sort
# output: .bam, .bam.bai, and .bam.sbi file in aligned_reads folder. 

gatk MarkDuplicatesSpark -I ${aligned_reads}/[DNA_ID].paired.sam -O ${aligned_reads}/[DNA_ID]_sorted_dedup_reads.bam



# ------4. GATK------ base quality recalibation
# output: .table and bqsr.bam file 

gatk BaseRecalibrator -I ${aligned_reads}/[DNA_ID]_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
gatk ApplyBQSR -I ${aligned_reads}/[DNA_ID]_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/[DNA_ID]_sorted_dedup_bqsr_reads.bam



# ------5. GATK------ collect alignment and insert size metrics
# output: .txt and .pdf file in aligned reads folder.

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/[DNA_ID]_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

# ====================================variant calling==================================== #

# ------1. GATK------ call vairants "HaplotypeCaller for chr9"
# output: .g.vcf.gz file in result

gatk HaplotypeCaller \
    -R ${ref} \
    -I ${aligned_reads}/[DNA_ID]_sorted_dedup_bqsr_reads.bam \
    -L chr9 \
    -ERC GVCF \
    -O ${results}raw_variants_chr9.g.vcf.gz



# ------2. GATK------ GenotypeGVCFs
# output: .vcf.gz file in result

gatk GenotypeGVCFs \
    -R ${ref} \
    -V ${results}raw_variants_chr9.g.vcf.gz \
    -L chr9 \
    -O ${results}/probability_variants_chr9.vcf.gz



# ------3. GATK------ Filter and prepare analysis ready variants "Variant Quality Score Recalibration"
# output: .recal and .tranches file in result

gatk VariantRecalibrator \
    -R ${ref} \
    -V ${results}/teacher_AB/probability_variants_chr9_AB_KT.vcf.gz \
    --max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 [PATH OF FOLDER]/reference/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 [PATH OF FOLDER]/reference/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 [PATH OF FOLDER]/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 [PATH OF FOLDER]/reference/dbsnp_138.hg38.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    --trust-all-polymorphic \
    -mode BOTH \
    -O ${results}/teacher_AB/variants_chr9.recal \
    --tranches-file ${results}/teacher_AB/variants_chr9.tranches



# ------4. GATK------ ApplyVQSR
# output: .vqsr.vcf file in result

gatk ApplyVQSR \
    -R ${ref} \
    -V ${results}/probability_variants_chr9.vcf.gz \
    -O ${results}/variants_chr9.vqsr.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file ${results}/teacher_AB/variants_chr9.tranches\
    --recal-file ${results}/teacher_AB/variants_chr9.recal \
    -mode BOTH



# ------4. GATK------ VariantFiltration
# output: .vqsr.varfilter.vcf file in result

gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/variants_chr9.vqsr.vcf \
    -O ${results}/variants_chr9.vqsr.varfilter.vcf \
    --filter-expression "DP < 10" --filter-name "Low_depth10" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 



# ------4. GATK------ VariantsToTable
# output: .vqsr.varfilter.pass.tsv file in result

gatk VariantsToTable \
    -V ${results}/variants_chr9.vqsr.varfilter.vcf \
    -F CHROM -F POS -F TYPE -F REF -F ALT -GF GT \
    -O ${results}/tearcher/variants_chr9.vqsr.varfilter.pass.tsv

# ====================================finish==================================== #
