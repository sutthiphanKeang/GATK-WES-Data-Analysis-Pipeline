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
# > mkdir aligned_reads reads scripts results data
# check directory using
# > ls *
# keep this file in scripts folder.



# 2. download sample data file using
# > wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
# > wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz 
# check correctness using
# > ls *
# sample files have to be in reads folder. 



# 3. download reference file using
# > wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# you will get a zipped file. To unzip file using
# > gunzip /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/hg38.fa.gz

# for filter
# wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz
# wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
# wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
# wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz


# 4. index refernce file -> .fai file
# > samtools faidx /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/hg38.fa



# 5. -> .dict file
# > gatk CreateSequenceDictionary R=/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/hg38.fa O=/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/hg38.dict



# 6. download known sites files for BQSR from GATK resource bundle (use in recalibation)
# > wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
# > wget -P /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# ====================================set PATH==================================== #

ref="/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/hg38.fa"
known_sites="/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/aligned_reads"
reads="/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reads"
results="/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/results"
data="/Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/data"

# ====================================alignment==================================== #

# ------1. Run fast qc------ check quility
# output: fastqc.html, fastqc.zip in reads folder.

# fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
# fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/



# ------2. BWA-MEM------ map to reference
# output: .paired.sam file in aligned_reads folder.
# *if use others sample data, you have to change tID, tSM, and .paired.sam name.

# bwa index ${ref}
# bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

# reads .paired.sam file using > samtools view SRR062634.paired.sam | less



# ------3. GATK------ mark duplicates and sort
# output: .bam, .bam.bai, and .bam.sbi file in aligned_reads folder. 

# gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# ------4. GATK------ base quality recalibation
# output: .table and bqsr.bam file 

# gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
# gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam



# ------5. GATK------ collect alignment and insert size metrics
# output: .txt and .pdf file in aligned reads folder.

# gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

# ====================================variant calling==================================== #

# ------1. GATK------ call vairants
# output: .vcf file in result

# gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# ------2. GATK------ snps and indels
# output snps.vcf and indels.vcf

# gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
# gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf

# ------X. GATK------ HaplotypeCaller for chr9
# gatk HaplotypeCaller \
#     -R ${ref} \
#     -I ${aligned_reads}/teacher_AB/gr_AB_KT.final.bam \
#     -L chr9 \
#     -ERC GVCF \
#     -O ${results}/teacher_AB/hete_test/raw_variants_chr9_AB_KT.g.vcf.gz


# ------X. GATK------ GenotypeGVCFs
# gatk GenotypeGVCFs \
#     -R ${ref} \
#     -V ${results}/teacher_AB/raw_variants_chr9_AB_KT.g.vcf.gz \
#     -L chr9 \
#     -O ${results}/teacher_AB/probability_variants_chr9_AB_KT.vcf.gz

#  gatk IndexFeatureFile \
#      -I /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/hapmap_3.3.hg38.vcf.gz 

#  gatk IndexFeatureFile \
#      -I /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/1000G_omni2.5.hg38.vcf.gz

# gatk IndexFeatureFile \
#      -I /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/1000G_phase1.snps.high_confidence.hg38.vcf.gz

# gatk IndexFeatureFile \
#      -I /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/dbsnp_138.hg38.vcf.gz


# ------X. GATK------ Filter and prepare analysis ready variants "Variant Quality Score Recalibration"
# gatk VariantRecalibrator \
#     -R ${ref} \
#     -V ${results}/teacher_AB/probability_variants_chr9_AB_KT.vcf.gz \
#     --max-gaussians 6 \
#     --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/hapmap_3.3.hg38.vcf.gz \
#     --resource:omni,known=false,training=true,truth=true,prior=12.0 /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/1000G_omni2.5.hg38.vcf.gz \
#     --resource:1000G,known=false,training=true,truth=false,prior=10.0 /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#     --resource:dbsnp,known=true,training=false,truth=false,prior=7.0 /Users/sutthiphanprananpaeng/CMU/204490Resarch/demo/reference/filter/dbsnp_138.hg38.vcf.gz \
#     -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
#     --trust-all-polymorphic \
#     -mode BOTH \
#     -O ${results}/teacher_AB/variants_chr9_AB_KT_snps.recal \
#     --tranches-file ${results}/teacher_AB/variants_chr9_AB_KT_snps.tranches

# gatk ApplyVQSR \
#     -R ${ref} \
#     -V ${results}/teacher_AB/probability_variants_chr9_AB_KT.vcf.gz \
#     -O ${results}/teacher_AB/variants_chr9_AB_KT.vqsr.vcf \
#     --truth-sensitivity-filter-level 99.0 \
#     --tranches-file ${results}/teacher_AB/variants_chr9_AB_KT_snps.tranches \
#     --recal-file ${results}/teacher_AB/variants_chr9_AB_KT_snps.recal \
#     -mode BOTH

# gatk VariantFiltration \
#     -R ${ref} \
#     -V /Users/sutthiphanprananpaeng/CMU/204490Resarch/test_map_blood/probability_variants_chr9_AB_KT.vcf  \
#     -O /Users/sutthiphanprananpaeng/CMU/204490Resarch/test_map_blood/variants_chr9_AB_KT.vqsr.varfilter.vcf \
#     # --filter-name "Low_depth10" \
#     # --filter-expression "DP < 10" 
#     --filter-expression "QD < 2.0" --filter-name "QD2" \
#     --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
#     --filter-expression "SOR > 3.0" --filter-name "SOR3" \
#     --filter-expression "FS > 60.0" --filter-name "FS60" \
#     --filter-expression "MQ < 40.0" --filter-name "MQ40" \
#     --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
#     --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" 



gatk VariantsToTable \
    -V ${results}/tearcher/variants_chr9_O_NA.vqsr.varfilter.vcf  \
    -F CHROM -F POS -F TYPE -F REF -F ALT -GF GT \
    -O ${results}/tearcher/variants_chr9_O_NA.vqsr.varfilter.pass.tsv

# ====================================finish==================================== #
