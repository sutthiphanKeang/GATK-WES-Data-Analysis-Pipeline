# GATK-WES-Data-Analysis-Pipeline 
This process is part of the framework of the Whole-Exome Sequencing (WES) Analysis for ABO Subgroups Identification research project to analyze DNA and obtain results in the form of a tsv file containing necessary data for subsequent processes, as shown in the research report, specified in the Proposed Framework.

## Installation

### System requirements
* Unix-style OS (Linux, MacOSX, etc)
* Openjdk 8

### Tools
* [GATK version 4.3.0.0](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
* [BWA version 0.7.17](https://bio-bwa.sourceforge.net/)
* [SAMTOOLS version 1.16.1](http://www.htslib.org/download/)
* [FASTQC version 0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Setup

Go to `<PATH OF DIRECTORY>` and create directory
```shell
mkdir aligned_reads reads scripts results data reference
```

Download data file using 
(you can skip this step if you already have input files. It must be contained in the reads folder.)
```shell
wget -P <PATH OF DIRECTORY>/reads [DNA_FORWARD].fastq.gz
wget -P <PATH OF DIRECTORY>/reads [DNA_REVERSE].fastq.gz 
```
Download reference file
```shell
wget -P <PATH OF DIRECTORY>/reference https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip <PATH OF DIRECTORY>/reference/hg38.fa.gz
```

Index reference file
```shell
samtools faidx <PATH OF DIRECTORY>/reference/reference/hg38.fa
```

Create sequence dictionary for refernce file
```shell
gatk CreateSequenceDictionary R=<PATH OF DIRECTORY>/reference/hg38.fa O=[PATH OF FOLDER]/reference/hg38.dict
```

Download reference file for VariantRecalibrator
```shell
wget -P <PATH OF DIRECTORY>/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hapmap_3.3.hg38.vcf.gz
wget -P <PATH OF DIRECTORY>/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz
wget -P <PATH OF DIRECTORY>/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -P <PATH OF DIRECTORY>/reference ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_138.hg38.vcf.gz

gatk IndexFeatureFile -I <PATH OF DIRECTORY>/reference/hapmap_3.3.hg38.vcf.gz 
gatk IndexFeatureFile -I <PATH OF DIRECTORY>/reference/1000G_omni2.5.hg38.vcf.gz
gatk IndexFeatureFile -I <PATH OF DIRECTORY>/reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz
gatk IndexFeatureFile -I <PATH OF DIRECTORY>/reference/dbsnp_138.hg38.vcf.gz
```

## Usage

Run `Analysis_FastQ_To_TSV.sh`
```shell
cd scripts
./Analysis_FastQ_To_TSV.sh
```
## Reference

* [WGS Variant Calling](https://www.youtube.com/watch?v=iHkiQvxyr5c)
* [Variant calling using GATK4](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/)
