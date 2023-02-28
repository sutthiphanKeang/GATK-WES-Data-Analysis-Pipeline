# GATK-WES-Data-Analysis-Pipeline (*NOT AVAILABLE)
Project Description

## Installation

### System requirements
* Unix-style OS (Linux, MacOSX, etc)
* Java 1.8

### Tools
* [GATK version 4.3.0.0](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
* [BWA version 0.7.17](https://bio-bwa.sourceforge.net/)
* [SAMTOOLS version 1.16.1](http://www.htslib.org/download/)
* [FASTQC version 0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## Setup

Create directory
```shell
cd $HOME
mkdir aligned_reads reads scripts results data reference
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
