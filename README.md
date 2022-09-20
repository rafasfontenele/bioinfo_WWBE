# Bioinformatics Framework for Wastewater-based Surveillance of SARS-CoV-2

## 1. Environment Setup
### Option 1: Install all dependencies using conda

```
git clone https://github.com/nlm-irp-jianglab/bioinfo_WWBE.git
cd bioinfo_WWBE
conda env update -f wwbe-pipeline/environment.yml --prune
conda activate wwbe
```

### Option 2: Using docker image

```
git clone https://github.com/nlm-irp-jianglab/bioinfo_WWBE.git
cd bioinfo_WWBE
docker build -t wwbe .
```

## 2. Individual sample processing

Place paired-end fastq files in the folder `raw_data/`

important: The script is set up for raw data files names `<sample_name>_1.fastq.gz`

Process each sample with the following command line:

This script is set up for samples sequenced with Swift amplicon panel. For sequencing data with Artic protocol prepare primer files as the ones listed in the files folder for Swift.


### Option 1: Using conda environment
```
/bin/bash ./scripts/process_sample.sh <sample_name>
```

### Option 2: Using docker image
```
docker run -v $(pwd):/workdir -it wwbe ./scripts/process_sample.sh <sample_name>
```

This will output the following files

* `<output>/bamfiles/<sample_name>.trimmed.realigned.indelqual.readsrem.bam` - final alingment bam file
* `<output>/vcffiles/<sample_name>.vcf` - final annotated VCF files
* `<output>/vcffiles/<sample_name>.depth` - depth file
* `<output>/tsvfiles/<sample_name>_trimmed_union_snpEff_final.tsv` - tsv file formatted from final VCF files with variants associated with primer bias removed
* `<output>/fastqfiles/<sample_name>_1.fq.gz` - fastq files generated from the final alingment bam file - only aligned reads
* `<output>/fastqfiles/<sample_name>_2.fq.gz`
* `<output>/freyja/<sample_name>.freyja`

## 3. Recover relative lineage abundances using [Freyja](https://github.com/andersen-lab/Freyja)

After running alls samples aggregate the date using the freyja command

```
freyja aggregate freyja/ --output [aggregated-filename.tsv]
```
This will output the following files
* `[aggregated-filename.tsv]` - freyja aggregated output

## 4. Ordination analysis with [CoV-Dist](https://github.com/nlm-irp-jianglab/CoV-Dist)

