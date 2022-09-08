# Bioinformatics Framework for Wastewater-based Surveillance of SARS-CoV-2

### Raw sequencing data processing 

Dependencies:
* [iVAR](https://github.com/andersen-lab/ivar)
* [snpEff](http://pcingola.github.io/SnpEff/download/)
* [lofreq](https://csb5.github.io/lofreq/installation/)
* [bwa](https://github.com/lh3/bwa)
* [vcflib](https://github.com/vcflib/vcflib#INSTALL)
* [samtools](http://www.htslib.org/)
* Python modules: pandas, vcf, sys

1. Environment Setup
```{bash}
module load samtools/1.15
module load bwa/0.7.17
module load vcflib/1.0.1
module load snpEff/5.0
module load ivar/1.3.1
module load lofreq/2.1.5
sed version 4.8

base="./samples_process/"
analysis=${base}analysis/

mdkir -p ${base}raw_data
mkdir -p $analysis
mkdir -p "${analysis}tsvfiles"
mkdir -p "${analysis}vcffiles"
mkdir -p "${analysis}bamfiles"
mkdir -p "${analysis}depth"
mkdir -p "${analysis}fastqfiles"
```

2. Individual sample processing

Place paired-end fastq files in the folder `samples_process/raw_data/`

important: The script is set up for raw data files names `<sample_name>_1.fastq.gz`

Process each sample with the following command line:

This script is set up for samples sequenced with Swift amplicon panel. For sequencing data with Artic protocol prepare primer files as the ones listed in the files folder for Swift.

```{bash}
/bin/bash ./scripts/sample_processing.sh <sample_name>
```

This will output the following files

* `analysis/bamfiles/<sample_name>.trimmed.realigned.indelqual.readsrem.bam` - final alingment bam file
* `analysis/depth/<sample_name>.depth` - depth file
* `analysis/vcffiles/<sample_name>_trimmed_union_snpEff.vcf` - final annotated VCF files
* `analysis/tsvfiles/<sample_name>_trimmed_union_snpEff_final.tsv` - tsv file formatted from final VCF files with variants associated with primer bias removed
* `analysis/fastqfiles/<sample_name>_1.fq.gz` - fastq files generated from the final alingment bam file - only aligned reads
* `analysis/fastqfiles/<sample_name>_2.fq.gz`

### GISAID data processing

Dependecies:

* [seqtk](https://github.com/lh3/seqtk)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/user/nextclade-cli.html)
* [Pangolin](https://github.com/cov-lineages/pangolin)
* [datamash](https://www.gnu.org/software/datamash/)
* Python modules: pandas, sys

1. Environment Setup
```{bash}
module load seqtk seqkit
module load datamash

mkdir -p gisaid_process
base="./gisaid_process/"
mkdir -p ${base}raw_data/
```

2. Processing

Place downloaded GISAID `sequence.fasta` and `metadata.tsv` files in the folder `gisaid_process/raw_data/`

Process GISAID data with the following scripts:

..1.Pre-processing
```{bash}
/bin/bash ./scripts/gisaid_processing_1.sh <date_GISAID_download>
```
After this step, go to GISAID and download the sequences present in the `gisaid_process/<date_GISAID_download>/01.preprocessing/redo.list` output file and save them as `gisaid_process/<date_GISAID_download>/01.preprocessing/redo_sequences.fasta`

..2. Running Nextclade and Pangolin - scripts setup for docker versions and to run using singularity
```{bash}
/bin/bash ./scripts/gisaid_processing_2.sh <date_GISAID_download>
```

..3.Processing Nexclade and Pangolin outputs to obtain insertions, deletions and SNVs associated with each pangolin lineages based on a frequency threshold (used 75.0 (75%))
```{bash}
/bin/bash ./scripts/gisaid_processing_3.sh <date_GISAID_download> <frequency_cuttoff>
```

This will output the following files

* `gisaid_process/<date_GISAID_download>/02.nextclade/nextclade.selected.tsv` - output of nextclade for all GISAID sequences
* `gisaid_process/<date_GISAID_download>/03.pangolin/lineage_report.csv` - output of pangoling for all GISAID sequences
* `gisaid_process/<date_GISAID_download>/04.refine/full_table_nodup.tsv` - merged output from nextclade and pangolin (only sequences that passed quality control)
* `gisaid_process/<date_GISAID_download>/05.sort_SNVs/full_table_info.tsv` - table with info for each GISAID sequence ID: pangolin lineage, SNPs, insertions and deletions
* `gisaid_process/<date_GISAID_download>/05.sort_SNVs/full_table_ins_del_final.tsv` - table containing insertion and deletions with their associated prevalence in GISAID data
* `gisaid_process/<date_GISAID_download>/05.sort_SNVs/full_table_ins_del_final_gtr<frequency_cutoff>.tsv` - table containing insertion and deletions that are present in more then <frequency_cuttoff> of GISAID sequences
* `gisaid_process/<date_GISAID_download>/05.sort_SNVs/full_table_ins_del_final_gtr<frequency_cutoff>_listsnvs.tsv` - insertion and deletions that are present in more then <frequency_cuttoff> of GISAID sequences listed by pangolin lineage

### Usher tree data processing

Dependecies:

* [matUtils](https://usher-wiki.readthedocs.io/en/latest/matUtils.html)
* [ete3](http://etetoolkit.org/download/)
* Python modules: pandas, sys, os, subprocess, urllib.request


1. Environment Setup
```{bash}
mkdir -p usher_process
```

2. Processing

Process Usher data with the following scripts:

* Every time this script is ran it will download the latest Usher tree so use the date to keep track of versions
```{bash}
python ./scripts/usher_processing.py <date_Usher_processing>
```

This will output the following files

* `usher_process/<date_Usher_processing>/all_Paths.txt` - SNVs associated with each node of tree
* `usher_process/<date_Usher_processing>/lineage_nodes` - list of all nodes in the tree
* `usher_process/<date_Usher_processing>/lineage_tree.nwk` - nwk tree containing only parental nodes of pangolin lineages
* `usher_process/<date_Usher_processing>/lineagePaths_edited_<date_Usher_processing>.pkl` - `lineagePaths.txt` file edited to include SNVs in a list format without nextrain clade info
* `usher_process/<date_Usher_processing>/lineagePaths_edited_<date_Usher_processing>.tsv` - `lineagePaths.txt` file edited to include SNVs in a list format without nextrain clade info
* `usher_process/<date_Usher_processing>/lineagePaths_edited_clades_<date_Usher_processing>.pkl` - `lineagePaths.txt` file edited to include SNVs in a list format
* `usher_process/<date_Usher_processing>/lineagePaths_edited_clades_<date_Usher_processing>.tsv` - `lineagePaths.txt` file edited to include SNVs in a list format
* `usher_process/<date_Usher_processing>/lineagePaths.txt` - path of SNVs associated with each pangolin lineage and nextrain clade
* `usher_process/<date_Usher_processing>/node_snvs_all.tsv` - formatted `all_Paths.txt` file
* `usher_process/<date_Usher_processing>/Usher_fulltree.nwk` - downloaded usher tree in nwk format

### Obtaining derived SNVs per lineage

This step is dependent on the GISAID and Usher processing.

Dependecies:

* [Pango designation repo](https://github.com/cov-lineages/pango-designation)
* [ete3](http://etetoolkit.org/download/)
* Python modules: pandas, sys, pickle, numpy.lib.utils, ete3, json


1. Environment Setup
```{bash}
mkdir -p derived_snvs

#Get Pango designation Repo
git clone https://github.com/cov-lineages/pango-designation.git
```

2. Processing

Process derived snvs with the following script:

* This bash scripts will clone the pango-designation repo to get lineage and aliases information and process it with a in house python script
```{bash}
/bin/bash ./scripts/derived_snvs_preprocessing.sh <date_GISAID_download> <date_Usher_processing>
```

This will output the following files

Table containing all derived SNVs associated with each pangolin lineage (tsv and pkl formats):
* `derived_snvs/child_parent_info_<date_GISAID_download>_all.pkl` 
* `derived_snvs/child_parent_info_<date_GISAID_download>_all.tsv` 

### Running Freyja

Dependecies:

* [Freyja](https://github.com/andersen-lab/Freyja)

1. Environment Setup
```{bash}
mkdir -p freyja_process
```

2. Processing

Process each sample with Freyja:

* The samples files will contain one sample name per line. Example in files folder (`freyja_processing_input.txt`).

```{bash}
/bin/bash ./scripts/freyja_processing.sh <samples file>
```
* After running alls samples aggregate the date using the freyja command

```{bash}
freyja aggregate freyja_process/ --output freyja_process/[aggregated-filename.tsv]
```
This will output the following files

* `freyja_process/<sample_name>.freyja` - freyja output for each sample
* `freyja_process/[aggregated-filename.tsv]` - freyja aggregated output


### Population diversity analysis with Cov-Dist

Dependecies:

* [Cov-Dist](https://github.com/ncbi/CoV-Dist)
* [ete3](http://etetoolkit.org/download/)
* Python modules: pandas, sys, pickle, numpy.lib.utils, ete3, json


1. Environment Setup
```{bash}
mkdir -p cov-dist_process


git clone https://github.com/nlm-irp-jianglab/CoV-Dist.git

# set conda environment
conda create --name cov-dist python=3.7 biopython seaborn scikit-bio pysam
conda activate cov-dist
```

2. Processing

Obtain SNPs from VOC (Omicron and Delta):

```{bash}
python ./scripts/voc_snvs_covdist.py <date_Usher_processing>
```

This will output the following files

* `cov-dist_process/metadata_voc_snvs_<date_Usher_download>.tsv` - VOC metadata to be merged with you sample metadata (maintain same collumn structure) to be inputed into cov-dist
* `cov-dist_process/pcoa_snvs_formatted_<date_Usher_download>.tsv` - table with VOC associated SNVs to be inputed into cov-dist

Running Cov-Dist

```{bash}
usage: covdist.py [-h] -f  [-q] [-c] [-r] [-t] -o

Compute pairwise metagenome distance for SARS-CoV-2 samples

optional arguments:
  -h, --help       show this help message and exit
  -f , --file      file that lists the path to sorted bam files (required)
  -q , --qual      Only reads with mapping quality equal to or greater than
                   this value will be counted (0 by default).
  -c , --cov       Only samples with reads mapped to equal to or greater than
                   this fraction of the genome will be used for PcoA analysis
                   (0.5 by default).
  -v , --voc       table with SNVs associated with VOCs of interest (optional)
  -r , --ref       input reference file (optional)
  -t , --threads   number of threads used for computing (optional)
  -o , --out       output folder name (required)
```
Run cov-dist with the following command:

Examples uder files folder for -f `bamfiles_path.txt` and -m `metadata_covdist.tsv`

```{bash}
python3 Cov-Dist/covdist.py -f files/<file_path_bam_files> -t 72 -o cov-dit_process/<output_folder_name> -m files/<metadata_format_cov-dist> -v cov-dist_process/pcoa_snvs_formatted_<date_Usher_download>.tsv
```

### Analysis to identifying novel SNVs 

Dependecies:

* Python modules: pandas, os


1. Environment Setup
```{bash}
mkdir -p novel_snvs
```

2. Processing

To obtain novel SNVs information ran the following script:

Important: this analysis is dependent on the GISAID processing and derivided SNVs processing.

* The bash script will obtain all the SNVs from sequences available in GISAID, associate with the metadata and compare to input samples. A python code is called within the bash script. 

* The metatada file needs to have at least a collum with Sample_name but it should also include if possible City, State, Sampler_Start_Date columns. See example in files folder `WWBE_metadata.tsv`

```{bash}
/bin/bash ./scripts/gisaid_meta_novelsnvs.sh <date_GISAID_processing> <path_to_metadata_file>
```

This will output the following files

* `novel_snvs_${date}_all/*.tsv` - folder with the output for each sample
* `novel_snvs/novel_snvs_gisaid_<date_GISAID_download>.tsv` - output of all SNVs in the samples analyzed - Variant_info column (P) present in GISAID and (NP) not present in GISAID
* `novel_snvs/novel_snvs_gisaid_<date_GISAID_download>_NPonly.tsv` - output with only novel snvs (NP) of the analyzed samples
