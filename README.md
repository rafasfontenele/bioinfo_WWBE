# bioinfo_WWBE
Bioinformatic pipeline for WWBE project

This repository contains the codes use to process and analyze wastewater samples for the presence of SARS-CoV-2. 

1. WWBE_sample_processing

Process each batch of samples
input: Text file containing the name of the samples
relevant outputs: bam file, vcf edited file, depth file, fastq file (mapped reads only)


2.Process GISAID data - updated frequently
input: metadata and fasta sequences downloaded from GISAID
relevant outputs: nextclade table (nextclade.tsv), pangolin table (lineage_report.csv), merged nextclade/pangolin that passed QC (full_table_nodup.tsv), table with insertions and deletion (full_table_info_ins_del_final_gtr75_listsnvs.tsv).

3.Pocess Usher data - update frequently
relevant outputs:lineagePaths_edited_{date}.pkl, lineagePaths_edited_clade_{date}.pkl, lineage_tree.nwk, pcoa_snvs_formatted_{date}_updated.tsv, metadata_voc_snvs_{date}_updated.tsv

4.Obtaing derived SNVs for analysis
input:full_table_info_ins_del_final_gtr75_listsnvs.tsv, lineagePaths_edited_{date}.pkl
relevant ouputs: child_parent_info_{date}_all.pkl

5.Process freyja - also dependent on Usher tree - updated frequently
input: vcf file ({sample}_trimmed_union.vcf) and depth files from step 1
relevant output: edited output tsv file and figures (stacked bar plots per location)

6.figures - this notebook contains the scripts to generate all figures used in the manuscript (except Freyja)
inputs and outputs are highlighted in the notebook for each figure


Dependencies:
samtools/1.15
bwa/0.7.17
vcflib/1.0.1
snpEff/5.0
ivar/1.3.1
lofreq/2.1.5
matUtils
ete3
freyja
nexclade
pangolin
datamash
pandas
