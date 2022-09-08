#!/bin/bash

date=$1 #Date of GISAID download which is the name of the folder containing the fasta and metadata file

dir=${PWD}
base="${dir}/gisaid_process/raw_data"

fasta="${base}/sequences.fasta"
metadata="${base}/metadata.tsv"
mkdir -p "${dir}/gisaid_process/${date}"
mkdir -p "${dir}/gisaid_process/${date}/01.preprocessing"

cd "${dir}/gisaid_process/${date}/01.preprocessing"
module load seqtk seqkit
seqkit fx2tab -j 56 $fasta > sequences.tsv
sort -k1,1 sequences.tsv > sequences_sort.tsv
awk '{print $1"|"$4"|"$16"\t"$3}' FS="\t" $metadata |sed '1d'|sort -k1,1 -S 80% --parallel 56 > name.tsv
cut -f1 name.tsv|uniq -d >name.d
fgrep -f name.d -w name.tsv|cut -f 2 > redo.list
paste name.tsv sequences_sort.tsv|fgrep -f redo.list -w -v |cut -f 2,4|seqkit -j 56 tab2fx >sequences.rename.fa
rm -rf sequences.tsv name.tsv name.d sequences_sort.tsv



