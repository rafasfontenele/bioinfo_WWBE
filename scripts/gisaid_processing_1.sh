#!/bin/bash

dataset=$1 #Date of GISAID download which is the name of the folder containing the fasta and metadata file
step=$2

fasta="/gpfs/gsfs12/users/Irp-jiang/share/covid_data/GISAID/${1}/sequences.fasta"
metadata="/gpfs/gsfs12/users/Irp-jiang/share/covid_data/GISAID/${1}/metadata.tsv"

if [ $step == "1" ]
then
    mkdir "/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/$1"
    mkdir "/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/$1/01.processing"
    cd "/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/$1/01.processing"
    module load seqtk seqkit
    seqkit fx2tab -j 56 $fasta > sequences.tsv
    sort -k1,1 sequences.tsv > sequences_sort.tsv
    awk '{print $1"|"$4"|"$16"\t"$3}' FS="\t" $metadata |sed '1d'|sort -k1,1 -S 80% --parallel 56 > name.tsv
    cut -f1 name.tsv|uniq -d >name.d
    fgrep -f name.d -w name.tsv|cut -f 2 > redo.list
    paste name.tsv sequences_sort.tsv|fgrep -f redo.list -w -v |cut -f 2,4|seqkit -j 56 tab2fx >sequences.rename.fa
fi
if [ $step == "2" ]
then 
    # Download sequences from the redo.list and save as redo_sequences.fasta before running step 2
    cd "/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/$1/01.processing"
    module load seqtk seqkit
    seqkit replace -j 56 -p "^(.+)\|(.+)\|(.+)$" -r "\${2}" redo_sequences.fasta > sequences.add.fa
    cat sequences.rename.fa sequences.add.fa >sequences.acc.fa
    seqkit fx2tab sequences.add.fa -n|cat - redo.list |sort |uniq -u > remove.list
    module load seqkit
    rm sequences.rename.fa sequences.add.fa
    rm -rf sequences.tsv name.tsv name.d sequences_sort.tsv
    seqkit split sequences.acc.fa -s 50000 -O split
fi