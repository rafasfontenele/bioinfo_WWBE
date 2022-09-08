#!bin/bash
#This bash script will set the folders and files to run Nextclade and Pangolin at the same time through swarm

#variables
date=$1 #Date of GISAID download used to create an output folder for processed filess
dir=${PWD}
base="${dir}/gisaid_process/${date}"

#finish pre-processing
#Download sequences from the redo.list and save as redo_sequences.fasta before running step 2

cd "${dir}/gisaid_process/${date}/01.preprocessing"
module load seqtk seqkit
seqkit replace -j 56 -p "^(.+)\|(.+)\|(.+)$" -r "\${2}" redo_sequences.fasta > sequences.add.fa
cat sequences.rename.fa sequences.add.fa >sequences.acc.fa
seqkit fx2tab sequences.add.fa -n|cat - redo.list |sort |uniq -u > remove.list
module load seqkit
rm sequences.add.fa sequences.rename.fa
seqkit split sequences.acc.fa -s 50000 -O split


# going into the versiond (date)
cd $base

#02.Nextclade
#create folder

mkdir 02.nextclade
cd 02.nextclade

#copy the updated run_nextclade.sh to be used 
cp $dir/scripts/run_nextclade.sh $base/02.nextclade

#run nextclade in split files - creates output file for each split file 
# create swarm file using split X 
for file in $base/01.preprocessing/split/*; 
do
    searchstring="split/"
    file_name=${file#*$searchstring}
    file_num=$(echo $file_name | cut -d "." -f 1,2,3)
    echo "bash run_nextclade.sh "$file_num" "$date >> run.swarm
done

# run nextclade swarm
swarm -f run.swarm -t 32 -g 80 --time 72:00:00

cd ..

#03.Pangolin
# run pangolin

mkdir 03.pangolin
cd 03.pangolin

#copy the updated run_pangolin.sh to be used
cp $dir/scripts/run_pangolin.sh $base/03.pangolin

#this step needs to be done in split files - use the same strategy as nextclade
#create swarm file using split X 
for file in $base/01.preprocessing/split/*; 
do
    searchstring="split/"
    file_name=${file#*$searchstring}
    file_num=$(echo $file_name | cut -d "." -f 1,2,3)
    echo "bash run_pangolin.sh "$file_num" "$date >> run.swarm
done

# run pangolin swarm
swarm -f run.swarm -t 32 -g 80 --time 72:00:00
