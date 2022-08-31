#!bin/bash
#This bash script will set the folders and files to run Nextclade and Pangolin at the same time through swarm
#variables
wd=$(pwd)
version=$1 #Date of GISAID download used to create an output folder for processed filess
path="/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/"

# going into the versiond (date)
cd $path$version

#02.Nextclade
#create folder

mkdir 02.nextclade
cd 02.nextclade

#copy the updated run_nextclade.sh to be used 
cp $wd/run_nextclade.sh $path$version/02.nextclade

#run nextclade in split files - creates output file for each split file 
# create swarm file using split X 
for file in $path$version/01.processing/split/*; 
do
    base=$(echo $file | cut -d "/" -f 12 | cut -d "." -f 1,2,3)
    echo "bash run_nextclade.sh "$base" "$version >> run.swarm
done

# run nextclade swarm
swarm -f run.swarm -t 32 -g 80

cd ..

#03.Pangolin
# run pangolin

mkdir 03.pangolin
cd 03.pangolin

#copy the updated run_pangolin.sh to be used
cp $wd/run_pangolin.sh $path$version/03.pangolin

#this step needs to be done in split files - use the same strategy as nextclade
#create swarm file using split X 
for file in $path$version/01.processing/split/*; 
do
    base=$(echo $file | cut -d "/" -f 12 | cut -d "." -f 1,2,3)
    echo "bash run_pangolin.sh "$base" "$version >> run.swarm
done

# run pangolin swarm
swarm -f run.swarm -t 32 -g 80