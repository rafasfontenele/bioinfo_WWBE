num=$1
version=$2 #Date of GISAID download
wd=$(pwd)
module load singularity
singularity run /gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/auxiliar_files/pangolin/pangolin_4.1.2-pdata-1.13.sif pangolin --outfile ${num}.lineage_report.csv -t 32 /gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/${version}/01.processing/split/${num}.fa
sleep 10