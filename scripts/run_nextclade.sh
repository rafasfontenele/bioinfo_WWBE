num=$1
version=$2
wd=$(pwd)
module load singularity
singularity run /gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/auxiliar_files/nextclade/nextclade_2.4.0.sif nextclade run \
   --output-all=${wd}/output/ \
   --in-order \
   --input-dataset= gisaid_process/data/sars-cov-2/ \
   --output-selection=tsv,csv \
   --output-basename=${num}.nextclade \
   gisaid_process/${version}/01.processing/split/${num}.fa
sleep 10