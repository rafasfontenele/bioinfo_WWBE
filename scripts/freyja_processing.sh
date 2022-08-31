#Freyja processing step 1 - one per sample

#freyja update - run this to update the Usher tree for freyja
filename=$1
wd=$(pwd)

lines=$(cat $filename)
refs="/gpfs/gsfs12/users/Irp-jiang/share/covid_data/WWBE/info/"

for line in $lines;
do
    prefix=$(echo $line | cut -d "/" -f 12 | cut -d "_" -f 1)
    freyja demix ${line}_trimmed_union.vcf depth/${prefix}.depth --output freyja_output/${prefix}_freyja
done