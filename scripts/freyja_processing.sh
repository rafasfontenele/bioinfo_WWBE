#Freyja processing step 1 - one per sample

#freyja update - run this to update the Usher tree for freyja
filename=$1
dir=$(pwd)

lines=$(cat $filename)
base="${dir}/samples_process/analysis/"

for line in $lines;
do
    prefix=$(echo $line | cut -d "/" -f 12 | cut -d "_" -f 1)
    freyja demix ${base}vcffiles/${line}_trimmed_union_snpEff.vcf ${base}depth/${prefix}.depth --output ${dir}/freyja_process/${prefix}_freyja
done