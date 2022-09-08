module load datamash

#variables
dir=$(pwd)
date=$1 #Date of GISAID download and the output folder
metadata=$2
base="${dir}/gisaid_process/${date}"

# this merges nextclade and pangolin results regardless of passing quality controls
cat $base/02.nextclade/output/nextclade.tsv | sort -k1,1 > $base/02.nextclade/output/nextclade_sorted.tsv
awk '{print $1"\t"$2}' FS="," $base/03.pangolin/lineage_report.csv |sort -k1,1 > $base/03.pangolin/acc_lineage_sorted.tsv

awk 'NR>1' $base/02.nextclade/output/nextclade_sorted.tsv |sort -k1,1 |join -1 1 -2 1 $base/03.pangolin/acc_lineage_sorted.tsv - -t $'\t'  > $dir/novel_snvs/table_allgenomes.tsv
cat $dir/novel_snvs/table_allgenomes.tsv | sort -k1,1 -u > $dir/novel_snvs/table_allgenomes_nodup.tsv

cat $dir/gisaid_process/raw_data/metadata.tsv | sort -t$'\t' -k3,3 > $dir/novel_snvs/metadata_sorted.tsv

awk 'NR>1' $dir/novel_snvs/table_allgenomes_nodup.tsv | sort -k1,1 | join -1 3 -2 1 $dir/novel_snvs/metadata_sorted.tsv - -t $'\t' > $dir/novel_snvs/table_allgenomes_metadata.tsv

rm $dir/novel_snvs/metadata_sorted.tsv

# Get all the snvs, insertions and deletions associated with all GISAID sequences
cut -f 1,38,39,40 $dir/novel_snvs/table_allgenomes_metadata.tsv |tr "\t" ","|awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}' FS="," > $dir/novel_snvs/acc_nuc_ins_del.tsv
awk -F '\t' '$2 != ""' $dir/novel_snvs/acc_nuc_ins_del.tsv > $dir/novel_snvs/acc_nuc_ins_del_clean.tsv
cat $dir/novel_snvs/acc_nuc_clean.tsv | datamash --header-out --sort groupby 2 count 2 > $dir/novel_snvs/snvs_ins_del_all.tsv
python3 ${dir}/scripts/check_snvs_present.py -f ${dir}/derived_snvs/child_parent_info_${date}_all.pkl -m $metadata -o ${dir}/novel_snvs/novel_snvs_${date}_all -t $dir/novel_snvs/snvs_ins_del_all.tsv
awk 'FNR==1 && NR!=1{next;}{print}' ${dir}/novel_snvs/novel_snvs_${date}_all/*.tsv > ${dir}/novel_snvs/novel_snvs_${date}_all/novel_snvs_gisaid_${date}_ins_del_all.tsv
cat ${dir}/novel_snvs/novel_snvs_${date}_all/novel_snvs_gisaid_${date}_ins_del_all.tsv | awk -F"\t" 'NR==1||$20 == "NP"' > ${dir}/novel_snvs/novel_snvs_${date}_all/novel_snvs_gisaid_${date}_ins_del_all_NPonly.tsv
cat ${dir}/novel_snvs/novel_snvs_${date}_all/novel_snvs_gisaid_${date}_ins_del_all_NPonly.tsv | datamash --header-out --sort groupby 16 count 16 unique 14,18,19,20,23,24 > ${dir}/novel_snvs/novel_snvs_${date}_all/novel_snvs_gisaid_${date}_ins_del_all_NPonly_perSNV.tsv

# Cleaning

rm $dir/novel_snvs/acc_nuc_ins_del.tsv