module load datamash

# variables
dir=$(pwd)
date=$1 #Date of GISAID download and the output folder
frequency=$2 #frequency of the insertions, deletions and SNVs
base="${dir}/gisaid_process/${date}"


# organize log files
mkdir $base/02.nextclade/log
mv $base/02.nextclade/swarm* $base/02.nextclade/log/

mkdir $base/03.pangolin/log
mv $base/03.pangolin/swarm* $base/03.pangolin/log/

# merge files from Nextclade and Pangolin output

# merge and select Nextclade
awk 'FNR==1 && NR!=1{next;}{print}' $base/02.nextclade/output/*.tsv > $base/02.nextclade/output/nextclade.tsv

#get the samples that passed filters
dos2unix $base/02.nextclade/output/nextclade.tsv

# This was valid before the Nextclade v.1.10.2
#awk 'NR==1||$4=="good" && $28=="good" && $32=="good" && $37=="good" && $41=="good" && $48=="good" && $52=="good"' FS=$'\t' $base/02.nextclade/output/nextclade.tsv > $base/02.nextclade/nextclade.selected.tsv

# new collumns after update to nextclade v.1.11.0
awk 'NR==1||$5=="good" && $38=="good" && $42=="good" && $47=="good" && $51=="good" && $58=="good" && $62=="good"' FS=$'\t' $base/02.nextclade/output/nextclade.tsv > $base/02.nextclade/nextclade.selected.tsv

# Concatanate de pangolin results:
awk '(NR == 1) || (FNR > 1)' $base/03.pangolin/*.csv > $base/03.pangolin/lineage_report.csv

# get the info from samples that passed filter - OLD PANGOLIN - before May 2022 updates
#awk '$12=="passed_qc"{print $1"\t"$2}' FS="," $base/03.pangolin/lineage_report.csv |sort -k1,1 > $base/03.pangolin/acc_pangolin.tsv

# get the info from samples that passed filter
awk '$14=="pass"{print $1"\t"$2}' FS="," $base/03.pangolin/lineage_report.csv |sort -k1,1 > $base/03.pangolin/acc_pangolin.tsv

#04.refine

mkdir $base/04.refine

awk 'NR>1' $base/02.nextclade/nextclade.selected.tsv |sort -k1,1 |join -1 1 -2 1 $base/03.pangolin/acc_pangolin.tsv - -t $'\t'  > $base/04.refine/full_table.tsv
cat $base/04.refine/full_table.tsv | sort -k1,1 -u > $base/04.refine/full_table_nodup.tsv


#05.sort_SNVs

mkdir $base/05.sort_SNVs


echo -e "sequence\tlineage\tsubstitutions\tdeletions\tinsertions\tprivateNucMutations.reversionSubstitutions\tprivateNucMutations.labeledSubstitutions\tprivateNucMutations.unlabeledSubstitutions" > $base/05.sort_SNVs/full_table_info.tsv
cat $base/04.refine/full_table_nodup.tsv | cut -f 1,2,17,18,19,20,21,22 >> $base/05.sort_SNVs/full_table_info.tsv

#take out any empy SNPs lines
#add one nucleotide per line
#get snps info
awk '{print $1,$2","$3}' FS="\t" OFS="\t" $base/05.sort_SNVs/full_table_info.tsv | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}' FS="," OFS="\t" > $base/05.sort_SNVs/full_table_info_nuc.tsv
awk -F '\t' '$3 != ""' $base/05.sort_SNVs/full_table_info_nuc.tsv > $base/05.sort_SNVs/full_table_info_nuc_nonNAN.tsv
#get deletion info
awk '{print $1,$2","$4}' FS="\t" OFS="\t" $base/05.sort_SNVs/full_table_info.tsv | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}' FS="," OFS="\t" > $base/05.sort_SNVs/full_table_info_del.tsv
awk -F '\t' '$3 != ""' $base/05.sort_SNVs/full_table_info_del.tsv > $base/05.sort_SNVs/full_table_info_del_nonNAN.tsv
#get insertion info
awk '{print $1,$2","$5}' FS="\t" OFS="\t" $base/05.sort_SNVs/full_table_info.tsv | awk '{for(i=2;i<=NF;i++){print $1"\t"$i}}' FS="," OFS="\t" > $base/05.sort_SNVs/full_table_info_ins.tsv
awk -F '\t' '$3 != ""' $base/05.sort_SNVs/full_table_info_ins.tsv > $base/05.sort_SNVs/full_table_info_ins_nonNAN.tsv

#merge snps, deltetions and insertios
cat $base/05.sort_SNVs/full_table_info_nuc_nonNAN.tsv $base/05.sort_SNVs/full_table_info_del_nonNAN.tsv $base/05.sort_SNVs/full_table_info_ins_nonNAN.tsv > $base/05.sort_SNVs/full_table_info_all.tsv

#merge deletions and insertion to use SNVS from USHER tree
cat $base/05.sort_SNVs/full_table_info_del_nonNAN.tsv $base/05.sort_SNVs/full_table_info_ins_nonNAN.tsv > $base/05.sort_SNVs/full_table_info_ins_del.tsv

#get count of variants per lineage
cat $base/05.sort_SNVs/full_table_info_all.tsv | datamash --header-out --sort groupby 3,2 count 1 > $base/05.sort_SNVs/full_table_info_all_permut.tsv
cat $base/05.sort_SNVs/full_table_info_ins_del.tsv | datamash --header-out --sort groupby 3,2 count 1 > $base/05.sort_SNVs/full_table_info_ins_del_permut.tsv

#get count of total sequences per lineage
cat $base/05.sort_SNVs/full_table_info.tsv | datamash --header-out --sort groupby 2 countunique 1 > $base/05.sort_SNVs/full_table_info_countseqs.tsv

#get lineage prevalence using above files
python $dir/scripts/get_prevalence.py $base/05.sort_SNVs/full_table_info_countseqs.tsv $base/05.sort_SNVs/full_table_info_all_permut.tsv $base/05.sort_SNVs/full_table_info_all_final.tsv
python $dir/scripts/get_prevalence.py $base/05.sort_SNVs/full_table_info_countseqs.tsv $base/05.sort_SNVs/full_table_info_ins_del_permut.tsv $base/05.sort_SNVs/full_table_info_ins_del_final.tsv

#filter data based on frequency within the lineage - only did it for snvs info since deletions and isertion are already rare
cat $base/05.sort_SNVs/full_table_info_all_final.tsv | awk -v start="${frequency}" -v end="100.0" -F"\t" '$5>=start && $5<=end' > $base/05.sort_SNVs/full_table_info_all_final_gtr${frequency:0:2}.tsv
cat $base/05.sort_SNVs/full_table_info_ins_del_final.tsv | awk -v start="${frequency}" -v end="100.0" -F"\t" '$5>=start && $5<=end' > $base/05.sort_SNVs/full_table_info_ins_del_final_gtr${frequency:0:2}.tsv

#put lineage SNPS per line
cat $base/05.sort_SNVs/full_table_info_all_final_gtr${frequency:0:2}.tsv | datamash --sort groupby 1 count 2 collapse 2 > $base/05.sort_SNVs/full_table_info_all_final_gtr${frequency:0:2}_perlin.tsv
cat $base/05.sort_SNVs/full_table_info_all_final_gtr${frequency:0:2}.tsv | datamash --sort groupby 2 count 1 collapse 1 > $base/05.sort_SNVs/full_table_info_all_final_gtr${frequency:0:2}_listsnvs.tsv
cat $base/05.sort_SNVs/full_table_info_ins_del_final_gtr${frequency:0:2}.tsv | datamash --sort groupby 1 count 2 collapse 2 > $base/05.sort_SNVs/full_table_info_ins_del_final_gtr${frequency:0:2}_perlin.tsv
cat $base/05.sort_SNVs/full_table_info_ins_del_final_gtr${frequency:0:2}.tsv | datamash --sort groupby 2 count 1 collapse 1 > $base/05.sort_SNVs/full_table_info_ins_del_final_gtr${frequency:0:2}_listsnvs.tsv


#Cleaning up
rm $base/03.pangolin/acc_pangolin.tsv
rm $base/04.refine/full_table.tsv
rm $base/05.sort_SNVs/full_table_info_nuc.tsv
rm $base/05.sort_SNVs/full_table_info_nuc_nonNAN.tsv
rm $base/05.sort_SNVs/full_table_info_del.tsv
rm $base/05.sort_SNVs/full_table_info_del_nonNAN.tsv
rm $base/05.sort_SNVs/full_table_info_ins.tsv
rm $base/05.sort_SNVs/full_table_info_ins_nonNAN.tsv
rm $base/05.sort_SNVs/full_table_info_countseqs.tsv
rm $base/05.sort_SNVs/full_table_info_ins_del_permut.tsv
rm $base/05.sort_SNVs/full_table_info_all_permut.tsv