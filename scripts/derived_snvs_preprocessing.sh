wd=$(pwd)
version=$1 # GISAID download date
usher_date=$2 # Usher download date
path="/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/"

#update lineages
cd /data/salgadofontenr2/pango-designation/ #repo that cotains the curated information on pangolin lineages
git pull
cd  $wd

cp /data/salgadofontenr2/pango-designation/lineage_notes.txt $path$version/06.lineage_assingment
cp /data/salgadofontenr2/pango-designation/pango_designation/alias_key.json $path$version/06.lineage_assingment
cp /gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/${version}/05.sort_SNVs/full_table_info_ins_del_final_gtr75_listsnvs.tsv $path$version/06.lineage_assingment
cp /gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/Usher_processing/${usher_date}/lineagePaths_edited_${usher_date}.pkl $path$version/06.lineage_assingment
cp /gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/Usher_processing/${usher_date}/lineagePaths_edited_clades_${usher_date}.pkl $path$version/06.lineage_assingment

# conda activate cov-dist

cd $path$version/06.lineage_assingment
cat lineage_notes.txt | cut -f 1 | sed '/^*/d' > lineage_list.txt
python $wd/get_derivedsnvs.py $version full_table_info_ins_del_final_gtr75_listsnvs.tsv lineagePaths_edited_${usher_date}.pkl