dir=$(pwd)
gisaid_date=$1 # GISAID download date
usher_date=$2 # Usher download date
path="/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/GISAID_processing_WWBE/"
base_gisaid="${dir}/gisaid_process/${gisaid_date}"
base_usher="${dir}/usher_process/${usher_date}"
base="${dir}/derived_snvs/"

cp ${dir}/pango-designation/lineage_notes.txt $base
cp ${dir}/pango-designation/pango_designation/alias_key.json $base
cp ${base_gisaid}/05.sort_SNVs/full_table_info_ins_del_final_gtr75_listsnvs.tsv $base
cp ${base_usher}/lineagePaths_edited_${usher_date}.pkl $base
cp ${base_usher}/lineagePaths_edited_clades_${usher_date}.pkl $base


cat ${base}/lineage_notes.txt | cut -f 1 | sed '/^*/d' > ${base}/lineage_list.txt
python $dir/scripts/get_derivedsnvs.py $gisaid_date ${base}/full_table_info_ins_del_final_gtr75_listsnvs.tsv ${base}/lineagePaths_edited_${usher_date}.pkl