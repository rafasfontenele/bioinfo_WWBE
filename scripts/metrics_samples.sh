#Get metrics of the sequencing batch

module load samtools
echo -e "File \tMean \tBreadth \tRead_count" >> /data/salgadofontenr2/ww_cov_2021/metrics_galaxy_${1}.txt

for file in /gpfs/gsfs12/users/Irp-jiang/share/covid_data/WWBE/${1}/analysis_galaxy/bamfiles/*.trimmed.realigned.indelqual.readsrem.bam;
do
	prefix=$(echo $file | cut -d "/" -f 12 | cut -d "." -f 1)
	count=$(samtools view -c -F 4 "${file}")
	mean=$(samtools depth -a "${file}" | awk '{c++;s+=$3}END{print s/c}')
	breadth=$(samtools depth -a "${file}" | awk '{c++; if($3>10) total+=1}END{print (total/(c))*100}')
	echo -e "${prefix} \t${mean} \t${breadth} \t${count}" >> /data/salgadofontenr2/ww_cov_2021/metrics_galaxy_${1}.txt
done