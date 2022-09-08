module load samtools/1.15
module load bwa/0.7.17
module load vcflib/1.0.1
module load snpEff/5.0

sample=$1

dir=${PWD}
refs="${dir}/files/"
base="${dir}/samples_process/"


analysis="${base}analysis/"

#aligning
bwa mem -t 32 -v 1 "${refs}"NC_045512.2.fasta "${base}raw_data/${sample}"_1.fastq.gz "${base}raw_data/${sample}"_2.fastq.gz | samtools view -@ 32 -b -f 1 -F 268 -q 20 -s 1.0 | samtools sort -@ 32 -o "${analysis}bamfiles/${sample}".sorted.bam
echo "Sample was initially aligned with BWA"

#trimming primers and base quality
module load ivar/1.3.1
ivar trim -e -m 1 -q 0 -b "${refs}"swift_primersv2.bed -p "${analysis}bamfiles/${sample}".trimmed -i "${analysis}bamfiles/${sample}".sorted.bam
echo "Amplification primers were removed"

#lofreq realign - may need to separate samtools to controls the version
module load lofreq/2.1.5
lofreq viterbi --defqual 2 -f "${refs}"NC_045512.2.fasta "${analysis}bamfiles/${sample}".trimmed.bam | samtools sort -@ 32 -o "${analysis}bamfiles/${sample}".trimmed.realigned.sorted.bam
echo "Bam file was re-aligned using lofreq"

#lofreq indelquality and indexing
lofreq indelqual --dindel -f "${refs}"NC_045512.2.fasta -o "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.bam "${analysis}bamfiles/${sample}".trimmed.realigned.sorted.bam 
echo "Bam file annotated with indel information"

#lofreq variant call #1
lofreq call --min-cov 5 --max-depth 1000000 -q 30 -Q 30 -e --min-mq 20 --sig 0.0005 --bonf dynamic --no-default-filter --no-default-filter -f "${refs}"NC_045512.2.fasta --call-indels -o "${analysis}vcffiles/${sample}"_trimmed.vcf "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.bam
echo "First round of variant call with lofreq finished"

#lofreq filter
lofreq filter -V 0 -v 5 -a 0.05 -A 0.95 -i "${analysis}vcffiles/${sample}"_trimmed.vcf -o "${analysis}vcffiles/${sample}"_trimmed_filtered.vcf
echo "Fisrt round of variant call filtered"

#step to revome primer bias reads
module load ivar/1.3.1
ivar getmasked -i "${analysis}vcffiles/${sample}"_trimmed_filtered.vcf -b "${refs}"swift_primersv2.bed -f "${refs}"swift_primersv2_info.tsv -p "${analysis}${sample}"_primer_mismatchers_indices
python "${dir}/scripts/"completemask.py "${analysis}${sample}"_primer_mismatchers_indices.txt "${refs}"swift_primersv2_info.tsv "${analysis}${sample}"_primer_mismatchers_indices_v2.txt
ivar removereads -i "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.bam -b "${refs}"swift_primersv2.bed -t "${analysis}${sample}"_primer_mismatchers_indices_v2.txt -p "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.readsrem.bam
echo "Bias reads removed with ivar"

# Get depth files to be used by Freyja
samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f "${refs}"NC_045512.2.fasta "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.readsrem.bam | cut -f1-4 > "${analysis}depth/${sample}".depth

#lofreq variant call #2
module load lofreq/2.1.5
lofreq call --min-cov 5 --max-depth 1000000 -q 30 -Q 30 -e --min-mq 20 --sig 0.0005 --bonf dynamic --no-default-filter --no-default-filter -f "${refs}"NC_045512.2.fasta --call-indels -o "${analysis}vcffiles/${sample}"_trimmed_v2.vcf "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.readsrem.bam
echo "Second round of variant call with lofreq finished - using removed reads bam"

#vcf files intersect #1
vcfintersect -r "${refs}"NC_045512.2.fasta -v -w 0  -i "${analysis}vcffiles/${sample}"_trimmed_v2.vcf "${analysis}vcffiles/${sample}"_trimmed_filtered.vcf > "${analysis}vcffiles/${sample}"_trimmed_intersect.vcf
echo "VCF_trimmed_v2 intersected with vcf_trimmed_filtered"

#re-naming
sed -r --sandbox -e 's/^(#CHROM.+)$/##FILTER=<ID=AmpliconRemoval,Description="Variant removed upon removal of amplicon">\n\1/g' -e 's/(.+\t)PASS(\t.+)/\1AmpliconRemoval\2/g' "${analysis}vcffiles/${sample}"_trimmed_intersect.vcf > "${analysis}vcffiles/${sample}"_trimmed_renamed.vcf
echo "Intersected vcf re-named"

#vcf files union
vcfintersect -r "${refs}"NC_045512.2.fasta -w 0 -u "${analysis}vcffiles/${sample}"_trimmed_renamed.vcf  "${analysis}vcffiles/${sample}"_trimmed_v2.vcf > "${analysis}vcffiles/${sample}"_trimmed_union.vcf
lofreq filter -V 0 -v 0 -a 0.0 -A 0.0 -b fdr -c 0.001 --print-all -i "${analysis}vcffiles/${sample}"_trimmed_union.vcf -o "${analysis}vcffiles/${sample}"_trimmed_union_filtered.vcf
echo "Final vcf created - union of vcf_trimmed_renamed with vcf_trimmed_v2"

#annotation
#wget https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_NC_045512.2.zip
#unzip snpEff_v5_0_NC_045512.2.zip will create data folder - This is under files folder. In case you ran into issues follow those steps to downlaod again.
snpEff eff -nodownload -dataDir ${refs}data -i vcf -o vcf -formatEff -classic -no-downstream -no-intergenic -no-upstream -ud 0 -stats "${analysis}vcffiles/${sample}"stats.html -noLog NC_045512.2 "${analysis}vcffiles/${sample}"_trimmed_union_filtered.vcf > "${analysis}vcffiles/${sample}"_trimmed_union_snpEff.vcf
echo "Final VCF annotated"

/data/salgadofontenr2/conda/envs/cov-dist/bin/python "${dir}/scripts/"add_variants_tsv.py  "${analysis}vcffiles/${sample}"_trimmed_union_snpEff.vcf  "${analysis}tsvfiles/${sample}"_trimmed_union_snpEff_final.tsv

module load samtools/1.15
samtools fastq -1 "${analysis}fastqfiles/${sample}"_1.fq.gz -2 "${analysis}fastqfiles/${sample}"_2.fq.gz -@ 32 -n "${analysis}bamfiles/${sample}".trimmed.realigned.indelqual.bam
echo "Generated FASTQ files for SRA submission"

echo "Cleaning up intermidiete files"
rm "${analysis}bamfiles/${1}".sorted.bam
rm "${analysis}bamfiles/${1}".sorted.bam.bai
rm "${analysis}bamfiles/${1}".trimmed.bam
rm "${analysis}bamfiles/${1}".trimmed.realigned.sorted.bam
rm "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.bam
rm "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.bam.bai
rm "${analysis}vcffiles/${1}"_trimmed.vcf
rm "${analysis}vcffiles/${1}"_trimmed_filtered.vcf
rm "${analysis}vcffiles/${1}"_trimmed_v2.vcf
rm "${analysis}vcffiles/${1}"_trimmed_intersect.vcf
rm "${analysis}vcffiles/${1}"_trimmed_union.vcf
rm "${analysis}vcffiles/${1}"_trimmed_union_filtered.vcf
rm "${analysis}${1}"_primer_mismatchers_indices.txt
rm "${analysis}${1}"_primer_mismatchers_indices_v2.txt
rm "${analysis}vcffiles/${1}"_trimmed_renamed.vcf
