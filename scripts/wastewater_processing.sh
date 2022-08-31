#Bash script - processing of wastewater sequencing raw_data

#!/bin/bash

module load samtools/1.15
module load bwa/0.7.17

module load vcflib/1.0.1
module load snpEff/5.0

batch=$2
refs="/gpfs/gsfs12/users/Irp-jiang/share/covid_data/WWBE/info/"
base="/gpfs/gsfs12/users/Irp-jiang/share/covid_data/WWBE/${batch}/"

analysis="${base}analysis_galaxy/"

#aligning
bwa mem -t 32 -v 1 "${refs}"NC_045512.2.fasta "${base}raw_data/${1}"_1.fastq.gz "${base}raw_data/${1}"_2.fastq.gz | samtools view -@ 32 -b -f 1 -F 268 -q 20 -s 1.0 | samtools sort -@ 32 -o "${analysis}bamfiles/${1}".sorted.bam
echo "Sample was initially aligned with BWA"

#trimming primers and base quality
module load ivar/1.3.1
ivar trim -e -m 1 -q 0 -b "${refs}"swift_primersv2.bed -p "${analysis}bamfiles/${1}".trimmed -i "${analysis}bamfiles/${1}".sorted.bam
echo "Amplification primers were removed"

#lofreq realign - may need to separate samtools to controls the version
module load lofreq/2.1.5
lofreq viterbi --defqual 2 -f "${refs}"NC_045512.2.fasta "${analysis}bamfiles/${1}".trimmed.bam | samtools sort -@ 32 -o "${analysis}bamfiles/${1}".trimmed.realigned.sorted.bam
echo "Bam file was re-aligned using lofreq"

#lofreq indelquality and indexing
lofreq indelqual --dindel -f "${refs}"NC_045512.2.fasta -o "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.bam "${analysis}bamfiles/${1}".trimmed.realigned.sorted.bam 
echo "Bam file annotated with indel information"

#lofreq variant call #1
lofreq call --min-cov 5 --max-depth 1000000 -q 30 -Q 30 -e --min-mq 20 --sig 0.0005 --bonf dynamic --no-default-filter --no-default-filter -f "${refs}"NC_045512.2.fasta --call-indels -o "${analysis}vcffiles/${1}"_trimmed.vcf "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.bam
echo "First round of variant call with lofreq finished"

#lofreq call --min-cov 20  -q 20 -Q 20 -e --min-mq 20 -f "${refs}"NC_045512.2.fasta --call-indels -o "${analysis}vcffiles/${1}"_trimmed.vcf "${analysis}vcffiles/${1}".trimmed.realigned.indelqual.bam
lofreq filter -V 0 -v 5 -a 0.05 -A 0.95 -i "${analysis}vcffiles/${1}"_trimmed.vcf -o "${analysis}vcffiles/${1}"_trimmed_filtered.vcf
echo "Fisrt round of variant call filtered"

#step to revome primer bias reads
module load ivar/1.3.1
ivar getmasked -i "${analysis}vcffiles/${1}"_trimmed_filtered.vcf -b "${refs}"swift_primersv2.bed -f "${refs}"swift_primersv2_info.tsv -p "${analysis}${1}"_primer_mismatchers_indices
python completemask.py "${analysis}${1}"_primer_mismatchers_indices.txt "${refs}"swift_primersv2_info.tsv "${analysis}${1}"_primer_mismatchers_indices_v2.txt
ivar removereads -i "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.bam -b "${refs}"swift_primersv2.bed -t "${analysis}${1}"_primer_mismatchers_indices_v2.txt -p "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.readsrem.bam
echo "Bias reads removed with ivar"

# Get depth files to be used by Freyja
samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f "${refs}"NC_045512.2.fasta "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.readsrem.bam | cut -f1-4 > "/gpfs/gsfs12/users/Irp-jiang/share/rafa_data/freyja_data/depth/${1}".depth

#lofreq variant call #2
module load lofreq/2.1.5
lofreq call --min-cov 5 --max-depth 1000000 -q 30 -Q 30 -e --min-mq 20 --sig 0.0005 --bonf dynamic --no-default-filter --no-default-filter -f "${refs}"NC_045512.2.fasta --call-indels -o "${analysis}vcffiles/${1}"_trimmed_v2.vcf "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.readsrem.bam
echo "Second round of variant call with lofreq finished - using removed reads bam"

#vcf files intersect #1
vcfintersect -r "${refs}"NC_045512.2.fasta -v -w 0  -i "${analysis}vcffiles/${1}"_trimmed_v2.vcf "${analysis}vcffiles/${1}"_trimmed_filtered.vcf > "${analysis}vcffiles/${1}"_trimmed_intersect.vcf
echo "VCF_trimmed_v2 intersected with vcf_trimmed_filtered"

#re-naming
sed -r --sandbox -e 's/^(#CHROM.+)$/##FILTER=<ID=AmpliconRemoval,Description="Variant removed upon removal of amplicon">\n\1/g' -e 's/(.+\t)PASS(\t.+)/\1AmpliconRemoval\2/g' "${analysis}vcffiles/${1}"_trimmed_intersect.vcf > "${analysis}vcffiles/${1}"_trimmed_renamed.vcf
echo "Intersected vcf re-named"

#vcf files union
vcfintersect -r "${refs}"NC_045512.2.fasta -w 0 -u "${analysis}vcffiles/${1}"_trimmed_renamed.vcf  "${analysis}vcffiles/${1}"_trimmed_v2.vcf > "${analysis}vcffiles/${1}"_trimmed_union.vcf
lofreq filter -V 0 -v 0 -a 0.0 -A 0.0 -b fdr -c 0.001 --print-all -i "${analysis}vcffiles/${1}"_trimmed_union.vcf -o "${analysis}vcffiles/${1}"_trimmed_union_filtered.vcf
echo "Final vcf created - union of vcf_trimmed_renamed with vcf_trimmed_v2"

#annotation
#wget https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_NC_045512.2.zip
#unzip snpEff_v5_0_NC_045512.2.zip will create data folder
snpEff eff -nodownload -dataDir ${refs}data -i vcf -o vcf -formatEff -classic -no-downstream -no-intergenic -no-upstream -ud 0 -stats stats.html -noLog NC_045512.2 "${analysis}vcffiles/${1}"_trimmed_union_filtered.vcf > "${analysis}vcffiles/${1}"_trimmed_union_snpEff.vcf
echo "Final VCF annotated"

python "${refs}"add_variant_tsv.py  "${analysis}vcffiles/${1}"_trimmed_union_snpEff.vcf  "${analysis}vcffiles/${1}"_trimmed_union_snpEff_final.tsv

module load samtools/1.15
samtools fastq -1 "${analysis}fastqfiles/${1}"_1.fq.gz -2 "${analysis}fastqfiles/${1}"_2.fq.gz -@ 32 -n "${analysis}bamfiles/${1}".trimmed.realigned.indelqual.bam
echo "Generated FASTQ files for SRA submission"