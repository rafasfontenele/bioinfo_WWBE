#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

sample=$1
output=${2:-"output/"}
primerBed=${3:-"primer_info/swift_primers_v2.bed"}
primerInfo=${4:-"primer_info/swift_primers_v2_info.tsv"}
cpu=${5:-16}


genome=data/NC_045512.2.fasta
mkdir -p ${output}/bamfiles ${output}/vcffiles ${output}/tsvfiles ${output}/fastqfiles ${output}/freyja

##########################################################
#      Align sample reads to SARS-CoV-2 reference        #
##########################################################

bwa mem -t ${cpu} -v 1 ${genome} raw_data/${sample}_1.fastq.gz raw_data/${sample}_2.fastq.gz | samtools view -@ ${cpu} -b -f 1 -F 268 -q 20 -s 1.0 | samtools sort -@ ${cpu} -o ${output}/bamfiles/${sample}.sorted.bam
echo "Sample was initially aligned with BWA"

if [ ! -f ${output}/bamfiles/${sample}.sorted.bam ]; then
    echo -e "\nFile ${sample}.sorted.bam not found! Are there reads in your fastq files?\n"
    exit 0
fi

##########################################################
#   Trim the amplification primers from alignment file   #
##########################################################

ivar trim -e -m 1 -q 0 -b ${primerBed} -p ${output}/bamfiles/${sample}.trimmed -i ${output}/bamfiles/${sample}.sorted.bam
echo "Amplification primers were removed"
rm ${output}/bamfiles/${sample}.sorted.bam
rm ${output}/bamfiles/${sample}.sorted.bam.bai

##########################################################
#    Reaglined trimmed file using lofreq viterbi         #
##########################################################

lofreq viterbi --defqual 2 -f ${genome} ${output}/bamfiles/${sample}.trimmed.bam | samtools sort -@ ${cpu} -o ${output}/bamfiles/${sample}.trimmed.realigned.sorted.bam
echo "Bam file was re-aligned using lofreq"
rm ${output}/bamfiles/${sample}.trimmed.bam

##########################################################
#        Adding indel quality to alignment file          #
##########################################################

lofreq indelqual --dindel -f ${genome} -o ${output}/bamfiles/${sample}.trimmed.realigned.indelqual.bam ${output}/bamfiles/${sample}.trimmed.realigned.sorted.bam 
echo "Bam file annotated with indel information"
rm ${output}/bamfiles/${sample}.trimmed.realigned.sorted.bam

##########################################################
# Initial variant calling - used to identify primer bias #
##########################################################

lofreq call --min-cov 5 --max-depth 1000000 -q 30 -Q 30 -e --min-mq 20 --sig 0.0005 --bonf dynamic --no-default-filter -f ${genome} --call-indels -o ${output}/vcffiles/${sample}_trimmed.vcf ${output}/bamfiles/${sample}.trimmed.realigned.indelqual.bam
echo "First round of variant call with lofreq finished"

if [ -z "$(grep -v '^#' ${output}/vcffiles/${sample}_trimmed.vcf)" ]; then
    echo -e "\nNo variants found in file ${sample}_trimmed.vcf!\n"
    exit 0
fi

lofreq filter -V 0 -v 5 -a 0.05 -A 0.95 -i ${output}/vcffiles/${sample}_trimmed.vcf -o ${output}/vcffiles/${sample}_trimmed_filtered.vcf
echo "Fisrt round of variant call filtered"
rm ${output}/vcffiles/${sample}_trimmed.vcf

##########################################################
#      Removing reads associated with primer bias -      #
#      obtained based on initial variant calling         #
##########################################################

ivar getmasked -i ${output}/vcffiles/${sample}_trimmed_filtered.vcf -b ${primerBed} -f ${primerInfo} -p ${output}/${sample}_primer_mismatchers_indices
python scripts/completemask.py ${output}/${sample}_primer_mismatchers_indices.txt ${primerInfo} ${output}/${sample}_primer_mismatchers_indices_v2.txt
ivar removereads -i ${output}/bamfiles/${sample}.trimmed.realigned.indelqual.bam -b ${primerBed} -t ${output}/${sample}_primer_mismatchers_indices_v2.txt -p ${output}/bamfiles/${sample}.bam
echo "Bias reads removed with ivar"

if [ ! -s ${output}/${sample}_primer_mismatchers_indices.txt ]; then
    echo -e "\nNo primer bias - empty primer_mismatcher_indices file!\n"
fi

rm ${output}/${sample}_primer_mismatchers_indices.txt
rm ${output}/${sample}_primer_mismatchers_indices_v2.txt

##########################################################
#     Final variant call after removal of bias reads     #
##########################################################

lofreq call --min-cov 5 --max-depth 1000000 -q 30 -Q 30 -e --min-mq 20 --sig 0.0005 --bonf dynamic --no-default-filter --no-default-filter -f ${genome} --call-indels -o ${output}/vcffiles/${sample}_trimmed_v2.vcf ${output}/bamfiles/${sample}.bam
echo "Second round of variant call with lofreq finished - using removed reads bam"

if [ -z "$(grep -v '^#' ${output}/vcffiles/${sample}_trimmed_v2.vcf)" ]; then
    echo -e "\nNo variants found after removing bias reads! check file ${sample}_trimmed_v2.vcf\n"
    exit 0
fi

##########################################################
# Intersection between initial and final variant call to #
#    obtain the variants associated with primer bias     #
##########################################################

vcfintersect -r ${genome} -v -w 0  -i ${output}/vcffiles/${sample}_trimmed_filtered.vcf ${output}/vcffiles/${sample}_trimmed_v2.vcf > ${output}/vcffiles/${sample}_trimmed_intersect.vcf
echo "VCF_trimmed_v2 intersected with vcf_trimmed_filtered"

if [ -z "$(grep -v '^#' ${output}/vcffiles/${sample}_trimmed_intersect.vcf)" ]; then
    echo -e "\nNo variants associated with primer bias!\n"
fi

rm ${output}/vcffiles/${sample}_trimmed_filtered.vcf

##########################################################
#   Adding the "AmpliconRemoval" filter tag to identify  #
#   the variants that were associated with primer bias   #
##########################################################

sed -r -e 's/(.+\t)PASS(\t.+)/\1AmpliconRemoval\2/g' ${output}/vcffiles/${sample}_trimmed_intersect.vcf > ${output}/vcffiles/${sample}_trimmed_renamed.vcf
echo "Intersected vcf re-named"
rm ${output}/vcffiles/${sample}_trimmed_intersect.vcf

##########################################################
#   Union of the final variant call with the vcf file    #
#    containing the "AmpliconRemoval" tagged variants    #
##########################################################

vcfintersect -r ${genome} -w 0 -u ${output}/vcffiles/${sample}_trimmed_renamed.vcf ${output}/vcffiles/${sample}_trimmed_v2.vcf > ${output}/vcffiles/${sample}_trimmed_union.vcf
lofreq filter -V 0 -v 0 -a 0.0 -A 0.0 -b fdr -c 0.001 --print-all -i ${output}/vcffiles/${sample}_trimmed_union.vcf -o ${output}/vcffiles/${sample}_trimmed_union_filtered.vcf
echo "Final vcf created - union of vcf_trimmed_renamed with vcf_trimmed_v2"
rm ${output}/vcffiles/${sample}_trimmed_renamed.vcf
rm ${output}/vcffiles/${sample}_trimmed_v2.vcf

##########################################################
#        Annotation of the union file by snpEff          #
##########################################################

#wget https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_NC_045512.2.zip
#unzip snpEff_v5_0_NC_045512.2.zip will create data folder - This is under files folder. In case you ran into issues follow those steps to downlaod again.

snpEffpath=$(readlink -f  $(which snpEff))
if  [ ! -f ${snpEffpath/snpEff}data/NC_045513.2/snpEffectPredictor.bin ]; then
    mkdir -p ${snpEffpath/snpEff}/data/
    cp -r data/NC_045512.2  ${snpEffpath/snpEff}/data/
fi

snpEff eff -nodownload  -i vcf -o vcf -formatEff -classic -no-downstream -no-intergenic -no-upstream -ud 0 -noStats -noLog NC_045512.2 ${output}/vcffiles/${sample}_trimmed_union_filtered.vcf > ${output}/vcffiles/${sample}_trimmed_union_snpEff.vcf
echo "Final VCF annotated"
rm ${output}/vcffiles/${sample}_trimmed_union.vcf
rm ${output}/vcffiles/${sample}_trimmed_union_filtered.vcf

##########################################################
#  Generating an annotated version of final VCF wihtout  #
#           the "AmpliconRemoval" variants               #
##########################################################

sed -r -e 's/^(#CHROM.+)$/##FILTER=<ID=AmpliconRemoval,Description="Variant removed upon removal of amplicon">\n\1/g' ${output}/vcffiles/${sample}_trimmed_union_snpEff.vcf > ${output}/vcffiles/${sample}_trimmed_union_snpEff_renamed.vcf
vcftools --vcf ${output}/vcffiles/${sample}_trimmed_union_snpEff_renamed.vcf --remove-filtered AmpliconRemoval --recode --recode-INFO-all --stdout -c > ${output}/vcffiles/${sample}.vcf
if [ -z "$(grep -v '^#' ${output}/vcffiles/${sample}.vcf)" ]; then
    echo -e "\nNo variants left after removing primers biased variants - check ${sample}.vcf file\n"
    exit 0
fi
rm ${output}/vcffiles/${sample}_trimmed_union_snpEff_renamed.vcf

##########################################################
#  Transforming final VCF file in a TSV format wihtout   #
#           the "AmpliconRemoval" variants               #
##########################################################

python scripts/add_variants_tsv.py ${output}/vcffiles/${sample}_trimmed_union_snpEff.vcf  ${output}/tsvfiles/${sample}.tsv

##########################################################
#      Generating an fastq files from mapped reads       #
##########################################################

samtools fastq -1 ${output}/fastqfiles/${sample}_1.fq.gz -2 ${output}/fastqfiles/${sample}_2.fq.gz -@ ${cpu} -n ${output}/bamfiles/${sample}.trimmed.realigned.indelqual.bam
echo "Generated FASTQ files for SRA submission"
rm ${output}/bamfiles/${sample}.trimmed.realigned.indelqual.bam
rm ${output}/bamfiles/${sample}.trimmed.realigned.indelqual.bam.bai

##########################################################
#  Run Freyja on the sample using the final VCF files    #
#       wihtout the "AmpliconRemoval" variants           #
##########################################################

samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f ${genome} ${output}/bamfiles/${sample}.bam | cut -f1-4 > ${output}/vcffiles/${sample}.depth
mkdir -p ${output}/freyja
freyja demix ${output}/vcffiles/${sample}.vcf ${output}/vcffiles/${sample}.depth --output ${output}/freyja/${sample}.freyja
echo "Freyja analysis on individual sample"

##################################################
