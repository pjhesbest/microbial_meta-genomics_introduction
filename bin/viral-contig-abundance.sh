#!/bin/bash

# Requirements:
#	seqkit https://bioinf.shenwei.me/seqkit/
#	samtools https://www.htslib.org/
# 	minimap2 https://github.com/lh3/minimap2

##########################################################################################################################################
##########################################################################################################################################


echo -e "\e[34mCOMMAND LINE JOB SUBMISSSION:\n\tbash mapping_and_counts.sh $@\e[0m"
echo -e "\e[32m\nThis script performs the following:\n\n\t1. Filters out reads with a minimum (-m) and maximum (-M) length\n\t2. Filtered reads are mapped to reference (-r)\n\t3. Read counts are extracted from mapping file and contigs with no reads mapping are filtered out\n\e[0m\nArguments:"

Help()
{
echo -e "This script will perform the following:\n\t\t\t(1) Filter reads (ONT for the original purpose of this script) base on size (-m, --min/-M , --max)\n\t\t\t(2) Perform alignment of reads to reference using minimap2\n\t\t\t(3) Use samtools to convert sam to bam and output tsv files of the number of reads mapped per contig and the depths."   
echo -e "-i, --input\tfasta file of reads for mapping \e[31m[Required]\e[0m"
echo -e "-p, --prefix\tprefix for read file \e[31m[Optional]\e[0m"
echo -e "-r, --reference\treference file for mapping \e[31m[Required]\e[0m"
echo -e "-m, --min\tMinimum read length for mapping [default: -m 200]"
echo -e "-M, --max\tMaximum read length for mapping \e[31m[Optional]\e[0m\n"
echo -e "-t, --threads\tNumber of threads \e[31m[Default: -t 4]\e[0m\n"
}

##########################################################################################################################################
##########################################################################################################################################

while getopts i:p:r:m:M:t:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
        p)prefix=${OPTARG};;
        r)reference=${OPTARG};;
        m)min=${OPTARG};;
        M)max=${OPTARG};;
        t)threads=${OPTARG};;
    h)Help; exit;;
    esac
done

##########################################################################################################################################
##########################################################################################################################################

if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${prefix}" ]]; then 
    prefix=$(basename ${input} .gz | sed 's/\.fasta\?//g; s/\.fastq\?//g; s/\.fa\?//g; s/\.fq\?//g'); fi
if [[ -z "${reference}" ]]; then echo "-r, --reference REQUIRED"; Help; exit; fi
if [[ -z "${min}" ]]; then min=200; fi
if [[ -z "${threads}" ]]; then threads=4; fi

##########################################################################################################################################
##########################################################################################################################################

# Perform the mapping
if [[ -z "${max}" ]]; then
    seqkit seq --min-len ${min} ${input} | minimap2 -ax splice -uf -t ${threads} --secondary=no ${reference} - | samtools sort -O BAM -  > ${prefix}.bam #|samtools view -h -F 0x900
    samtools index -@ ${threads} -b ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
	samtools depth ${prefix}.bam > ${prefix}.depth

    totalReads=$(seqkit seq --min-len ${min} ${input} | seqkit stats -T - | awk '{print $4}' | sed '1d')

elif [[ ! -z "${max}" ]]; then
    seqkit seq --min-len $min --max-len $max ${input} | minimap2 -ax splice -uf -t ${threads} --secondary=no ${reference} - | samtools sort -O BAM - > ${prefix}.bam #|samtools view -h -F 0x900
    samtools index -@ ${threads} -b ${prefix}.bam
    samtools idxstats ${prefix}.bam > ${prefix}.idxstats
	samtools depth ${prefix}.bam > ${prefix}.depth

    totalReads=$(seqkit seq --min-len ${min} --max-len ${max} ${input} | seqkit stats -T - | awk '{print $4}' | sed '1d')

fi

##########################################################################################################################################
##########################################################################################################################################

    # Modify the idxstats file to contains the prefix and remove suffix
awk -v prefix="$prefix" '{print prefix "\t" $0}' ${prefix}.idxstats | sed 's/\.fasta\?//g; s/\.fastq\?//g; s/\.fa\?//g; s/\.fq\?//g' > ${prefix}.idxstats.tmp1
awk '{print $1,$2,$3,$4}' ${prefix}.idxstats.tmp1 > ${prefix}.idxstats.tmp2
	
	# Filter out mapping with 0 reads mapped to contigs
		awk '{if ($4 > 0) print}' ${prefix}.idxstats.tmp3 > ${prefix}.idxstats.tmp4
		mv ${prefix}.idxstats.tmp4 ${prefix}.idxstats
awk -v totalReads="$totalReads" '{print $0 "\t" totalReads}' ${prefix}.idxstats.tmp2 > ${prefix}.idxstats.tmp3

# Modify the depth file to contains the prefix and remove suffix
awk -v prefix="$prefix" '{print prefix "\t" $0}' ${prefix}.depth | sed 's/\.fasta\?//g; s/\.fastq\?//g; s/\.fa\?//g; s/\.fq\?//g' > ${prefix}.depth.tmp1
mv ${prefix}.depth.tmp1 ${prefix}.depth

awk '{if ($4 > 0) print}' ${prefix}.idxstats.tmp3 > ${prefix}.idxstats.tmp4
mv ${prefix}.idxstats.tmp4 ${prefix}.idxstats

# Clean up intermediate files
rm ${prefix}.idxstats.tmp*
rm ${prefix}.depth.tmp*

sed -i 's/ /\t/g' ${prefix}.idxstats
sed -i 's/ /\t/g' ${prefix}.depth
