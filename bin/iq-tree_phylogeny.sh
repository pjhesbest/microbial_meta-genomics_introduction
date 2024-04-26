#!/bin/bash

# Requirements:
# IQ-tree http://www.iqtree.org/
# MAFFT https://mafft.cbrc.jp/alignment/software/
# trimAL http://trimal.cgenomics.org/

##################################################################################################
##############################################################################################

echo -e "\e[34mCOMMAND LINE JOB SUBMISSSION:\n\tbash iq-tree_phylogeny.sh $@\e[0m\n\nArguments:"

Help()  
{
echo -e "This script is to perform the following:\n\t\t\t(1) Multiple sequence alignment with MAFFT\n\t\t\t(2) Trim alignment using TrimAL\n\t\t\t(3) Phylogenetic reconstruction with IQ-TREE\n"
echo -e "-i, --input\t/full/path/to/file.fasta for MSA and phylogeny \e[31m[Required]\e[0m"
echo -e "-b, --bootstraps\tNumber of boottraps for tree \e[31m[Default: -b 200]\e[0m"
echo -e "-R, --Retree\tRerun only the tree using existing alignment \e[31m[Optional]\e[0m"
echo -e "-t, --threads\tNumber of threads \e[31m[Default: -t 4]\e[0m"
}

##################################################################################################
##################################################################################################

while getopts i:b:R:t:h option
do 
    case "${option}" in
        i)input=${OPTARG};;
	b)bootstraps=${OPTARG};;
	R)retree=true;;
  	t)threads=${OPTARG};;
    h)Help; exit;;
    esac

done

##################################################################################################
##################################################################################################

	# check flags are used or set defaults
if [[ -z "${input}" ]]; then echo "-i, --input REQUIRED"; Help; exit; fi
if [[ -z "${bootstraps}" ]]; then bootstraps="200"; fi
if [[ -z "${threads}" ]]; then threads="4"; fi

	# create prefix variable for naming outputs
prefix=$(basename ${input} .gz | sed 's/\.fasta//g; s/\.fa//g; s/\.fq//g; s/\.fastq//g; s/\.fna//g')

##################################################################################################
##################################################################################################

	# Check if the input file has a .gz extension and unzip it if necessary
if [[ "${input}" == *.gz ]]; then
  gunzip -c "${input}" >"${prefix}.fasta"
  input="${prefix}.fasta"
  echo "Unzipped file: ${input}"
fi

##################################################################################################
##################################################################################################

	# Check if the -R flag is used to determine whether to run mafft and iqtree
if [[ -n "${retree}" ]]; then
  
  sed -i 's/>_R_/>/g' ${prefix}.aln.trim.fasta
  iqtree -s "${prefix}.aln.fasta" -m GTR -T AUTO -ntmax "${threads}" -b "${bootstraps}" --prefix "${prefix}" -redo

else
	
  mafft --thread "${SLURM_CPUS_PER_TASK}" --adjustdirectionaccurately --auto "${input}" > "${prefix}.aln.fasta"
  trimal -in "${prefix}.aln.fasta" -out "${prefix}.aln.trim.fasta" -gt 0.9 -cons 60
  sed -i 's/>_R_/>/g' ${prefix}.aln.trim.fasta
  iqtree -s "${prefix}.aln.trim.fasta" -m GTR -T AUTO -ntmax "${threads}" -b "${bootstraps}" --prefix "${prefix}" -redo

fi
