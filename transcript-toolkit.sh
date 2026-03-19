#!/bin/bash
#SBATCH -J Run_SD
#SBATCH -c 45
#SBATCH --mem=60G
#SBATCH -o ./%u-%x-%j
#SBATCH -e ./%u-%x-%j.err
#SBATCH --time=72:00:00

## Load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate splice-decoder

config=$1
source ${config} 

## Prepare working dir
mkdir -p ${input}
cd ${input}

## Prepare log
timestamp=$(date +%Y%m%d_%H%M%S)
logfile="${input}/SD_${timestamp}.log"
exec > >(tee -a "$logfile") 2>&1

## Genome ver setting
if [[ ${species} == "human" ]]; then
	genomefa=${Main}"/dat/reference/hg38.fa"  # Reference genome version
else
	genomefa=${Main}"/dat/reference/GRCm38.fa"  # Reference genome version
fi

## GTF file setting
if [[ ! -L "./main.gtf" ]]; then
	ln -s ${Your_GTF} ./main.gtf
fi


## SpliceDecoder Main
echo "Running input processing"

## CLEAN UP pre-existed gtf
ls ${input}*gtf | grep exon | while read file
do
rm ${file}
done

## Make subset of exon gtf from the given gtf
python ${code}gtf_filter.py -i ${input} -g ${gene_list} -p ${Main} -c ${canonical}

## Check ths existence of "exon_number" tag
exon_n=`less ${input}exon_only.gtf | grep -v exon_number | wc -l`

## Check exon number tag
if [ ${seq_type} == "LR" ]
then
	cat ${input}exon_only.gtf | grep -v exon_number > ${input}LR_exon_only.gtf
	cat ${input}exon_only.gtf | grep exon_number > ${input}non_pacbio_exon_only.gtf

	if [ ${exon_n} -ne 0 ]	# If the exon gtf has exon number
	then
		echo "add exon_number"
	        python ${code}/00-1_add_exon_n.py -i ${input} -s ${seq_type}
	else	# If the exon gtf does not have exon number
		cat ${input}non_pacbio_exon_only.gtf ${input}LR_exon_only.gtf > ${input}exon_only.gtf
	fi
fi

if [ ${seq_type} != "LR" ]
then
	if [ ${exon_n} != 0 ]
	then
                echo "add exon_number"
		python ${code}/00-1_add_exon_n.py -i ${input} -s ${seq_type}
	fi
fi

## Final input processing
python ${code}/00-2_processing_gtf.py -e ${input}exon_only.gtf -o ${input}
echo "Finished input processing"

## Make an processed input
python ${code}00_Make_comp_input.py -i ${input} -t ${njobs} -cp ${cpat} -cpdb ${cpatdb} -b ${bedtools} -f ${genomefa} -p ${species} -q ${query_list}

## Make a result table
python ${code}/01_comp_domain.py -i ${input} -q ${query_list} -p ${Main}${config} -c ${code}
