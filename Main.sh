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

input_var=$1
case ${input_var} in
    help|-h|--help)
	echo "#################################################"
        echo "You should specify your job: Make_input / DS_mapping / ORF_mapping / Simulation / Scoring / all AND \${Your_config}"
	echo
        echo "[Usage example]"
	echo "bash $0 Make_input \${Your_config}"
        echo "#################################################"
        exit 1
        ;;
esac

config=$2
source ${config} 

## Prepare working dir
mkdir -p ${input}
mkdir -p ${input}/tpm/

## Prepare log
timestamp=$(date +%Y%m%d_%H%M%S)
logfile="${input}/SD_${timestamp}.log"
exec > >(tee -a "$logfile") 2>&1

## Check rmat.csv
if [ -f "${input}/rmat.csv" ]; then
        echo "rmats file already exists"
else
        # echo "Processed rmats file does not exist"
        echo "Making a processed rmats file..."
        python ${code}/NEW_make_input_from_rmats.py -r ${Your_rMATS} -i ${input} -fdr ${FDR} -psi ${PSI}
fi

cd ${input}

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

## TPM file setting
if [[ ${tpm} == "Y" ]]; then 
	bash ${code}/stringtie.sh ./main.gtf ${bam_list} ${input}/tpm/
	tab_files=`ls ${input}/tpm | grep .tab  | grep -v ctab`	# Udated 25.01.19
	for tab in ${tab_files}
	do
		python ${code}/tpm_process.py -i ${input}/tpm/ -n ${tab}
	done

elif [[ ${tpm} == "Y_own" ]]; then
	cd ${input}/tpm/
	if [[ ! -L "./tpm_matrix.tsv" ]]; then
	    cat ${Your_TPM} | sed 's/,/\t/g' > ./tpm_matrix.tsv # Prevent to formatting issue
        # ln -s ${Your_TPM} ./tpm_matrix.tsv
	fi 
fi

## SpliceDecoder Main
case ${input_var} in
    Make_input)
	echo "Running input processing"

	## CLEAN UP pre-existed gtf
	ls ${input}*gtf | grep exon | while read file
	do
	rm ${file}
	done

	## Make exon gtf from the given gtf
	cat ${input}*gtf | awk -F "\t" '{if ($3 == "exon") print }' > ${input}exon_only.gtf

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
			python ${code}/00-1_add_exon_n.py -i ${input} -s ${seq_type}
		fi
	fi

	## Final input processing
	python ${code}/00-2_processing_gtf.py -e ${input}exon_only.gtf -o ${input}
	python ${code}/00-3_make_dict.py -i ${input}
 
	echo "Finished input processing"
	;;

    DS_mapping)
	echo "Running the DS Mapping"
	## Mapping
	python -W ignore ${code}/01-1_exon_coordinate_v5.py -i ${input} -t ${target_gene}
	python -W ignore ${code}/01-2_stat_exon_coor_result.py -i ${input} -t ${target_gene}
	
	echo "Finished the DS Mapping"
	;;

    ORF_mapping)
	echo "Running ORF mapping"
	## ORF mapping
	rm -r ${input}/sim_bed
	rm -r ${input}/get_Fasta
	rm -r ${input}/cpat
    mkdir -p ${input}/sim_bed
    mkdir -p ${input}/get_Fasta
    mkdir -p ${input}/cpat
	python ${code}02-1_ORF_mapping.py -i ${input} -t ${njobs} -cp ${cpat} -cpdb ${cpatdb} -b ${bedtools} -f ${genomefa} -p ${species}
	echo "Finished ORF mapping"
	;;

    Simulation)
	echo "Running Splicing simulation"
	## Splicing simulation
	mkdir -p ${input}/result
	ls ${input}/sim_bed/ | grep merged | cut -f 2 -d "_" | sed 's/[0-9]\+$//g' | sort -u | \
	parallel --jobs 5 '
	    echo "Start {} Simulation"
	    python -W ignore '${code}'02-2_SimTX.py -s {} -i '${input}' -t "all" -p '${Main}${config}' -c '${code}' -n '${nmd_met}'
	'
	echo "Finished Splicing simulation"
	;;

    Scoring)
	if [ ${get_score} == "Yes" ] || [ ${get_score} == "yes" ]
	then
	    echo "Calculated effect score"
		python ${code}03-2_overview_change_rate.py -i ${input}
		python ${code}03-3_scoring_function.py -i ${input} -t ${tpm}
		python ${code}Make_summary.py -i ${input}
	    echo "Finished"
	    if [ -f "${input}/result/Effect_score.tsv" ] && [ $(wc -l < ${input}/result/Effect_score.tsv) -ge 5 ]; then
				bash ${code}Clean_up.sh ${input}
        fi
    else
        bash ${code}Clean_up.sh ${input}
	fi
	;;

    all)
        echo "Running input processing"

        ## CLEAN UP pre-existed gtf
        ls ${input}*gtf | grep exon | while read file
        do
        rm ${file}
        done

        ## Make exon gtf from the given gtf
        cat ${input}*gtf | awk -F "\t" '{if ($3 == "exon") print }' > ${input}exon_only.gtf

        ## Check ths existence of "exon_number" tag
        exon_n=`less ${input}exon_only.gtf | grep -v exon_number | wc -l`

        ## Check exon number tag
	## Long-read RNAseq
        if [ ${seq_type} == "LR" ]
        then
		cat ${input}exon_only.gtf | grep -v exon_number > ${input}LR_exon_only.gtf
                cat ${input}exon_only.gtf | grep exon_number > ${input}non_pacbio_exon_only.gtf

                if [ ${exon_n} != 0 ]   # If the exon gtf has exon number
                then
                        python ${code}/00-1_add_exon_n.py -i ${input} -s ${seq_type}
                else    # If the exon gtf does not have exon number
                        cat ${input}non_pacbio_exon_only.gtf ${input}LR_exon_only.gtf > ${input}exon_only.gtf
                fi
        fi
	
	## Short-read RNAseq
        if [ ${seq_type} != "LR" ]
        then
                if [ ${exon_n} != 0 ]
                then
                        python ${code}/00-1_add_exon_n.py -i ${input} -s ${seq_type}
                fi
        fi

        ## Final input processing
        python ${code}/00-2_processing_gtf.py -e ${input}exon_only.gtf -o ${input}
		python ${code}/00-3_make_dict.py -i ${input}
 
        echo "Finished input processing"

       
        echo "Running the DS Mapping"
        ## Mapping
        python -W ignore ${code}/01-1_exon_coordinate_v5.py -i ${input} -t ${target_gene}
        python -W ignore ${code}/01-2_stat_exon_coor_result.py -i ${input} -t ${target_gene}
        echo "Finished the DS Mapping"


        echo "Running ORF mapping"
        ## ORF mapping
		rm -r ${input}/sim_bed
        rm -r ${input}/get_Fasta
        rm -r ${input}/cpat
        mkdir -p ${input}/sim_bed
        mkdir -p ${input}/get_Fasta
        mkdir -p ${input}/cpat
        python ${code}02-1_ORF_mapping.py -i ${input} -t ${njobs} -cp ${cpat} -cpdb ${cpatdb} -b ${bedtools} -f ${genomefa} -p ${species}
        echo "Finished ORF mapping"

 
        echo "Running Splicing simulation"
        ## Splicing simulation
		mkdir -p ${input}/result
        ls ${input}/sim_bed/ | grep merged | cut -f 2 -d "_" | sed 's/[0-9]\+$//g' | sort -u | \
        parallel --jobs 5 '
            echo "Start {} Simulation"
            python -W ignore '${code}'02-2_SimTX.py -s {} -i '${input}' -t "all" -p '${Main}${config}' -c '${code}' -n '${nmd_met}'
        '
        echo "Finished Splicing simulation"

        if [ ${get_score} == "Yes" ] || [ ${get_score} == "yes" ]
        then
	        echo "Calculated effect score"
	        python ${code}03-2_overview_change_rate.py -i ${input}
	        python ${code}03-3_scoring_function.py -i ${input} -t ${tpm}
	        python ${code}Make_summary.py -i ${input}
	        echo "Finished"

            if [ -f "${input}/result/Effect_score.tsv" ] && [ $(wc -l < ${input}/result/Effect_score.tsv) -ge 5 ]; then
                    bash ${code}Clean_up.sh ${input}
            fi
        else
            bash ${code}Clean_up.sh ${input}
		fi
	;;

    *)
	echo "##################################################################################################################"
	echo "[ERROR]"
        echo "You should specify your job: Make_input / DS_mapping / ORF_mapping / Simulation / Scoring / all AND \${Your_config}"
        echo
        echo "[Usage example]"
        echo "bash $0 Make_input \${Your_config}"
        echo "##################################################################################################################"
        exit 1
	;;
esac
