#!/usr/bin/env python3
import questionary
import os
import sys

config = {}

output_prefix = questionary.text("Specify your config file name (e.g. HGjob)").ask()
config['Main'] = questionary.text("Enter the path of SpliceDecoder (e.g. /User/usr/Tool/Splice-decoder-main/)").ask()+"/"
config['code'] = os.path.join(config['Main'], "code/")
config['cpatdb'] = os.path.join(config['Main'], "dat/")
config['input'] = questionary.text("Enter your working directory (e.g. /User/usr/Tool/Splice-decoder-main/project1)").ask()+"/"
config['Your_rMATS'] = questionary.text("Enter your rMATS output path (e.g. /User/usr/Tool/Splice-decoder-main/toy_data)").ask()+"/"
config['target_gene'] = questionary.text("Enter your target gene list (e.g. /User/usr/Tool/Splice-decoder-main/target_genes.tsv)").ask()
config['Your_GTF'] = questionary.text("Enter your GTF file that you used in rMATS with its full path\n(e.g. /User/usr/Tool/Splice-decoder-main/toy_data/toy.gtf or /User/usr/Tool/Splice-decoder-main/toy_data/*.gtf))").ask()

## Scoring
if sys.argv[2] != "toy":
	config['get_score'] = questionary.text("Do you want to calculated the effect score? [yes/no]").ask()
	if config['get_score'] == "no" or config['get_score'] == "No":
		config['Your_TPM'] = "N"
		config['bam_list'] = "N"
		config['tpm'] = "N"
	else:
		config['Your_TPM'] = questionary.text("Enter your TPM matrix with full path (e.g. /User/usr/Tool/Splice-decoder-main/toy_data/tpm.tsv or N)").ask()
		if config['Your_TPM'] == "N":	# SpliceDecoder will calculate TPM by using stringtie with the given bamfiles
			config['bam_list'] = questionary.text("Enter your bamlist which should contains bamfile with their full path in each line\n(e.g. /User/usr/Tool/Splice-decoder-main/toy_data/bam_list.txt or N)").ask()
			config['tpm'] = "Y"
		else:	# SpliceDecoder will use the given TPM matrix
			config['bam_list'] = "N"
			config['tpm'] = "Y_own"
	config['species'] = questionary.text("Enter a species of your data (e.g. human or mouse)").ask()
	config['seq_type'] = questionary.text("Enter a type of GTF (e.g., SR (GENCODE GTF) or LR (Custom GTF) )").ask()


else:	# For toy data set
	config['get_score'] = "no"
	config['Your_TPM'] = "N"
	config['bam_list'] = "N"
	config['tpm'] = "N"
	config["species"] = "mouse"
	config['seq_type'] = "LR"

config['nmd_met'] = questionary.text("Specify a NMD definition method (e.g., default (55rule) or sensitive)").ask()
fdr_input = questionary.text("Enter a FDR cut off for your rMATS (float [0-1], default 0.05)").ask()
config['FDR'] = float(fdr_input) if fdr_input else 0.05
psi_input = questionary.text("Enter a |dPSI| cut off for your rMATS (float [0-1], default 0.1)").ask()
config['PSI'] = float(psi_input) if psi_input else 0.1
config['njobs'] = questionary.text("Enter a number of cpu in splice-decoder job (int [5-?], default 5").ask()

## Specify 3rd party tool path
config['conda'] = sys.argv[1]
config['cpat'] = os.path.join(config['conda'], "bin/cpat")
config['bedtools'] = os.path.join(config['conda'], "bin/bedtools")
with open(f"{config['Main']}{output_prefix}.config", "w") as f:
    for key, value in config.items():
        f.write(f'{key}="{value}"\n')

print(f"### Your {output_prefix}.config is available here: {config['Main']} \nYou can use `bash Main.sh $Function {output_prefix}.config` or `sbatch Main.sh $Function {output_prefix}.config` in your terminal.")
print(f"### You can get more information for the Function variable of Main.sh by using bash `Main.sh -h`")
