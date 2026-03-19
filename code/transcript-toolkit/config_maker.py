#!/usr/bin/env python3
import questionary
import os
import sys

config = {}

output_prefix = questionary.text("Specify your config file name (e.g. HGjob)").ask()
config['Main'] = questionary.text("Enter the path of Splice-decoder (e.g. /User/usr/Tool/Splice-decoder-main/)").ask()+"/"
config['code'] = os.path.join(config['Main'], "code/transcript-toolkit/")
config['cpatdb'] = os.path.join(config['Main'], "dat/")
config['input'] = questionary.text("Enter your working directory (e.g. /User/usr/Tool/Splice-decoder-main/project1)").ask()+"/"
config['canonical'] = questionary.text("Specify a strategy to define canonical transcripts (it can be 'default' or 'user'").ask()
config['gene_list'] = questionary.text("Enter your gene_list. If you answered 'user' for previous question,please put your gene_list, it must contain gene and canonical transcript IDs (tab separated)\n(e.g. /User/usr/Tool/Splice-decoder-main/toy_data/gene_list.txt))").ask()
config['Your_GTF'] = questionary.text("Enter your query GTF file with its full path\n(e.g. /User/usr/Tool/TX-comp-main/toy_data/toy.gtf or /User/usr/Tool/Splice-decoder-main/toy_data/*.gtf))").ask()
config['species'] = questionary.text("Enter a species of your data (e.g. human or mouse)").ask()
config['seq_type'] = questionary.text("Enter a sequencing type of your data (e.g., SR (short-read) or LR (long-read) )").ask()
config['njobs'] = questionary.text("Enter a number of cpu in TX-comp job (int [0-?]").ask()

## Specify 3rd party tool path
config['conda'] = sys.argv[1]
config['cpat'] = os.path.join(config['conda'], "bin/cpat")
config['bedtools'] = os.path.join(config['conda'], "bin/bedtools")
config['query_list'] = os.path.join(config['input'], "query_list.txt")
with open(f"{config['Main']}{output_prefix}.config", "w") as f:
    for key, value in config.items():
        f.write(f'{key}="{value}"\n')

print(f"{config['Main']}{output_prefix}.config has been created! \nYou can use `bash Main.sh $Function {output_prefix}.config` or `sbatch Main.sh $Function {output_prefix}.config` in your terminal.")
print(f"You can get more information for the Function variable of Main.sh by using bash `Main.sh -h`")
