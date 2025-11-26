#!/bin/bash

yml=pkg.yml
source "$(conda info --base)/etc/profile.d/conda.sh"
conda config --add channels conda-forge
conda config --add channels bioconda
conda install -n base -c conda-forge mamba
mamba env create --solver libmamba --name splice-decoder -f ${yml}
conda activate splice-decoder
unzip toy_data.zip

## Prepare data source
unzip dat.zip
cd dat/
cd reference/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz
gunzip -c GRCh38.primary_assembly.genome.fa.gz > hg38.fa
samtools faidx hg38.fa
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
gunzip -c GRCm38.primary_assembly.genome.fa.gz > GRCm38.fa
samtools faidx GRCm38.fa
