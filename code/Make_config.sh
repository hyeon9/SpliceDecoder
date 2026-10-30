#!/bin/bash

echo "To make right config file, you should give right answer for each question!"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate splice-decoder
conda_path="$CONDA_PREFIX"
export PATH="$(pwd):$PATH"
toy_key="${1:-data}"

config_maker.py "$conda_path" "$toy_key"
