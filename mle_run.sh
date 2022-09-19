#!/bin/bash
# 2022-06-10

# Use "debug" for debugging and "tier3" for the actual run
#SBATCH --partition debug  ################ set to tier3 for actual run!

# Establish various settings
#SBATCH --time 5:00:00:00  # job timeout as HH:MM:SS  ################ set to 23:59:59 for actual run!
#SBATCH --account=grn  # RC account name
#SBATCH --mem=2048  # job memory limit in MB

python3 main.py $1 > $2
