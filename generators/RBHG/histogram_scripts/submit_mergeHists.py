import os
import subprocess

# List of directories

directories = """
/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/CPP_PS1/real_FFN/CPP_2205_ee_Berlin_v110/135DEG
/work/halld/home/acschick/channels/epemmissprot/newRBHG/output/CPP_PS1/real_FFN/CPP_2205_ee_Berlin_v110/45DEG
"""
directories = directories.strip().splitlines()


# Path to the merge_histograms.py script
merge_script_path = '/work/halld/home/acschick/channels/epemmissprot/newRBHG/mergeHists/mergeHists.py'

# Template for the sbatch script
sbatch_script_template = """#!/bin/bash
#SBATCH --job-name=mergeHists
#SBATCH --output=mergeHists_%j.out
#SBATCH --error=mergeHists_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=500mb
#SBATCH --cpus-per-task=1
#SBATCH --partition=production

# Source the gxenv command
source /group/halld/Software/build_scripts/gluex_env_boot_jlab.sh

# Set up the environment that has ROOT built with python
gxenv 

# Execute the Python script
python {} {}
"""

for directory in directories:
    # Create sbatch script for each directory
    sbatch_script = sbatch_script_template.format(merge_script_path, directory)

    # Write the sbatch script to a temporary file
    with open('temp_sbatch_script.sh', 'w') as file:
        file.write(sbatch_script)

    # Submit the job
    try:
        subprocess.run(['sbatch', 'temp_sbatch_script.sh'], check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while submitting the job for directory {directory}: {e}")

    # Optional: Remove the temporary sbatch script file
    os.remove('temp_sbatch_script.sh')
