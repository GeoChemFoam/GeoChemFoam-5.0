#!/bin/bash
#SBATCH --job-name=initCaseFlow
#SBATCH --time=1:00:0
#SBATCH --ntasks=1
#SBATCH --account=ecseaj02
#SBATCH --partition=serial
#SBATCH --qos=serial

# Define memory required for this job. By default, you would
# get just under 2 GB, but you can ask for up to 125 GB.
#SBATCH --mem=4G

# Set the number of threads to 1 to avoid auto-threading
export OMP_NUM_THREADS=1

# Configure GCF
module load openfoam/com/v2212
source /work/ecseaj02/ecseaj02/gavingcf/works/GeoChemFoam-5.0/etc/bashrc

# Configure Python
module load cray-python
source /work/ecseaj02/ecseaj02/gavingcf/myvenv/bin/activate

# Ensure we are in the correct working directory
cd /work/ecseaj02/ecseaj02/gavingcf/works/GeoChemFoam-5.0/runs/Ketton

# initCaseFlow needs time=1:00:0
./initCaseFlow.sh  

sbatch workflow_step_4.bat

