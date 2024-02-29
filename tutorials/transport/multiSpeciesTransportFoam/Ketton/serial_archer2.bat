#!/bin/bash
#SBATCH --job-name=serial_jobs
#SBATCH --time=0:20:0
#SBATCH --ntasks=1
#SBATCH --account=ecseaj02
#SBATCH --partition=serial
#SBATCH --qos=serial

# Define memory required for this job. By default, you would
# get just under 2 GB, but you can ask for up to 125 GB.
#SBATCH --mem=4G

# Set the number of threads to 1 to avoid auto-threading
export OMP_NUM_THREADS=1

# Configure GCF using Archer2 module
module load gcfoam/5.0
source $GCFOAM_DIR/etc/bashrc_archer2

# Configure Python
module load cray-python
source /work/ecseaj02/ecseaj02/gavingcf/myvenv/bin/activate

# Choose just one of the following serial scripts
./createMesh.sh  
#./initCaseFlow.sh  
#./initCaseTransport.sh 

