#!/bin/bash
#SBATCH --job-name=parallel
#SBATCH --time=0:20:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=ecseaj02
#SBATCH --partition=standard
#SBATCH --qos=standard

# Set the number of threads to 1 to avoid auto-threading
export OMP_NUM_THREADS=1

# Propagate the cpus-per-task setting from script to srun commands
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Configure GCF using Archer2 module
module load gcfoam/5.0
source $GCFOAM_DIR/etc/bashrc_archer2

# Configure Python
module load cray-python
source /work/ecseaj02/ecseaj02/gavingcf/myvenv/bin/activate

# Exit if the total number of tasks do not equal the total number of processor directories.
export NP="$(find processor* -maxdepth 0 -type d -print| wc -l)"
if [ "$NP" -ne "$((SLURM_TASKS_PER_NODE*SLURM_JOB_NUM_NODES))" ]; then
     echo "Error: number of tasks does not equal the number of processor directories"
     echo "SBATCH ntasks-per-node times nodes must equal $NP"
     echo "However they are $SLURM_TASKS_PER_NODE and $SLURM_JOB_NUM_NODES, respectively"
     exit 1
fi
echo -e "Run job_name in parallel on $NP processors"

# Choose just one of the following parallel scripts.
./runSnappyHexMesh.sh 
#./runCaseFlow.sh 
#./processFlow.sh 
#./runCaseTransport.sh 
#./processTransport.sh 


