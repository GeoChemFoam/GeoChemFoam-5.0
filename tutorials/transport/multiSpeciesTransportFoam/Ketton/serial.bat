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

# Configure GCF using own installation
#module load openfoam/com/v2212
#source /work/ecseaj02/ecseaj02/gavingcf/works/GeoChemFoam-5.0/etc/bashrc

# Configure GCF using Archer2 module
module load gcfoam/5.0
source $FOAM_INSTALL_DIR/etc/bashrc
source $GCFOAM_DIR/ThirdParty/bashrc

# Configure Python
module load cray-python
source /work/ecseaj02/ecseaj02/gavingcf/myvenv/bin/activate

# createMesh needs time=0:20:0
./createMesh.sh  

# initCaseFlow needs time=1:00:0
#./initCaseFlow.sh  

# initCaseTransport needs time=0:20:0
#./initCaseTransport.sh 

