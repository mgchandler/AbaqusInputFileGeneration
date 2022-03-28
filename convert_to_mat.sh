#!/bin/bash
#
#SBATCH --job-name=I_45npw_02mm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5500mb
#SBATCH --time=00:14:59

JobName="${SLURM_JOB_NAME}_1"

# Load modules
module load apps/matlab/r2019a
module load apps/abaqus/2018
module load lang/python/anaconda/3.8-2020.07

# Change to working directory. This should be the same location as the control Matlab script.
cd /user/work/mc16535/Abaqus_Sidewalls/${SLURM_JOB_NAME}/${JobName}

# Execute code

# EXE1="abaqus python write_single_frame.py ${JobName}_16.odb 160"
# echo ${EXE1}
# ${EXE1}

matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "FMC_read('/user/work/mc16535/Abaqus_Sidewalls/${SLURM_JOB_NAME}/${JobName}', '${JobName}_BP');quit;"
