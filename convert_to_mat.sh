#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=1000mb
#PBS -l walltime=00:10:00

# Load modules
module load apps/matlab/r2019a

# Change to working directory. This should be the same location as the control Matlab script.
cd /work/mc16535/Abaqus_Sidewalls/L_Vis

# Execute code

matlab -nodisplay -nodesktop -singleCompThread -nosplash -r "FMC_read('/work/mc16535/Abaqus_Sidewalls/L_Vis', 'L_Vis');quit;"