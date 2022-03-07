#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=4000mb
#PBS -l walltime=05:59:59

# Load modules
module load apps/abaqus/2018

# Change to working directory. This should be the same location as the control Matlab script.
cd /work/mc16535/Abaqus_Sidewalls/L_Vis

# Get useful information

echo JOB ID: ${PBS_JOBID}
echo ARRAY ID: ${PBS_ARRAY_INDEX}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname
Job=${PBS_ARRAY_INDEX}

echo Rectangle Geometry

# Execute code
echo abaqus job=FMC_32els_L_SDH_vis_1 inp=FMC_32els_L_SDH_vis_1.inp
abaqus job=FMC_32els_L_SDH_vis_1 inp=FMC_32els_L_SDH_vis_1.inp

echo End Time: $(date)