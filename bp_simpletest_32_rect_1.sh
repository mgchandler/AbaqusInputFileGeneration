#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=1000mb
#PBS -l walltime=00:19:59
#PBS -J 1-32

# Load modules
module load apps/abaqus/2018

# Change to working directory. This should be the same location as the control Matlab script.
cd /work/mc16535/Abaqus_Sidewalls/Rectangle_Geom

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
echo abaqus job=FMC_32els_rect_SDH_${Job} inp=FMC_32els_rect_SDH_${Job}.inp
abaqus job=FMC_32els_rect_SDH_${Job} inp=FMC_32els_rect_SDH_${Job}.inp

echo End Time: $(date)