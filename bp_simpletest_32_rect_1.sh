#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=1000mb
#PBS -l walltime=00:19:59

# Load modules
module load apps/abaqus/2018

# Change to working directory. This should be the same location as the control Matlab script.
cd ${PBS_O_WORKDIR}

# Get useful information

echo JOB ID: ${PBS_JOBID}
echo ARRAY ID: ${PBS_ARRAY_INDEX}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname

echo Valid Path Testing - Direct SDH

# Execute code
echo abaqus job=TrussAnalysisJob inp=TrussAnalysisJob.inp
abaqus job=TrussAnalysisJob inp=TrussAnalysisJob.inp

echo End Time: $(date)