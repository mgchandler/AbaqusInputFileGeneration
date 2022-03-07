#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=4000mb
#PBS -l walltime=05:59:59

cd /home/mc16535/Abaqus_Sidewalls/Initial_Testing

# Load modules
module load apps/abaqus/2018

# Get useful information

echo JOB ID: ${PBS_JOBID}
echo ARRAY ID: ${PBS_ARRAY_INDEX}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname
Job=${PBS_ARRAY_INDEX}

echo Truss Analysis Test

# Execute code
echo abaqus job=TrussAnalysis inp=TrussAnalysisJob.inp
abaqus job=TrussAnalysis inp=TrussAnalysisJob.inp

echo End Time: $(date)