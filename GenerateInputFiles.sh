#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=1000mb
#PBS -l walltime=00:04:59

# Load modules
module load lang/python/anaconda/3.8-2020.07

# Make a directory to store .inp files.
export JOBNAME="L_Obs"

cd /work/mc16535/Abaqus_Sidewalls/
mkdir ${JOBNAME}

# Change to the submission directory, and copy everything to the final dir location.
cd ${PBS_O_WORKDIR}
cp * /work/mc16535/Abaqus_Sidewalls/${JOBNAME}

cd /work/mc16535/Abaqus_Sidewalls/${JOBNAME}

# Get useful information

echo JOB ID: ${PBS_JOBID}
echo ARRAY ID: ${PBS_ARRAY_INDEX}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname

echo ${JOBNAME} .inp generation from python

# Execute code
echo python AbaqusInputGeneration_SRM.py ${JOBNAME}
python AbaqusInputGeneration_SRM.py ${JOBNAME}

echo End Time: $(date)