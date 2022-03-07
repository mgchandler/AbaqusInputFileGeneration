#!/bin/bash
#
#
#PBS -l select=1:ncpus=1:mem=8250mb
#PBS -l walltime=09:59:59
#PBS -J 1-32

export ARRAYJOBNAME="L_Obs"

cd /work/mc16535/Abaqus_Sidewalls/${ARRAYJOBNAME}

# Load modules
module load apps/abaqus/2018
module load lang/python/anaconda/3.8-2020.07

# Get useful information
echo JOB ID: ${PBS_JOBID}
echo ARRAY ID: ${PBS_ARRAY_INDEX}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname
Job=${PBS_ARRAY_INDEX}

echo ${ARRAYJOBNAME} Job no ${Job}

# Define executable so the function is written only once
export THISJOBNAME="${ARRAYJOBNAME}_${Job}"
mkdir ${THISJOBNAME}
cp ${THISJOBNAME}.inp ./${THISJOBNAME}
cd ./${THISJOBNAME}

export EXE1="abaqus job=${THISJOBNAME} input=${THISJOBNAME}.inp mem=8000MB scratch=. interactive"
echo ${EXE1}

# Execute
${EXE1}

echo Abaqus completed analysis. Starting .dat extraction.

cp ${THISJOBNAME}.odb /work/mc16535/Abaqus_Sidewalls/${ARRAYJOBNAME}
cd /work/mc16535/Abaqus_Sidewalls/${ARRAYJOBNAME}
rm -f ${THISJOBNAME}/*
rm -d ${THISJOBNAME}

export EXE2="abaqus python write_single_timetrace.py ${THISJOBNAME}.odb"
echo ${EXE2}
${EXE2}

# rm -f ${THISJOBNAME}.odb

echo Timetraces extracted. Compiling into FMC data.
export EXE3="python compile_timetraces.py ${THISJOBNAME}"
echo ${EXE3}
${EXE3}

rm -f ${THISJOBNAME}.dat

echo End Time: $(date)