#!/bin/bash
#
#SBATCH --job-name=I_45npw_02mm_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16250mb
#SBATCH --time=99:59:59
#SBATCH --array=1-8

# Global variables used in all scripts
MaxArrayNo=$SLURM_ARRAY_TASK_COUNT # Read this from SLURM directive above
JobsPerArray=4
MaxNoJobs=$MaxArrayNo*$JobsPerArray # This should equal the total number of elements in the probe.
# Make sure that MaxNoJobs / MaxArrayNo is an integer.
JobName=${SLURM_JOB_NAME}

Memory=$(($SLURM_MEM_PER_NODE - 500))

# Global parameters based on variables
ThisArrayNo=${SLURM_ARRAY_TASK_ID}


# Load modules required
module load apps/abaqus/2018
module load lang/python/anaconda/3.8-2020.07
unset SLURM_GTIDS
set -e

# Change to storage location
cd /user/work/mc16535/Abaqus_Sidewalls/${JobName::-2}/${JobName}

# Print useful information to output file
echo \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
echo \# BP Job ID: ${SLURM_JOB_ID}
echo \# BP Array ID: ${SLURM_ARRAY_TASK_ID}
echo \# Working Directory: $(pwd)
echo \# Start Time: $(date)
echo \# Run location:
hostname
echo \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

# Begin loop over each Abaqus run
for Job in $(seq 1 $JobsPerArray); do
	
	# Get the index of the run we are going to perform next
	ThisRunIdx=$((JobsPerArray * (ThisArrayNo - 1) + Job))
	# Get the name of the input file which will be run
	ThisRunName="${JobName}_${ThisRunIdx}"
	echo \# ${ThisRunName}                  \#
	echo \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
	echo Run Start Time: $(date)
	# We are currently in the container directory for this job. Make a 
	# subdirectory to store Abaqus output files in. Copy input file.
	mkdir ${ThisRunName}
	cp ${ThisRunName}.inp ./${ThisRunName}
	cd ./${ThisRunName}
	
	# Define command as executable so that we only write the function once.
	# Do this so we can debug by writing the cmd to the output file
	echo Starting Abaqus execution.
	EXE1="abaqus job=${ThisRunName} input=${ThisRunName}.inp mem=${Memory}MB scratch=. interactive"
	echo ${EXE1}
	${EXE1}
	echo Abaqus completed analysis.
	
	# Copy output database to parent directory, and delete other output files
	cp ${ThisRunName}.odb /user/work/mc16535/Abaqus_Sidewalls/${JobName::-2}/${JobName}
	cd /user/work/mc16535/Abaqus_Sidewalls//${JobName::-2}/${JobName}
	rm -f ${ThisRunName}/*
	rm -d ${ThisRunName}
	
	# Define executable which will extract timetraces from .odb
	echo Starting Timetrace Extraction from .odb
	EXE2="abaqus python write_single_timetrace.py ${ThisRunName}.odb"
	echo ${EXE2}
	${EXE2}
	echo Timetraces extracted and stored as .dat files
	
	# Delete database if no longer required
	# if [ $ThisRunIdx -ne 16 ]; then
	# 	rm -f ${ThisRunName}.odb
	# fi
	
	# Compile timetraces at each node into timetraces for each element in the probe
	echo Starting Timetrace compilation.
	EXE3="python compile_timetraces.py ${ThisRunName}"
	echo ${EXE3}
	${EXE3}
	echo Timetraces compiled. Run FMC_read.m to assemble FMC into single .mat file.
	rm -f ${ThisRunName}.dat
	
	echo Run End Time: $(date)

done

echo End Time: $(date)