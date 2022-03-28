#!/bin/bash
#
#SBATCH --job-name=I_45npw_02mm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8000mb
#SBATCH --time=01:59:59

# Load modules
module load lang/python/anaconda/3.8-2020.07

# Get useful information

echo JOB ID: ${SLURM_JOB_ID}
echo ARRAY ID: ${SLURM_ARRAY_TASK_ID}
echo Working Directory: $(pwd)
echo Start Time: $(date)
echo Run location:
hostname

# Make a directory to store .inp files.
cd /user/work/mc16535/Abaqus_Sidewalls/
mkdir ${SLURM_JOB_NAME}

# Change to the submission directory, and copy everything to the final dir location.
cd ${SLURM_SUBMIT_DIR}
cp * /user/work/mc16535/Abaqus_Sidewalls/${SLURM_JOB_NAME}

# Execute code
cd /user/work/mc16535/Abaqus_Sidewalls/${SLURM_JOB_NAME}
echo python AbaqusInputGeneration.py ${SLURM_JOB_NAME}
python AbaqusInputGeneration.py ${SLURM_JOB_NAME}

# Copy everything into subdirectories
for Job in $(seq 1 5); do
	if [ $Job -ne 5 ]; then
		JobIdx="${Job}"
	else
		JobIdx="b"
	fi
	
	mkdir ${SLURM_JOB_NAME}_${JobIdx}
	mv ${SLURM_JOB_NAME}_${JobIdx}* ./${SLURM_JOB_NAME}_${JobIdx}
	cp ${SLURM_JOB_NAME}.yml ./${SLURM_JOB_NAME}_${JobIdx}
	cp compile_timetraces.py ./${SLURM_JOB_NAME}_${JobIdx}
	cp write_single_timetrace.py ./${SLURM_JOB_NAME}_${JobIdx}

	
done

echo End Time: $(date)
