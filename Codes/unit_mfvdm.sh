#!/bin/sh
#SBATCH -n 20            # cores
#SBATCH --mem 100GB     # memory
#SBATCH -t 02-00:00:00   # duration
#SBATCH -p compute      # partition name
#SBATCH -J sa_mfvdm    # job name
#SBATCH --output=job_output.log    # Output log file
#SBATCH --error=job_error.log      # Error log file

echo "Starting SLURM job at $(date)"
echo "Loading Python..."

module load python/3.11.6-gcc-13.1.0-rri7oiq
module load py-pip/23.1.2-gcc-13.1.0-k6tgxbx

echo "The virtual environment is being loaded"
source /home/users/brra/HPC/sa/bin/activate

echo "Using Python from: $(which python)"
echo "Python version: $(python --version)"
echo "Loaded Python packages:"
pip list

DENSITY=$1

echo "Running the python file"
python unit_mfvdm.py $DENSITY
