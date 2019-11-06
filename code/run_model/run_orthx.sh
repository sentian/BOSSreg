#!/bin/bash
#SBATCH --job-name=orthx
#SBATCH --mail-type=END
#SBATCH --mail-user=ssh.sentian@gmail.com
#SBATCH --output=/home/st1864/boss/output/orthx/%A_%a.out # master job id %A and array-task job id %a
#SBATCH --array=1-72
#SBATCH --time=48:00:00

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB


module purge
module load r/intel/3.6.0

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
cd /home/st1864/boss/code

Rscript run.R 'TRUE' $SLURM_ARRAY_TASK_ID