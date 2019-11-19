#!/bin/bash
#SBATCH --job-name=forestfire
#SBATCH --mail-type=END
#SBATCH --mail-user=ssh.sentian@gmail.com
#SBATCH --output=/home/st1864/boss/output/orthx/%A_%a.out # master job id %A and array-task job id %a
#SBATCH --array=1-10
#SBATCH --time=6:00:00

#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB


module purge
module load gurobi/8.1.1
module load r/intel/3.6.0

echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
cd /home/st1864/boss/code/run_model

Rscript run_forestfire.R $SLURM_ARRAY_TASK_ID