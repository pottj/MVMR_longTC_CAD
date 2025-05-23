#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J MVMRs

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=0:05:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 76 on icelake; 112 on sapphire
#SBATCH --cpus-per-task=1

## Each task is allocated 3.3G (icelake) or 6.7G (icelake-himem) or 4.6G (sapphire)
## If this is insufficient, uncomment and edit this line.
## Maximum value 256G (icelake/sapphire) or 512G (icelake-himem)
#SBATCH --mem=192G

## The system can send emails when your job starts and stops.
## Values include BEGIN, END, ALL, and TIME_LIMIT_80 and TIME_LIMIT_90 
## (reaching 80% or 90% of time limit.) Specify ARRAY_TASKS to receive
## a separate mail for each task. Multiple values can be given, separated by a comma.
#SBATCH --mail-type=FAIL

## The project account name.
## Use mrc-bsu-sl2-cpu for icelake and mrc-bsu-sl2-gpu for ampere
#SBATCH -A mrc-bsu-sl2-cpu

## The partition. Use icelake for normal jobs, or icelake-himem if needed.
#SBATCH -p icelake

## GPU jobs only:
## Uncomment and specify the number of GPUs required per node, maximum 4.
## Note that there is a maximum of 3 cores per GPU.
## The gpu partition is ampere.
## #SBATCH --gres=gpu:1

## Array jobs:
## Start multiple jobs at once.
## Note that resources (cores, memory, time) requested above are for each
## individual array task, NOT the total array.
## #SBATCH --array=1-4
#SBATCH --array=13-14

##  - - - - - - - - - - - - - -

## Section 2: Modules

# All scripts should include the first three lines.

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl            # REQUIRED - loads the basic environment

# Load the latest R version.
# Before running your code, you should run R and install any required packages.
module load R/4.3.1-icelake

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

# Step 0: run R scripts

RFILE=`ls ../scripts/04_MVMRs_*.R | sed "${SLURM_ARRAY_TASK_ID}q;d"`

R CMD BATCH --vanilla $RFILE $RFILE.out

RFILE=`ls ../scripts/06_MVMRs_*.R | sed "${SLURM_ARRAY_TASK_ID}q;d"`

R CMD BATCH --vanilla $RFILE $RFILE.out

###############################################################
### You should not have to change anything below this line ####
###############################################################

JOBID=$SLURM_JOB_ID

echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then
        echo "Running on nodes: $SLURM_JOB_NODELIST"
else
        echo "Running on node: `hostname`"
fi

echo "Current directory: `pwd`"
echo -e "\nNum tasks = $SLURM_NTASKS, Num nodes = $SLURM_JOB_NUM_NODES, OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo -e "\nExecuting command:\n==================\n$CMD\n"

eval $CMD
