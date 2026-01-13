#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J TrajGWAS_waldTest

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=5:00:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 76 on icelake; 112 on sapphire
#SBATCH --cpus-per-task=1

## Each task is allocated 3.3G (icelake) or 6.7G (icelake-himem) or 4.6G (sapphire)
## If this is insufficient, uncomment and edit this line.
## Maximum value 256G (icelake/sapphire) or 512G (icelake-himem)
#SBATCH --mem=64G

## The system can send emails when your job starts and stops.
## Values include BEGIN, END, ALL, and TIME_LIMIT_80 and TIME_LIMIT_90 
## (reaching 80% or 90% of time limit.) Specify ARRAY_TASKS to receive
## a separate mail for each task. Multiple values can be given, separated by a comma.
#SBATCH --mail-type=FAIL

## The project account name.
## Use mrc-bsu-sl2-cpu for icelake and mrc-bsu-sl2-gpu for ampere
#SBATCH -A mrc-bsu2-sl2-cpu

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
## #SBATCH --array=1-12

##  - - - - - - - - - - - - - -

## Section 2: Modules

# All scripts should include the first three lines.

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

# Load the latest R version.
# Before running your code, you should run R and install any required packages.
module load R/4.3.1-icelake
module load plink/2.00-alpha
module load julia/1.11.4

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

# Step 1: R script (approx 60 min)
# 
# Here I check the SPA results and extract significant SNPs for which I want to get the effect sizes

# R CMD BATCH --vanilla 03_PrepData_WaldTest.R 03_PrepData_WaldTest.R.out
# cp Rplots.pdf 03_PrepData_WaldTest_Rplots.pdf
# rm Rplots.pdf

# Step 2: PLINK & Julia (approx 30 min per setting)

AllSupPops=("men" "old" "young" "women" "post" "pre")

# Create the results folder
mkdir -p ../../rds/hpc-work/MVMR_longTC_CAD/UKB/trajGWAS_results_WaldTest/

#for SubPop in `seq 0 1 5`
for SubPop in `seq 5 1 5`
do 
    myPop=${AllSupPops[$SubPop]}
    #myPop=${AllSupPops[4]}
    echo $myPop
    
    # Copy phenotype file into the temp file 
    cp ../../rds/hpc-work/MVMR_longTC_CAD/UKB/01_Prep_TrajGWAS_$myPop.csv ../../rds/hpc-work/MVMR_longTC_CAD/UKB/temp_trajGWAS/BED/01_Prep_TrajGWAS.csv 
    
    # Get the BED format
    for CHR in `seq 1 1 22`
    do
        #CHR=22
        echo $CHR
        plink2 --pfile ~/rds/hpc-work/MVMR_longTC_CAD/UKB/temp_pgen/UKB_chr$CHR --keep-fam ~/rds/hpc-work/MVMR_longTC_CAD/UKB/01_Prep_TrajGWAS_SampleList_$myPop.txt --extract ~/rds/hpc-work/MVMR_longTC_CAD/UKB/03_candidateSNPs.txt --threads 20 --make-bed --out ~/rds/hpc-work/MVMR_longTC_CAD/UKB/temp_trajGWAS/BED/UKB_chr$CHR
        
    done
    
    # Run Julia script
    if [ $myPop == 'pre' ]
    then
      echo Use different code with selected PCs for variability
      julia 04_b_trajGWAS_WaldTest_youngerWomen.jl
    else
      echo Use standard julia code with all 10 PCs
      julia 04_a_trajGWAS_WaldTest.jl
    fi
    
    # Store result files somewhere else (default everything in the temp/BED)
    mkdir -p ../../rds/hpc-work/MVMR_longTC_CAD/UKB/trajGWAS_results_WaldTest/$myPop/
    mv ../../rds/hpc-work/MVMR_longTC_CAD/UKB/temp_trajGWAS/BED/*.pval.txt ../../rds/hpc-work/MVMR_longTC_CAD/UKB/trajGWAS_results_WaldTest/$myPop/ 
    mv ../../rds/hpc-work/MVMR_longTC_CAD/UKB/temp_trajGWAS/BED/*.null.txt ../../rds/hpc-work/MVMR_longTC_CAD/UKB/trajGWAS_results_WaldTest/$myPop/ 
    
    # Remove temp BED data and pheno file (to make sure nothing gets mixed up)
    rm -r ../../rds/hpc-work/MVMR_longTC_CAD/UKB/temp_trajGWAS/BED/*
        
done

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
