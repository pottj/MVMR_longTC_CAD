#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J TrajGWAS_genomewide

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=6:00:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 76 on icelake; 112 on sapphire
#SBATCH --cpus-per-task=1

## Each task is allocated 3.3G (icelake) or 6.7G (icelake-himem) or 4.6G (sapphire)
## If this is insufficient, uncomment and edit this line.
## Maximum value 256G (icelake/sapphire) or 512G (icelake-himem)
#SBATCH --mem=128G

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

# Step 1: R script (approx 1 min)
# 
# Here I generate the trajGWAS input files for 6 subpopulations (men and women, young and old): 
# 
# - 01_Prep_TrajGWAS_men.csv
# - 01_Prep_TrajGWAS_women.csv
# - 01_Prep_TrajGWAS_old.csv
# - 01_Prep_TrajGWAS_post.csv
# - 01_Prep_TrajGWAS_young.csv
# - 01_Prep_TrajGWAS_pre.csv
# 
# and one file containing all samples: 01_Prep_TrajGWAS_SampleList.txt 

R CMD BATCH --vanilla 01_PrepData.R 01_PrepData.R.out
cp Rplots.pdf 01_PrepData_Rplots.pdf
rm Rplots.pdf

# Step 2: PLINK & Julia
#
# Step 2.1: Create PGEN for each chromosome (keep for later reuse!) (~ 6:30h for all chromosomes)
# Get PGEN for all samples
mkdir -p ~/rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_pgen

for CHR in `seq 1 1 22`
do
    echo $CHR
    plink2 --pfile ~/rds/rds-post_qc_data-pNR2rM6BWWA/interval/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_${CHR}_interval --keep-fam ~/rds/hpc-work/MVMR_longTC_CAD/INTERVAL/01_Prep_TrajGWAS_SampleList.txt --extract ~/rds/hpc-work/MVMR_longTC_CAD/UKB/03_candidateSNPs.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --make-pgen --out ~/rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_pgen/INTERVAL_chr$CHR
done

# Step 2.2: Create BED per subpopulation

AllSupPops=("men" "old" "young" "women" "post" "pre")

# Create temporary subfolder for the BED files and the results folder

mkdir -p ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED
mkdir -p ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/trajGWAS_results_WaldTest/

for SubPop in `seq 0 1 5`
do 
    myPop=${AllSupPops[$SubPop]}
    #myPop=${AllSupPops[0]}
    echo $myPop
    
    # Copy phenotype file into the temp file 
    cp ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/01_Prep_TrajGWAS_$myPop.csv ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED/01_Prep_TrajGWAS.csv 
    
    # Get the BED format
    for CHR in `seq 1 1 22`
    do
        echo $CHR
        plink2 --pfile ~/rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_pgen/INTERVAL_chr$CHR --keep-fam ~/rds/hpc-work/MVMR_longTC_CAD/INTERVAL/01_Prep_TrajGWAS_SampleList_$myPop.txt --threads 20 --make-bed --out ~/rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED/INTERVAL_chr$CHR
    done
    
    # Run Julia script
    echo Use standard julia code with 10 PCs
    julia 02_trajGWAS_WaldTest.jl
    
    # Store result files somewhere else (default everything in the temp/BED)
    mkdir -p ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/trajGWAS_results_WaldTest/$myPop/
    mv ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED/*.pval.txt ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/trajGWAS_results_WaldTest/$myPop/ 
    mv ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED/*.null.txt ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/trajGWAS_results_WaldTest/$myPop/ 
    
    # Remove temp BED data and pheno file (to make sure nothing gets mixed up)
    rm -r ../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED/*
    
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
