#!/bin/bash

## Example SLURM script for BSU icelake jobs

## Section 1: SLURM Commands

## All SLURM commands must be placed at the start of the file
## Full documentation can be found here: https://slurm.schedmd.com/sbatch.html

## Enter a short name for the job, to be shown in SLURM output
#SBATCH -J 02_Prep_SNP_PLINK

## Enter the wall-clock time limit for your jobs.
## If jobs reach this limit they are automatically killed.
## Maximum value 36:00:00.
#SBATCH --time=06:00:00

## For single-core jobs, this number should be '1'. 
## If your job has built-in parallelism, eg using OpenMP or 
## R's foreach() and doParallel(), increase this number as desired.
## The maximum value is 76 on icelake; 112 on sapphire
#SBATCH --cpus-per-task=20

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
module load bgen/1.2

# If using the GPU cluster, replace the third line with the uncommented line:
# module load rhel8/default-amp

#! Insert additional module load commands after this line if needed:

## - - - - - - - - - - -

## Section 3: Run your application

# Step 0: run R scripts
R CMD BATCH --vanilla ../scripts/02_SNP_01_checkCandidateSNPs.R ../scripts/02_SNP_01_checkCandidateSNPs.R.out

# Step 1: create temporary folder for bgens
mkdir ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat/temp

# Step 2: extract SNPs and samples per chromosome and save as bgen
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c1_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c1_b0_v3_s487160.sample --chr 1 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr1
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c2_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c2_b0_v3_s487160.sample --chr 2 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr2
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c3_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c3_b0_v3_s487160.sample --chr 3 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr3
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c4_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c4_b0_v3_s487160.sample --chr 4 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr4
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c5_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c5_b0_v3_s487160.sample --chr 5 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr5
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c6_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c6_b0_v3_s487160.sample --chr 6 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr6
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c7_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c7_b0_v3_s487160.sample --chr 7 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr7
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c8_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c8_b0_v3_s487160.sample --chr 8 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr8
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c9_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c9_b0_v3_s487160.sample --chr 9 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr9
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c10_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c10_b0_v3_s487160.sample --chr 10 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr10
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c11_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c11_b0_v3_s487160.sample --chr 11 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr11
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c12_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c12_b0_v3_s487160.sample --chr 12 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr12
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c13_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c13_b0_v3_s487160.sample --chr 13 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr13
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c14_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c14_b0_v3_s487160.sample --chr 14 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr14
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c15_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c15_b0_v3_s487160.sample --chr 15 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr15
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c16_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c16_b0_v3_s487160.sample --chr 16 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr16
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c17_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c17_b0_v3_s487160.sample --chr 17 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr17
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c18_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c18_b0_v3_s487160.sample --chr 18 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr18
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c19_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c19_b0_v3_s487160.sample --chr 19 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr19
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c20_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c20_b0_v3_s487160.sample --chr 20 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr20
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c21_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c21_b0_v3_s487160.sample --chr 21 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr21
plink2 --bgen ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed//ukb22828_c22_b0_v3.bgen 'ref-last' --sample ~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/genotypes-imputed/ukb22828_c22_b0_v3_s487160.sample --chr 22 --keep-fam ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//01_Prep_01_SampleList_CAD.txt --extract ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_SNPList_GLGC.txt --mach-r2-filter 0.8 2 --maf 0.01 --threads 20 --export bgen-1.2 bits=8 id-delim='-' --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr22

# Step 3: merge bgens
cat-bgen -clobber -g ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr1.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr2.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr3.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr4.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr5.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr6.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr7.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr8.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr9.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr10.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr11.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr12.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr13.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr14.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr15.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr16.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr17.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr18.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr19.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr20.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr21.bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr22.bgen -og ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_merged.bgen

# Step 4: make pgen
plink2 --bgen ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_merged.bgen 'ref-last' --sample ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp/UKB_TC_GLGC_chr1.sample --make-pgen --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_UKB_TC_GLGC_merged

# Step 5: remove temporary files
rm -rf ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//temp

# Step 6: get allele frequencies in pgen data
plink2 --pfile ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_UKB_TC_GLGC_merged --freq --out ~/rds/hpc-work/data/MVMR_longTC_CAD_sexStrat//02_SNP_01_UKB_TC_GLGC_merged_AF

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
