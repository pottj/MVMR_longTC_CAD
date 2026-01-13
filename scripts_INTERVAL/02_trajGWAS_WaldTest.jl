# Julia Cheat Sheet: https://cheatsheet.juliadocs.org/
# 
# This is a code for the trajGWAS analysis
# Not yet completely sure how it will work

# Load necessary packages
using Pkg, TrajGWAS, CSV

# Define input path for genotypes and phenotypes
datadir = "../../rds/hpc-work/MVMR_longTC_CAD/INTERVAL/temp_trajGWAS/BED/"

# Step 1: fit the null model
nm = trajgwas(@formula(CHOL ~ 1 + age + age2 + med_lipid + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10),   
              @formula(CHOL ~ 1),								
              @formula(CHOL ~ 1 + age + age2 + med_lipid + PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10),	
              :genID,											
              datadir * "01_Prep_TrajGWAS.csv",		
              nothing,
              nullfile = datadir * "trajgwas.null.txt")
              
# Step 2: GWAS for each chromosome 
# this part can be submitted as separate jobs
for chr in 1:22
    #plinkfile = datadir * "INTERVAL_chr22" 
    plinkfile = datadir * "INTERVAL_chr" * string(chr)
    pvalfile = plinkfile * ".pval.txt" 
    trajgwas(nm, plinkfile, pvalfile=pvalfile,test = :wald)
end