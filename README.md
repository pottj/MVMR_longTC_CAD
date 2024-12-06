# MVMR of longitudinal TC on risk of CAD, sex-stratified

## Aim 

Repeat the multivariable Mendelian Randomization (MVMR) analyses of longitudinal total cholesterol (TC) on risk for coronary artery disease (CAD), but in a sex-stratified manner. In addition, I want to split the women by menopausal status: 

- set 1: women who report at the UK Biobank (UKB) baseline / follow-up that they have not yet had menopause and are at baseline younger than 61 years. I will use all observations of the electronic health records (EHR) dating before the UKB baseline examination
- set 2: women who report at the UKB baseline / follow-up that they have had menopause and are at baseline older than 50 years. I will use all observations of the EHR dating after the UKB baseline examination

## Data sets

- **exposure data**: UKB baseline, UKB follow-up, and EHR matched data on TC (UKB application number 98032)
- **outcome data**: [Agaram et al. (2022)](https://hugeamp.org/dinspector.html?dataset=Aragam2022_CAD_EU), CAD sex-stratified (EUR only) 
- **instrument selection**: [Kanoni et al. (2022)](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/sex_and_ancestry_specific_summary_stats/), TC sex-stratified (EUR only)

## Analysis plans

1) Data preparation: 
    - get the exposure phenotype data from the UKB and the EHR 
    - specify the subsets for pre- and postmenopausal women
    - test the none-genetic longitudinal model, adjusted for age and lipid-lowering medication, including random effects on the intercept, slope, and variability
2) SNP preparation: 
    - select genome-wide significant and independent SNPs per sex
    - check for overlapping loci and test if they have a shared signal (colocalization)
    - extract gene dosages of selected instruments from the UKB
3) GAMLSS: 
    - run regression model for men, and (post- and pre-menopausal) women
4) MVMR:
    - run MVMR and sensitivity checks 
5) Plots and figures: 
    - create all relevant plots and figures of the project
    
## Abbreviations

- UKB, UK Biobank
- CAD, Coronary Artery Disease
- CARDIoGRAMplusC4D, Coronary ARtery DIsease Genome wide Replication and Meta-analysis (CARDIoGRAM) plus The Coronary Artery Disease (C4D) Genetics consortium
- EHR, Electronic Health Records
- GAMLSS, Generalized Additive Models for Location, Scale and Shape
- GLGC, Global Lipids Genetics Consortium
- MVMR, MultiVariable Mendelian Randomization
- SNP, Single Nucleotide Polymorphism
- TC, Total Cholesterol
