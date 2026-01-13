# MVMR of longitudinal TC on risk of CAD, sex-stratified

## Aim 

Repeat the multivariable Mendelian Randomization (MVMR) analyses of longitudinal total cholesterol (TC) on risk for coronary artery disease (CAD), but in a sex- and age-stratified manner. The reason for this is the different direction of the slope of TC in pre-menopausal (increasing) and post-menopausal (decreasing) women. 

For women, I define age groups by their menopausal status at baseline: 

- **pre-menopausal**: women who report at the UK Biobank (UKB) baseline / follow-up that they have not yet had menopause and are at baseline younger than 61 years. I will use all observations of the electronic health records (EHR) dating before the UKB baseline examination
- **post-menopausal**: women who report at the UKB baseline / follow-up that they have had menopause and are at baseline older than 50 years. I will use all observations of the EHR dating after the UKB baseline examination

For men, I use the same age-threshold. For men between 50 and 60, I randomly put them in either group. 

## Data sets

- **exposure data**: 
    - discovery: UKB baseline, UKB follow-up, and EHR matched data on TC (UKB application number 98032)
    - replication: INTERVAL
- **outcome data**: [Agaram et al. (2022)](https://hugeamp.org/dinspector.html?dataset=Aragam2022_CAD_EU), CAD sex-stratified (EUR only) 
- **instrument selection**: 
    - trajGWAS using the UKB data
    - GLGC hits per sex [Kanoni et al. (2022)](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/results/sex_and_ancestry_specific_summary_stats/), TC sex-stratified (EUR only + only estimates for the mean, not the variability)

## Analysis plans

1) UKB analyses (discovery): 
    - prep data (get exposure, define groups)
    - trajGWAS screening (saddle point approximation, SPA, just p-values)
    - get candidate SNPs (significant in screening or GLGC hits)
    - trajGWAS effect estimation for candidate SNPs
    - GAMLSS effect estimation for candidate SNPs
2) INTERVAL analyses (replication): 
    - get data (get exposure, define groups)
    - trajGWAS effect estimation for candidate SNPs
    - GAMLSS effect estimation for candidate SNPs
3) MVMR:
    - MVMRs per group, using either trajGWAS or GLGC IVs and estimates by trajGWAS or GAMLSS
    - MVMRs per sex, combining the age groups, using either trajGWAS or GLGC IVs and estimates by trajGWAS or GAMLSS
4) Plots and figures: 
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
