# Freesurfer-based Alzheimer’s Disease Score and Brain Age Score

JPND-BRIDGET Working Package 1.2

S. Frenzel, University Medicine Greifswald, frenzels@uni-greifswald.de, 23/05/19

## Calculation of Freesurfer-based Alzheimer's disease score and brain age score

### Prerequisites
- Linux or Unix computer with Bash shell to run the shell-script. The Bash shell is available on most Linux computers.
- Installation of Freesurfer version 5.3. The corresponding environment variables need to be set up correctly.
- Installation of R with packages glmnet and splines.

### Prepare list of subjects with basic covariates age and sex

At first you'll need to create a list of subject IDs, named **_covariates.csv_**, and put it in the main directory. These IDs need to be the same as the directory names in the Freesurfer directory. Additionally, this file needs to contain the chronological age in the second colum, as well as sex in the third one (numerically encoded, **2=female, 1=male**). The first row should contain the column labels "SubjID", "Age", and "Sex". So, basically, this files looks like

```
SubjID,Age,Sex 
SubjID1,74,2 
SubjID2,55,2 
SubjID3,73,2 
SubjID4,28,1
```

We'll need the basic covariates for adjusting the Freesurfer data for age and sex in order to detect potential outliers.

### Extract data from Freesurfer

Please use the script **_step1_extract_fs_stats.sh_** to extract the Freesurfer data. It again contains the variable SUBJECTS_DIR, which need to be set appropriately. SUBJECTS_DIR is the directory containing the Freesurfer data with each subdirectory corresponding to a subject id as in *covariates.csv*.  After setting the variable you can run the script using the command

```
. step1_extract_fs_stats.sh
```

or

```
bash step1_extract_fs_stats.sh
```

The script extracts tables from the Freesurfer directory and puts them into the subdirectory ./data. After completion there should be the following four files in the data directory:

- lh.aparc.{thickness,volume}.stats
- rh.aparc.{thickness,volume}.stats
- aseg.stats
- wmparc.stats

### Calculate the scores

Start R and install the packages **glmnet** and **splines** if necessary. This can be done by executing the command
```
install.packages(c('splines', 'glmnet'))

```
Compute the scores by executing
```
source("step2_compute_scores.r")
```
 After completion there  should be two CSV files in the directory ./scores containing the Alzheimer's disease score **_./scores/fsad.csv_** and the brain age score **_scores/fsba.csv_**. These files also contain the basic covariates age and sex as provided by _covariates.csv_ and the estmated intracranial volume from Freesurfer.

### Statistical outlier detection

The output files _fsad.csv_ and _fsba.csv_ contain the binary variable **isOutlier** with 1 denoting an outlier. Outliers are defined by regressing each feature on age (non-linearly), sex and intracranial volume (ICV), and looking for residuals outside the four-sigma interval around the mean. We recommend removing samples with isOutlier equal to 1 in subsequent analyses.
       
## GWAS Analysis Plan

Endpoints (dependent variables):
- Alzheimer's disease score (_./scores/fsad.csv_) 
- Brain age score (_./scores/fsba.csv_) 

Exclusion criteria: 
- _Recommended:_ outliers as defined by the isOutlier variable in _./scores/fsad.csv_ and _./scores/fsba.csv_, respectively
- Morphological abnormalities (e.g. cysts, brain tumors etc.)
- Stroke (including lacunes)
- Brain surgery
- Epilepsy before 30 years of age
- Multiple sclerosis

Genotypes
 - genome-wide data imputed to the **HRC 1.1 panel** or **1000 Genomes, phase 1 version 3 (march 2012) ALL populations reference panel
 - imputed X chromosome included – **male alleles coded as 0/2, female 0/1/2**

Covariates
- Age
- Sex + interaction with age 
- Intracranial volume (ICV) as estimated by Freesurfer
- _Optional:_ cohort specific covariates (principal components for population stratification, MR scanner, ...)

GWA Analyses
- Genotype-phenotype association using an additive genetic model, taking genotype uncertainties of imputed SNPs into account
- Linear regressions of FSAD and FSBA 
- Models:
  1. FSAD ~ Age + Age² + Sex + Age\*Sex + ICV + _cohort_specifc_covariates_ + **SNP**
  2. Same as 1. but only cases with Age < 60
  3. Same as 1. but only cases with Age > 60
  4. FSBA ~ Age + Age² + Sex + Age\*Sex + ICV + _cohort_specifc_covariates_ + **SNP**
  5. FSAD ~ Age + Age² + Sex + Age\*Sex + ICV + _cohort_specifc_covariates_ + SNP + **SNP\*Age**
  6. FSBA ~ Age + Age² + Sex + Age\*Sex + ICV + _cohort_specifc_covariates_ + SNP + **SNP\*Age**

Data exchange format and result file upload
- Please do not remove any SNPs after GWAS e.g. based on MAF or imputation quality. This will be done centrally.
- Please name the result file according to the trait analyzed:
```
studyname_stratum_MMDDYYYY.txt
```
where stratum is one of the following: fsad, fsad_old, fsad_young, fsba, fsad_interaction, fsba_interaction
- the results file should be a comma- or tab-separated gzipped text file providing at least the following columns:
  - **SNPID** Please use the SNPID exactly as it is represented in the imputation output. We will convert these to a common ID during data cleaning and meta analysis
  - **chr** chromosome number. Use symbols X, XY, Y and mt for non-autosomal markers. Example: 10
  - **position** physical position for the reference sequence (only build 37/ hg19). Example: 30002
  - **coded_all** coded allele, also called modeled allele (in example of A/G SNP in which AA=0, AG=1 and GG=2, the coded allele is G). Example: C
  - **noncoded_all** the alternate allele. Example: G
  - **strand_genome** + or -, representing either the positive/forward strand or the negative/reverse strand of the human genome reference sequence. Example: +
  - **beta** beta estimate from genotype-phenotype association, at least 5 decimal places -- “NA” if not available. Example: 0.12225
  - **SE** standard error of beta estimate, to at least 5 decimal places -- “NA” if not available. Example: 0.00154
  - **pval** p-value of test statistic, here just as a double check -- “NA” if not available. Example: 0.58652
  - **AF_coded_all** allele frequency for the coded allele -- “NA” if not available. Example: 0.05
  - **n_total** total sample with phenotype and genotype for SNP. Example: 12000
  - **used_for_imp** 1/0 coding; 1=used for imputation, 0=not used for imputation. Example: 1
  - **oevar_imp** Imputation quality; observed divided by expected variance for imputed allele dosage -- NA otherwise
(r2 from MACH, info from IMPUTE2). Example: 0.9685
