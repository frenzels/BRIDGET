##Freesurfer-based Alzheimer’s Disease Score and Brain Age Score

JPND-BRIDGET Working Package 1.2

S. Frenzel, University Medicine Greifswald, frenzels@uni-greifswald.de
23/05/19


### Prerequisites
- Linux or Unix computer with Bash shell to run the shell-script. The Bash shell is available on most Linux computers.
- Installation of Freesurfer version 5.3. The corresponding environment variables need to be set up correctly.
- Installation of R with packages glmnet and splines.

### Prepare list of subjects with basic covariates age and sex
At first you'll need to create a list of subjects as csv file in the main directory, named covariates.csv.  These need to be the same as the directory names in the Freesurfer directory. Additionally, this file needs to contain the chronological age in the second colum, as well as sex in the third one (numerically encoded, 2=female, 1=male). The first row should contain the column labels SubjID, Age, and Sex. So, basically, this files looks like...
SubjID,Age,Sex 
SubjID1,74,2 
SubjID2,55,2 
SubjID3,73,2 
SubjID4,28,1
...

We'll need the basic covariates for adjusting the Freesurfer data for age and sex in order to detect possible outliers.
3. Calculation of Freesurfer-based AD score and brain age score
    1. Extract the Freesurfer tables:
       Please use the script step1_extract_fs_stats.sh to extract the Freesurfer data. It again contains the variable SUBJECTS_DIR, which need to be set to appropriately. SUBJECTS_DIR is the directory containing the Freesurfer data with each subdirectory corresponding to a subject id in covariates.csv.  After setting the variable you can run the script using the command 
       bash step1_extract_fs_stats.sh
       It extracts tables from the Freesurfer directory and puts them into the directory “data”. After completion there should be the following four files in the “data” directory:
    • lh.aparc.{thickness,volume}.stats
    • rh.aparc.{thickness,volume}.stats
    • aseg.stats
    • wmparc.stats
      
### Calculation of AD score and brain age score:
       Start R and install the packages glmnet and splines if necessary. This can be done by executing the command
       install.packages(c('splines', 'glmnet'))
       Compute the scores by executing 
       source("step2_compute_scores.r")
       After completion there  should be two CSV files in the directory “scores” containing the AD score (fsad.csv) and the brain age score (fsba.csv). These files again contain the basic covariates age and sex, as well as the binary variable isOutlier with 1 indicating a possible outlier.
       

