#!/bin/bash
# --------------------------------
# JPND-BRIDGET WP 1.2, NeuroCHARGE
# --------------------------------
# Freesurfer-based AD score and brain age score - Step 1
#
# extracts tables from $SUBJECT_DIR/*/stats for each subject and puts them
# together into a single table
#
# S. Frenzel, University Medicine Greifswald, 21/05/19
# frenzels@uni-greifswald.de


# directory containing freesurfer data
export SUBJECTS_DIR='???'

# CSV with the first column containing the subject ids, the second column containing
# chronologcial age, and the third column containing sex (2=female, 1=male). The first
# row contains the column labels.
export COV_FILE='covariates.csv'


##### in theory there should be no need to change anything below this line... #####


# perform some checks
[[ -d ${SUBJECTS_DIR} ]] \
|| ( echo 'The Freesurfer data directory you specified does not exist! Please'\
              'check the SUBJECTS_DIR variable.'
     exit )

[[ -f ${COV_FILE} ]] \
|| ( echo 'The covariate file you specified does not exist!'
     exit )

mkdir -p data

# build list of subjects from $COV_FILE
subj_list=""
for subj in $(cut -f1 -d, ${COV_FILE} | tail -n+2); do
	subj_list+=" ${subj}"
done


# create table containing volumes of subcortical structures
echo '----------------------------------------------------------'
echo 'Extract subcortical data...'
asegstats2table --subjects ${subj_list} --skip --tablefile data/aseg.stats


# create table containing cortical thicknesses and volumes
echo '----------------------------------------------------------'
echo 'Extract cortical data...'
for hemi in {lh,rh}; do
# thickness
aparcstats2table --subjects ${subj_list} --parc aparc --skip --hemi ${hemi} \
	--meas thickness --tablefile data/${hemi}.aparc.thickness.stats
# volume
aparcstats2table --subjects ${subj_list} --parc aparc --skip --hemi ${hemi} \
	--meas volume --tablefile data/${hemi}.aparc.volume.stats
done

# create table containing the standard (dmax=5mm) white matter parcellation
echo '----------------------------------------------------------'
echo 'Extract white matter parcellation with dmax=5mm (default)...'
asegstats2table --all-segs --subjects ${subj_list} --stats wmparc.stats --skip\
	--tablefile data/wmparc.stats
