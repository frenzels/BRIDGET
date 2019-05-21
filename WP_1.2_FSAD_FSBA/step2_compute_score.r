#!/usr/bin/env Rscript
# --------------------------------
# JPND-BRIDGET WP 1.2, NeuroCHARGE
# --------------------------------
# Freesurfer-based AD score and brain age score - Step 2
#
# S. Frenzel, University Medicine Greifswald, 21/05/19
# frenzels@uni-greifswald.de

library(glmnet)

library(splines)

loadData <- function (directory, measure, fullWMParc=F) {
	sanitizeColumnLabels <- function (x) {
		names(x)[1] <- 'SubjID'	
		names(x) <- gsub('.', '_', names(x), fixed=T)	
		names(x) <- sub('^lh_', 'L_', names(x))	 
		names(x) <- sub('^rh_', 'R_', names(x)) 
		names(x) <- sub('^lh', 'L_', names(x)) # in dieser Reihenfolge!	 
		names(x) <- sub('^rh', 'R_', names(x)) 
		names(x) <- sub('^wm_lh_', 'WM_L_', names(x))	 
		names(x) <- sub('^wm_rh_', 'WM_R_', names(x))
		names(x) <- sub('^Left_', 'L_', names(x))	 
		names(x) <- sub('^Right_', 'R_', names(x))
		return(x)	
	}
	## load cortical data ##
	xl <- read.table(paste(directory, '/lh.aparc.', measure, '.stats', sep=''), header=T)
	xl <- sanitizeColumnLabels(xl)
	xr <- read.table(paste(directory, '/rh.aparc.', measure, '.stats', sep=''), header=T)	
	xr <- sanitizeColumnLabels(xr)
	x <- merge(xl, xr, by='SubjID')
	## load white matter parcellation ##
	xx <- read.table(
		paste(directory, '/wmparc.', ifelse(fullWMParc,'dmax999.',''), 'stats', sep=''),
		header=T
	)
	xx <- sanitizeColumnLabels(xx)
	# already contained in aseg.stats:
	xx[,'MaskVol'] <- NULL
	xx[,'EstimatedTotalIntraCranialVol'] <- NULL
	# merge
	x <- merge(x, xx, by='SubjID')  # Daten mergen nach Probanden ID
	rm(xx)
	## load subcortical volumes ##
	xx <- read.table(paste(directory, '/aseg.stats', sep=''), header=T)
	xx <- sanitizeColumnLabels(xx)
	# already contained in wmparc.stats
	xx[,'L_CorticalWhiteMatterVol'] <- NULL
	xx[,'R_CorticalWhiteMatterVol'] <- NULL
	xx[,'CorticalWhiteMatterVol'] <- NULL
	# merge
	x <- merge(x, xx, by='SubjID')  # Daten mergen nach Probanden ID
	x <- sanitizeColumnLabels(x)
	rm(xx)
	#
	x[,'ICV']  <- x[,'EstimatedTotalIntraCranialVol']
	return(x)
}

detectOutliers <- function(x, featureLabels) {
	univariateOutliers <- function(x, colNames, form, C, logTrans=F) {
		# calculate residuals
		xx <- x
		if (logTrans)
			xx[,colNames] <- log(xx[,colNames])
		for (colName in colNames) {
			model <- lm(paste(colName, '~', form), xx)
		  	xx[,colName] <- residuals(model)
		}	
		# detect outliers
		outlier <- c()
		for (colName in colNames) {
			m <- mean(xx[, colName])
			sigma <- sd(xx[, colName]) #0%, 25%, 50%, 75%, 100%
			outlier <- c(outlier, which(xx[, colName] > (m + C*sigma)))
			outlier <- c(outlier, which(xx[, colName] < (m - C*sigma)))
		}
		return(unique(sort(outlier)))
	}
	vent <- which(featureLabels %in% c('L_Lateral_Ventricle','L_Inf_Lat_Vent','X3rd_Ventricle',
	'X4th_Ventricle', 'L_choroid_plexus', 'R_Lateral_Ventricle', 'R_Inf_Lat_Vent',
	'R_choroid_plexus'))
	form <- 'ns(Age,df=5)*factor(Sex) + ICV' 	
	threshold <- 4
	outliers <- univariateOutliers(x, featureLabels[vent], form, 4, logTrans=T)
	outliers <- c(outliers, univariateOutliers(x, featureLabels[-vent], form, 4))
	return(unique(sort(outliers)))

}

#########################################################################################
## FSAD
#########################################################################################

# load fsad model
fsad <- new.env()
load('models/fsad.rdata', envir=fsad)

# load data
x <- loadData('data', 'thickness')
covars <- read.csv('covariates.csv')
if (!identical(x$SubjID,covars$SubjID)) {
	message('ERROR: this should not happen :-( ...')
	return()
}
x <- merge(x, covars, by='SubjID')

# statistical outlier detection
x$isOutlier <- as.numeric((1:nrow(x)) %in% detectOutliers(x, fsad$featureLabels))

# standardize data
xx <- x[,fsad$featureLabels]
xx <- scale(xx, center=fsad$preproc$mean, scale=fsad$preproc$std)
x[,fsad$featureLabels] <- data.frame(xx)

# predict ad scores
x$FSAD <- -predict.glmnet(fsad$model, as.matrix(x[,fsad$featureLabels]), s=as.numeric(fsad$bestTune[2]))

# save it
write.csv(x[,c('SubjID','Age','Sex','ICV','FSAD','isOutlier')], file='scores/fsad.csv', row.names=F, quote=F)

# clean up
rm(fsad,x,covars,xx)

#########################################################################################
## FSBA
#########################################################################################

# load fsad model
fsba <- new.env()
load('models/fsba.rdata', envir=fsba)

# load data
x <- loadData('data', 'volume')
covars <- read.csv('covariates.csv')
if (!identical(x$SubjID,covars$SubjID)) {
	message('ERROR: this should not happen :-( ...')
	return()
}
x <- merge(x, covars, by='SubjID')

# statistical outlier detection
x$isOutlier <- as.numeric((1:nrow(x)) %in% detectOutliers(x, fsba$featureLabels))

# predict ba scores
x$FSBA <- rep(NA, nrow(x))
x$FSBA[covars$Sex==1] <- as.matrix(x[covars$Sex==1,featureLabels]) %*% fsba$modelMale[featureLabels] +
		fsba$modelMale["(Intercept)"]
x$FSBA[covars$Sex==2] <- as.matrix(x[covars$Sex==2,featureLabels]) %*% fsba$modelFemale[featureLabels] +
		fsba$modelFemale["(Intercept)"]

# save it
write.csv(x[,c('SubjID','Age','Sex','ICV','FSBA','isOutlier')], file='scores/fsba.csv', row.names=F, quote=F)

# clean up
rm(fsba,x,covars)

####
message('Done. Scores written to scores/fsad.csv and scores/fsba.csv :-)')
