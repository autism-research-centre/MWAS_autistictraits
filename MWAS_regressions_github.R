---
title: "Two_step_regression_MWAS"
author: "H Brunel"
date: "1/7/2018"
---

#The following script aims to run an MWAS with a two step strategy consisting on two regression models
#Email hb493@medschl.cam.ac.uk if you have any queries
#The whole procedure takes approximately 12 hours to run on a single-thread and requires a maximum of 20GB 
  
load("data_QC.Rdata") #load the ALSPAC data already pre-processed as data_QC.


# split the dataset into a dataset only containing the cpgs and a dataset with the additional info.
dinfo<-data_QC[,c("Sample_Name", names(data_QC)[482857:482885])]
cpg <- data_QC[, c(names(data_QC)[2:482856])]


#first regression model: 
#normalized methylation probes (betas) against against technical covariates: slide, sample type, and plates and cell counts (Bcell, CD4T, CD8T, Gran, Mono, NK). 
model <- apply(cpg, 2, function(x) 
  residuals(lm(x~ 
                 dinfo[,'sample_type'] + 
                 dinfo[,'Bcell'] +
                 dinfo[,'CD4T'] + 
                 dinfo[,'CD8T'] + 
                 dinfo[,'Gran'] +
                 dinfo[,'Mono'] +
                 dinfo[,'NK'] +
                 dinfo[,'BCD_plate'] +
                 dinfo[,'Slide'])))

#second regression model: scdc scores against residuals from first model and with sex and the first two genetic principal components as covariates.     
m2.model = apply(model, 2, function(x) summary(glm.nb(dinfo$scdc.score ~ x + dinfo$Sex.x + dinfo$PC1 + dinfo$PC2))$coefficients[2,])

