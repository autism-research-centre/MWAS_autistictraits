---
title: "Two_step_regression_MWAS"
author: "H Brunel"
date: "1/7/2017"
---

#The following script aims to run an MWAS with a two step strategy consisting on two regression models
#It also annotates the results and draw a Manhattan plot with the results
#Email hb493@medschl.cam.ac.uk if you have any queries
#The whole procedure takes approximately 14.5 hours to run on a single-thread and requires a maximum of 20GB 


library(MASS)
library(meffil)
library(qqman)

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

#Annotation of the results using the "meffil" cpg annotations. 
                 
df <- as.data.frame(m2.model)
res <- data.frame("name"=names(df), "pval" =as.vector(unlist(df[4,])), "est"= as.vector(unlist(df[1,])), "Z"=as.vector(unlist(df[3,])))
res_annot <- merge(res, y)
                 
 # Manhattan plot
                 
#Manhattan plot:

#It requires an intermediate data frame (x) containing cpg name, chromosome and position in a given format. 
x <- res_annot[,c("name", "pval", "chromosome", "position")]
x$chr <- do.call(rbind, strsplit(x$chromosome, split="r"))[,2]
x$pos <- as.numeric(x$position)

#CHR and BP must be numeric
x$chr[which(x$chr=="X")] <- "23"
x$chr[which(x$chr=="Y")] <- "24"
x$chr <-as.numeric(x$chr)
x$pos <-as.numeric(x$pos)

#Manhattan plot                 
png("ManhattanPlot.png")
manhattan(x, chr = "chr", bp = "pos", p = "pval", snp = "name", col = c("gray10", "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
highlight = NULL, logp = TRUE, annotatePval = NULL, annotateTop = TRUE)
dev.off()                 
