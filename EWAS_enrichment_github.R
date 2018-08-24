---
title: "MWAS_autism_SCDCanalysis"
author: "V_warrer"
date: "1/7/2018"
---

# The following code is run in 7 steps. 
# All steps must be repeated if you want to re-run the analysis with a new GWAS
# Steps 2, 4, 5, 6, and 7 must be repeated with different p-value threshold
# Email vw260@medschl.cam.ac.uk if you have any queries
# The permuatation takes approximately 5 hours to run on a machine with 8 GB RAM
  
  
#Step 1: Load the dataframes and merge
  
setwd(" ")#set your working directory

library(data.table)
library(dplyr)
library (plyr)

#ensure plyr is loaded after dplyr

load("SCDCnegbin.RData") #read the mwas data as merged3 
setnames(merged3, "name", "gene") #change name for further merging
merged1 = merged3[,c("gene", "beta_scdc", "pval_scdc", "se_scdc")] #keep only necessary columns

gwas = fread("danerformwaspruningclumped.txt.clumped") #read clumped GWAS data
gwas = gwas[,1:5] #keep only the necessary columns


mqtl = fread("cord.ALL.M.tab") #read cord-blood mqtls as mqtl

#merge all three files using two different merge steps.
merged2 = merge(merged1, mqtl, by = "gene")
merged2 = merge(merged2, gwas, by = "SNP")

```
#Step2:Basic t-test and mean check
#

pthreshold = 0.1 #modify this according to the threshold you need. This is the SCDC MWAS p-value threshold.

x1 = subset(merged2, pval_scdc < pthreshold)
y1 = subset(merged2, pval_scdc > pthreshold)
x2 = x1[!duplicated(x1$SNP),]
y2 = y1[!duplicated(y1$SNP),]

meandiff_check = mean(y2$P) - mean(x2$P)

meandiff_check

t.test(x2$P, y2$P, alternative = "less")



#Step 3a: Count the number of SNPs that map onto each CpG

count1 = count(merged2, "SNP")
merged2 = merge(merged2, count1, by = "SNP")
setnames(merged2, "freq", "SNPcount")


#Step 3b: create bins

first = subset(merged2, SNPcount < 6)
first$category = "1"
second = subset(merged2, SNPcount > 5  & SNPcount < 11)
second$category = "2"
third = subset(merged2, SNPcount > 10 & SNPcount < 16)
third$category = "3"
fourth = subset(merged2, SNPcount > 15 & SNPcount < 21)
fourth$category = "4"
fifth = subset(merged2, SNPcount > 20 & SNPcount < 26)
fifth$category = "5"
sixth = subset(merged2, SNPcount > 25 )
sixth$category = "6"



#Step 4: Count the number of SNPs in each category and merge

one_nominal = subset(first, pval_scdc < pthreshold)
one_notnominal = setdiff(first, one_nominal)

length_one = nrow(one_nominal)

two_nominal = subset(second, pval_scdc < pthreshold )
two_notnominal = setdiff(second, two_nominal)

length_two = nrow(two_nominal)

three_nominal = subset(third, pval_scdc < pthreshold)
three_notnominal = setdiff(third, three_nominal)

length_three = nrow(three_nominal)

four_nominal = subset(fourth, pval_scdc < pthreshold)
four_notnominal = setdiff(fourth, four_nominal)

length_four = nrow(four_nominal)

five_nominal = subset(fifth, pval_scdc < pthreshold)
five_notnominal = setdiff(fifth, five_nominal)

length_five = nrow(five_nominal)

six_nominal = subset(sixth, pval_scdc < pthreshold)
six_notnominal = setdiff(sixth, six_nominal)

length_six = nrow(six_nominal)

nominal_all = list(data.frame(one_nominal), data.frame(two_nominal), data.frame(three_nominal), data.frame(four_nominal),
                   data.frame(five_nominal), data.frame(six_nominal))

notnominal_all = list(data.frame(one_notnominal), data.frame(two_notnominal), data.frame(three_notnominal), data.frame(four_notnominal),
                      data.frame(five_notnominal), data.frame(six_notnominal))




nominal_all_2 = do.call(rbind, nominal_all)
notnominal_all_2 = do.call(rbind, notnominal_all)


#Step 5: Run the analysis
a = nominal_all_2[!duplicated(nominal_all_2$SNP),]
b = notnominal_all_2[!duplicated(notnominal_all_2$SNP),]

meandiff_original = mean(b$P) - mean(a$P)


paste("meandiff_original is", meandiff_original, sep = " ")
paste("meandiff_check is", meandiff_check, sep = " ")
meandiff_original == meandiff_check

##The above statement should be TRUE. If not TRUE, pause, stretch your legs, drink coffee (before 3 PM), and figure out where the bug is##


# Step 6: Permute


meandiff = 1

npermute = 10000
for (i in 1:npermute) {
  # get a random subsample from the background
  one_nominal = sample_n(first, length_one)
  one_notnominal = setdiff(first, one_nominal)
  
  two_nominal = sample_n(second, length_two)
  two_notnominal = setdiff(second, two_nominal)
  
  three_nominal = sample_n(third, length_three)
  three_notnominal = setdiff(third, three_nominal)
  
  four_nominal = sample_n(fourth, length_four)
  four_notnominal = setdiff(fourth, four_nominal)
  
  five_nominal = sample_n(fifth, length_five)
  five_notnominal = setdiff(fifth, five_nominal)
  
  six_nominal = sample_n(sixth, length_six)
  six_notnominal = setdiff(sixth, six_nominal)
  
  
  nominal_all = list(data.frame(one_nominal), data.frame(two_nominal), data.frame(three_nominal),    data.frame(four_nominal), data.frame(five_nominal), data.frame(six_nominal))
  
  notnominal_all = list(data.frame(one_notnominal), data.frame(two_notnominal),data.frame(three_notnominal), data.frame(four_notnominal), data.frame(five_notnominal), data.frame(six_notnominal))
  
  nominal_all_2 = do.call(rbind, nominal_all)
  notnominal_all_2 = do.call(rbind, notnominal_all)
  
  a = nominal_all_2[!duplicated(nominal_all_2$SNP),]
  b = notnominal_all_2[!duplicated(notnominal_all_2$SNP),]
  
  meandiff[i] = mean(b$P) - mean(a$P)
}

```
#Step 7: Run enrichment and generate plots

m = (sum(meandiff > meandiff_original) + 1) /npermute

hist(meandiff)
abline(v=meandiff_original, lwd=2, col="purple")

m
