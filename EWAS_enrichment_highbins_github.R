---
  title: "MWAS_autism_SCDCanalysis"
author: "V_warrer"
date: "1/7/2018"
output:
  pdf_document: default
html_document: default
---
  
  #Step 1: Load the dataframes and merge
  
  setwd(" ")#set your working directory

library(data.table)
library(dplyr)
library (plyr)

load("SCDCnegbin.RData") #read the mwas data as merged3 
setnames(merged3, "name", "gene") #change name for further merging
merged1 = merged3[,c("gene", "beta_scdc", "pval_scdc", "se_scdc")] #keep only necessary columns

gwas = fread("scdcclumped.clumped") #read clumped GWAS data
gwas = gwas[,1:5] #keep only the necessary columns


mqtl = fread("cord.ALL.M.tab") #read cord-blood mqtls as mqtl

#merge all three files using two different merge steps.
merged2 = merge(merged1, mqtl, by = "gene")
merged2 = merge(merged2, gwas, by = "SNP")

rm(merged3)
rm(mqtl)
rm(gwas)
rm(merged1)

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

first = subset(merged2, SNPcount == 1)
first$category = "1"
second = subset(merged2, SNPcount == 2)
second$category = "2"
third = subset(merged2, SNPcount == 3)
third$category = "3"
fourth = subset(merged2, SNPcount == 4)
fourth$category = "4"
fifth = subset(merged2, SNPcount == 5)
fifth$category = "5"
sixth = subset(merged2, SNPcount == 6 )
sixth$category = "6"
seventh = subset(merged2, SNPcount == 7 )
seventh$category = "7"
eighth = subset(merged2, SNPcount == 8 )
eighth$category = "8"
ninth = subset(merged2, SNPcount == 9 )
ninth$category = "9"
tenth = subset(merged2, SNPcount == 10 )
tenth$category = "10"
eleventh = subset(merged2, SNPcount == 11 )
eleventh$category = "11"
twelfth = subset(merged2, SNPcount == 12 )
twelfth$category = "12"
thirteenth = subset(merged2, SNPcount == 13 )
thirteenth$category = "13"
fourteenth = subset(merged2, SNPcount == 14 )
fourteenth$category = "14"
fifteenth = subset(merged2, SNPcount > 14 & SNPcount < 21 )
fifteenth$category = "15"
sixteenth = subset(merged2, SNPcount > 20 & SNPcount < 31 )
sixteenth$category = "16"
seventeenth = subset(merged2, SNPcount > 30 & SNPcount < 51 )
seventeenth$category = "17"
eighteenth = subset(merged2, SNPcount > 50 & SNPcount < 101 )
eighteenth$category = "18"
nineteenth = subset(merged2, SNPcount > 100 )
nineteenth$category = "19"

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

seventh_nominal = subset(seventh, pval_scdc < pthreshold)
seventh_notnominal = setdiff(seventh, seventh_nominal)

length_seventh = nrow(seventh_nominal)

eighth_nominal = subset(eighth, pval_scdc < pthreshold)
eighth_notnominal = setdiff(eighth, eighth_nominal)

length_eighth = nrow(eighth_nominal)

ninth_nominal = subset(ninth, pval_scdc < pthreshold)
ninth_notnominal = setdiff(ninth, ninth_nominal)

length_ninth = nrow(ninth_nominal)

tenth_nominal = subset(tenth, pval_scdc < pthreshold)
tenth_notnominal = setdiff(tenth, tenth_nominal)

length_tenth = nrow(tenth_nominal)

eleventh_nominal = subset(eleventh, pval_scdc < pthreshold)
eleventh_notnominal = setdiff(eleventh, eleventh_nominal)

length_eleventh = nrow(eleventh_nominal)

twelfth_nominal = subset(twelfth, pval_scdc < pthreshold)
twelfth_notnominal = setdiff(twelfth, twelfth_nominal)

length_twelfth = nrow(twelfth_nominal)

thirteenth_nominal = subset(thirteenth, pval_scdc < pthreshold)
thirteenth_notnominal = setdiff(thirteenth, thirteenth_nominal)

length_thirteenth = nrow(thirteenth_nominal)

fourteenth_nominal = subset(fourteenth, pval_scdc < pthreshold)
fourteenth_notnominal = setdiff(fourteenth, fourteenth_nominal)

length_fourteenth = nrow(fourteenth_nominal)

fifteenth_nominal = subset(fifteenth, pval_scdc < pthreshold)
fifteenth_notnominal = setdiff(fifteenth, fifteenth_nominal)

length_fifteenth = nrow(fifteenth_nominal)

sixteenth_nominal = subset(sixteenth, pval_scdc < pthreshold)
sixteenth_notnominal = setdiff(sixteenth, sixteenth_nominal)

length_sixteenth = nrow(sixteenth_nominal)

seventeenth_nominal = subset(seventeenth, pval_scdc < pthreshold)
seventeenth_notnominal = setdiff(seventeenth, seventeenth_nominal)

length_seventeenth = nrow(seventeenth_nominal)


eighteenth_nominal = subset(eighteenth, pval_scdc < pthreshold)
eighteenth_notnominal = setdiff(eighteenth, eighteenth_nominal)

length_eighteenth = nrow(eighteenth_nominal)

nineteenth_nominal = subset(nineteenth, pval_scdc < pthreshold)
nineteenth_notnominal = setdiff(nineteenth, nineteenth_nominal)

length_nineteenth = nrow(nineteenth_nominal)


nominal_all = list(data.frame(one_nominal), data.frame(two_nominal), data.frame(three_nominal), data.frame(four_nominal),
                   data.frame(five_nominal), data.frame(six_nominal), data.frame(seventh_nominal), data.frame(eighth_nominal),
                   data.frame(ninth_nominal), data.frame(tenth_nominal), data.frame(eleventh_nominal), data.frame(twelfth_nominal),
                   data.frame(thirteenth_nominal), data.frame(fourteenth_nominal), data.frame(fifteenth_nominal),
                   data.frame(sixteenth_nominal), data.frame(seventeenth_nominal), data.frame(eighteenth_nominal), data.frame(nineteenth_nominal))

notnominal_all = list(data.frame(one_notnominal), data.frame(two_notnominal), data.frame(three_notnominal), data.frame(four_notnominal),
                      data.frame(five_notnominal), data.frame(six_notnominal), data.frame(seventh_notnominal), data.frame(eighth_notnominal),
                      data.frame(ninth_notnominal), data.frame(tenth_notnominal), data.frame(eleventh_notnominal), data.frame(twelfth_notnominal),
                      data.frame(thirteenth_notnominal), data.frame(fourteenth_notnominal), data.frame(fifteenth_notnominal),
                      data.frame(sixteenth_notnominal), data.frame(seventeenth_notnominal), data.frame(eighteenth_notnominal), data.frame(nineteenth_notnominal))





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


#Permute


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
  
  seventh_nominal = sample_n(seventh, length_seventh)
  seventh_notnominal = setdiff(seventh, seventh_nominal)
  
  eighth_nominal = sample_n(eighth, length_eighth)
  eighth_notnominal = setdiff(eighth, eighth_nominal)
  
  ninth_nominal = sample_n(ninth, length_ninth)
  ninth_notnominal = setdiff(ninth, ninth_nominal)
  
  tenth_nominal = sample_n(tenth, length_tenth)
  tenth_notnominal = setdiff(tenth, tenth_nominal)
  
  eleventh_nominal = sample_n(eleventh, length_eleventh)
  eleventh_notnominal = setdiff(eleventh, eleventh_nominal)
  
  twelfth_nominal = sample_n(twelfth, length_twelfth)
  twelfth_notnominal = setdiff(twelfth, twelfth_nominal)
  
  
  thirteenth_nominal = sample_n(thirteenth, length_thirteenth)
  thirteenth_notnominal = setdiff(thirteenth, thirteenth_nominal)
  
  fourteenth_nominal = sample_n(fourteenth, length_fourteenth)
  fourteenth_notnominal = setdiff(fourteenth, fourteenth_nominal) 
  
  fifteenth_nominal = sample_n(fifteenth, length_fifteenth)
  fifteenth_notnominal = setdiff(fifteenth, fifteenth_nominal)
  
  sixteenth_nominal = sample_n(sixteenth, length_sixteenth)
  sixteenth_notnominal = setdiff(sixteenth, sixteenth_nominal)
  
  seventeenth_nominal = sample_n(seventeenth, length_seventeenth)
  seventeenth_notnominal = setdiff(seventeenth, seventeenth_nominal)
  
  eighteenth_nominal = sample_n(eighteenth, length_eighteenth)
  eighteenth_notnominal = setdiff(eighteenth, eighteenth_nominal)
  
  nineteenth_nominal = sample_n(nineteenth, length_nineteenth)
  nineteenth_notnominal = setdiff(nineteenth, nineteenth_nominal)
  
  nominal_all = list(data.frame(one_nominal), data.frame(two_nominal), data.frame(three_nominal), data.frame(four_nominal),
                     data.frame(five_nominal), data.frame(six_nominal), data.frame(seventh_nominal), data.frame(eighth_nominal),
                     data.frame(ninth_nominal), data.frame(tenth_nominal), data.frame(eleventh_nominal), data.frame(twelfth_nominal),
                     data.frame(thirteenth_nominal), data.frame(fourteenth_nominal), data.frame(fifteenth_nominal),
                     data.frame(sixteenth_nominal), data.frame(seventeenth_nominal), data.frame(eighteenth_nominal), data.frame(nineteenth_nominal))
  
  notnominal_all = list(data.frame(one_notnominal), data.frame(two_notnominal), data.frame(three_notnominal), data.frame(four_notnominal),
                        data.frame(five_notnominal), data.frame(six_notnominal), data.frame(seventh_notnominal), data.frame(eighth_notnominal),
                        data.frame(ninth_notnominal), data.frame(tenth_notnominal), data.frame(eleventh_notnominal), data.frame(twelfth_notnominal),
                        data.frame(thirteenth_notnominal), data.frame(fourteenth_notnominal), data.frame(fifteenth_notnominal),
                        data.frame(sixteenth_notnominal), data.frame(seventeenth_notnominal), data.frame(eighteenth_notnominal), data.frame(nineteenth_notnominal))
  
  nominal_all_2 = do.call(rbind, nominal_all)
  notnominal_all_2 = do.call(rbind, notnominal_all)
  
  a = nominal_all_2[!duplicated(nominal_all_2$SNP),]
  b = notnominal_all_2[!duplicated(notnominal_all_2$SNP),]
  
  meandiff[i] = mean(b$P) - mean(a$P)
}



#Run enrichment and generate plots

m = (sum(meandiff > meandiff_original) + 1) /npermute


hist(meandiff)
abline(v=meandiff_original, lwd=2, col="purple")

m