library(data.table)
library(tidyverse)

#load data
load("./blood.RData")
load("./aa.RData")
load("./lv.RData")
load("./muscle.RData")

pheno      <- read.delim("./phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt", sep = "\t",fill = TRUE)
samp <- fread('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')

cbPalette <- c("#009E73", "#D55E00", "#0072B2", "#CC79A7")

coding.region <- c(3307:4262, 4470:5511, 5904:7445, 7586:8269, 8366:8572, 8527:9207, 9207:9990,
                   10059:10404, 10470:10766, 10760:12137, 12337:14148, 14149:14673, 14747:15887)

all_tissue <- rbind(blood, AA, LV, muscle)

#####################
#correlation with age
#####################

age <- select(pheno, c("SUBJID","AGE","RACE","SEX"))
age$SEX <- factor(age$SEX, levels = c(1,2), labels = c("Male","Female"))
age$RACE <- factor(age$RACE, levels = c(1,2,3,4,98,99), labels = c("Asian", "Black or African American", "White", "American Indian or Alska Native", "Unknown","Unknown"))

risk <- select(pheno, c("BMI", "MHT1D", "MHT2D", "MHSMKSTS", "MHHTN", "SUBJID"))
risk$MHT1D <- factor(risk$MHT1D, levels = c(0,1), labels = c("No","Yes"))
risk$MHT2D <- factor(risk$MHT2D, levels = c(0,1), labels = c("No","Yes"))
risk$MHSMKSTS <- factor(risk$MHSMKSTS, levels = c("No","Yes"), labels = c("No","Yes"))
risk$MHHTN <- factor(risk$MHHTN, levels = c(0,1), labels = c("No","Yes"))

all_result <- read.table("MitImpact_result_all.tsv", header = T, sep = "\t")

all_result <- select(all_result, c("combine", "APOGEE_score","TISSUE"))
all_tissue <- merge(all_tissue, all_result, by = c("combine","TISSUE"))
all_tissue$mutation <- ifelse(all_tissue$APOGEE_score == ".", "synonymous","nonsynonymous")
all_tissue$mutation[str_count(all_tissue$ALT) != str_count(all_tissue$REF)] <- "frameshift"
all_tissue$mutation[!(all_tissue$POS %in% coding.region)] <- "non-coding"

age_correlation_base <- function(name){
  tissue <- subset(all_tissue, TISSUE == name)
  tissue <- subset(tissue, mutation %in% c("nonsynonymous", "synonymous","non-coding"))
  
  summary <- tissue %>% group_by(SUBJID, mutation) %>% summarise(count = n())
  summary <- merge(summary, age, by = "SUBJID")
  summary$mutation <- factor(summary$mutation, levels = c("nonsynonymous", "synonymous", "non-coding"))
  
  model <-lm(count ~ AGE, data=summary)
  output <- summary(model)
  return(output)
}

age_correlation_base_stratified <- function(name){
  tissue <- subset(all_tissue, TISSUE == name)
  tissue <- subset(tissue, mutation %in% c("nonsynonymous", "synonymous","non-coding"))
  
  summary <- tissue %>% group_by(SUBJID, mutation) %>% summarise(count = n())
  summary <- merge(summary, age, by = "SUBJID")
  summary$mutation <- factor(summary$mutation, levels = c("nonsynonymous", "synonymous", "non-coding"))
  
  model1 <-lm(count ~ AGE, data=subset(summary, mutation == "nonsynonymous"))
  model2 <-lm(count ~ AGE, data=subset(summary, mutation == "synonymous"))
  model3 <-lm(count ~ AGE, data=subset(summary, mutation == "non-coding"))
  
  output1 <- summary(model1)
  output2 <- summary(model2)
  output3 <- summary(model3)
  
  return(c(output1, output2, output3))
}

age_correlation_damaging <- function(name){
  tissue <- subset(all_tissue, TISSUE == name)
  tissue <- subset(tissue, mutation == "frameshift" | APOGEE_score > 0.5)
  
  summary <- tissue %>% group_by(SUBJID) %>% summarise(count = n())
  summary <- merge(summary, age, by = "SUBJID")
  summary <- merge(summary, risk, by = "SUBJID")

  model <-lm(count ~ AGE + BMI + MHT1D + MHT2D + MHSMKSTS + MHHTN, data=summary)
  output <- summary(model)
  return(output)
}

lv.age.base.model <- age_correlation_base("Left Ventricle")
lv.age.damaging.model <- age_correlation_damaging("Left Ventricle")

aa.age.base.model <- age_correlation_base("Atrial Appendage")
aa.age.damaging.model <- age_correlation_damaging("Atrial Appendage")

blood.age.base.model <- age_correlation_base("Blood")
blood.age.damaging.model <- age_correlation_damaging("Blood")

muscle.age.base.model <- age_correlation_base("Muscle")
muscle.age.damaging.model <- age_correlation_damaging("Muscle")

lv.stratified <- age_correlation_base_stratified("Left Ventricle")
aa.stratified <- age_correlation_base_stratified("Atrial Appendage")
muscle.stratified <- age_correlation_base_stratified("Muscle")
blood.stratified <- age_correlation_base_stratified("Blood")


#######################
#gene specific variants
#######################

gene_specific_variant <- function(tissue, start, end, total){
  gene <- subset(tissue, POS >= start & POS <= end)
  mutate.pos <- unique(gene$POS)
  test <- matrix(c(16569 - total - (end - start + 1 - length(mutate.pos)),
                   end - start + 1 - length(mutate.pos),
                   total - length(mutate.pos),
                   length(mutate.pos)), byrow = T, ncol = 2)
  output <- fisher.test(test)
  return(output)
}

gene_specific_variant(LV, 3307, 4262, 177)


##############################
#Wilcoxon test between tissues
##############################

#depth
wilcox.test(DP~TISSUE, data = subset(all_tissue, TISSUE == "Left Ventricle" | TISSUE == "Atrial Appendage" ), paried = T)

#VAF
wilcox.test(AF~TISSUE, data = subset(all_tissue, TISSUE == "Left Ventricle" | TISSUE == "Atrial Appendage" ), paried = T)

#number of heteroplasmic variants
wilcox.test(count~TISSUE, paired = TRUE, data = subset(summary_snp_count, TISSUE == "Left Ventricle" | TISSUE == "Atrial Appendage" ))

#APOGEE_score
all_result <- read.table("MitImpact_result_all.tsv", header = T, sep = "\t")
all_result <- select(all_result, c("combine", "APOGEE_score", "TISSUE"))
all_tissue <- merge(all_tissue, all_result, by = c("combine", "TISSUE"))

wilcox.test(APOGEE_score~TISSUE, data = subset(summary_snp_count, TISSUE == "Left Ventricle" | TISSUE == "Atrial Appendage"), paried = T)



