library(data.table)
library(tidyverse)

############
#input files
############

load_data <- function(filename){
  tissue <- read.table(filename, header = F, sep = "\t")
  names(tissue) = c("SAMPID","CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","GT","AD","AF","DP","F1R2","F2R1","SB")
  return(tissue)
}

data_qc <- function(tissue, drop, region, depth, name){
  tissue$SUBJID <- sapply(tissue$SAMPID %>% strsplit("-"), function(x) paste0(x[1:2], collapse = "-"))
  tissue$TISSUE <- name
  tissue <- subset(tissue, !(POS %in% c(301, 302, 310,316, 513,3107)))
  tissue <- subset(tissue, !((POS >= 66 & POS <= 71) | (POS >= 12418 & POS <= 12425) | (POS >= 16182 & POS <= 16194)))
  tissue <- subset(tissue, !((POS == 295 & REF == 'C' & ALT == 'T') | (POS == 295 & REF == 'CC' & ALT == 'TC') | (POS == 13710 & REF == 'A' & ALT == 'G') | (POS == 2617 & REF == 'A' & ALT == 'G') | (POS == 2617 & REF == 'AT' & ALT == 'GT')))
  tissue <- subset(tissue, !(POS %in% region))
  tissue <- subset(tissue, DP >= depth)
  tissue <- subset(tissue, !(SUBJID %in% drop))
  tissue <- subset(tissue, AF < 0.97 & AF > 0.03)
  return(tissue)
}

data_combine <- function(tissue){
  tissue$variant <- as.character(interaction(tissue$POS, tissue$REF, tissue$ALT))
  tissue$combine = as.character(interaction(tissue$POS, tissue$REF, tissue$ALT, tissue$SUBJID))
  return(tissue)
}

blood  <- load_data("./all_blood_no_filter.txt")
AA     <- load_data("./all_AA_no_filter.txt")
LV     <- load_data("./all_LV_no_filter.txt")
muscle <- load_data("./all_muscle_no_filter.txt")

tRNA <- c((577:647),(1602:1670),(3230:3304),(4263:4469),(5512:5891),(7446:7585),(8295:8364),(9991:10058),(10405:10469),(12138:12336),(14674:14742),(15888:16023))
drop <- c("GTEX-111FC", "GTEX-WRHU")

blood  <- data_qc(blood, drop, tRNA, 500, "Blood")
AA     <- data_qc(AA, drop, tRNA, 500, "Atrial Appendage")
LV     <- data_qc(LV, drop, tRNA, 500, "Left Ventricle")
muscle <- data_qc(muscle, drop, tRNA, 500, "Muscle")

blood  <- data_combine(blood)
AA     <- data_combine(AA)
LV     <- data_combine(LV)
muscle <- data_combine(muscle)

#############
#output files
#############

save(blood, file = "./blood.RData")
save(AA, file = "./aa.RData")
save(LV, file = "./lv.RData")
save(muscle, file = "./muscle.RData")
