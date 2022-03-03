library(ggplot2)
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

###########################
#violin plot for # variants
###########################

summary_snp_count <- all_tissue %>% group_by(TISSUE, SUBJID) %>% summarise(count = n())

blood_hetero_subj <- subset(summary_snp_count, TISSUE == 'Blood')$SUBJID
blood_hetero_samp <- subset(blood, SUBJID %in% blood_hetero_subj)$SAMPID
LV_hetero_subj    <- subset(summary_snp_count, TISSUE == 'Left Ventricle')$SUBJID
LV_hetero_samp <- subset(LV, SUBJID %in% LV_hetero_subj)$SAMPID
muscle_hetero_subj    <- subset(summary_snp_count, TISSUE == 'Muscle')$SUBJID
muscle_hetero_samp <- subset(muscle, SUBJID %in% muscle_hetero_subj)$SAMPID


blood_no_hetero_subj <- setdiff(LV_hetero_subj, blood_hetero_subj)
blood_no_hetero_samp <- setdiff(LV_hetero_samp, blood_hetero_samp)

muscle_no_hetero_subj <- setdiff(LV_hetero_subj, muscle_hetero_subj)
muscle_no_hetero_samp <- setdiff(LV_hetero_samp, muscle_hetero_samp)

summary_snp_count2 <- data.frame(matrix(ncol = 3, nrow = 88))
colnames(summary_snp_count2) <- c('TISSUE', 'SUBJID', 'count')
summary_snp_count2$TISSUE <- 'Blood'
summary_snp_count2$SUBJID <- blood_no_hetero_subj
summary_snp_count2$count <- 0
summary_snp_count2[nrow(summary_snp_count2) + 1,] = c("Muscle",muscle_no_hetero_subj,0)
summary_snp_count2$count <- as.numeric(summary_snp_count2$count)

summary_snp_count <- rbind(summary_snp_count, summary_snp_count2)

summary_snp_count$TISSUE <- factor(summary_snp_count$TISSUE, levels = c("Blood", "Atrial Appendage","Left Ventricle","Muscle"))


png("variants_num_vs_tissue_voilin_plot.png", width = 1000, height = 600)
ggplot(summary_snp_count, aes(x=TISSUE, y=count, fill=TISSUE)) + 
  geom_violin(color = "grey") +
  geom_boxplot(width=0.1, color = "black", fill = "white") +
  ylab("Number of variants per sample") +
  xlab("Tissue") +
  guides(fill = guide_legend(title="Tissue")) +
  scale_fill_manual(values = cbPalette) +
  ggtitle("Violin Plot of Number of Heteroplasmic Variants Per Sample") + 
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
dev.off()

###################
#violin plot for AF
###################

summary_AF <- all_tissue %>% group_by(TISSUE) %>% summarise(AF)
summary_AF$TISSUE <- factor(summary_AF$TISSUE, levels = c("Blood", "Atrial Appendage","Left Ventricle","Muscle"))

png("mean_AF_voilin_plot.png", width = 800, height = 600)
ggplot(summary_AF, aes(x=TISSUE, y=AF, fill=TISSUE)) + 
  geom_violin(color = "grey") +
  #geom_boxplot(width=0.1, color = "black", fill = "white") +
  scale_fill_manual(values = cbPalette) +
  ylab("VAF") +
  xlab("Tissue") +
  guides(fill = guide_legend(title="Tissue")) +
  ggtitle("Violin Plot of Variant Allele Frequency of Heteroplasmic Variants") + 
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
dev.off()

#########################################
#identify muscle tissue specific variants
#########################################

unique_LV <- subset(LV, !(LV$combine %in% blood$combine))
LV_only <- subset(unique_LV, !(unique_LV$combine %in% unique_AA$combine))
LV_only <- subset(LV_only, !(LV_only$combine %in% unique_muscle$combine))
LV_only$unique <- "Unique in Left Ventricle"

######################################################
#synonymous/nonsynonymous for tissue specific variants
######################################################

summary_mutation <- function(tissue, mutation_file, name){
  tissue <- tissue[order(tissue$POS),]
  result <- read.table(mutation_file, header = T, sep = "\t")
  tissue$mutation <- ifelse(result[12] == ".", "synonymous","nonsynonymous")
  tissue$mutation[str_count(tissue$ALT) != str_count(tissue$REF)] <- "frameshift"
  tissue$mutation[!(tissue$POS %in% coding.region)] <- "non-coding"
  tissue$APOGEE_score <- result$APOGEE_score
  tissue <- tissue[!duplicated(tissue$variant),]
  summary_tissue <- tissue %>% group_by(mutation) %>% summarise(count = n())
  summary_tissue$mutation <- factor(summary_tissue$mutation, levels = c("nonsynonymous", "synonymous","frameshift","non-coding","nonsense"))
  summary_tissue$tissue <- name
  return(summary_tissue)
}

summary.unique.LV <- summary_mutation(LV_only, "unique_lv.vcf", "left ventricle")

png("unique_LV_syn_barplot.png", width = 800, height = 600)
ggplot(data=summary.unique.LV, aes(x=mutation, y=count, fill = mutation)) +
  geom_bar(stat="identity", width=0.5) +
  ylab("Number of variants") +
  xlab("") +
  geom_text(size = 8, label = summary.unique_LV$count) +
  scale_fill_brewer(palette="Set2")+
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 15))
dev.off()

######################
#violin plot for depth
######################

all_tissue$TISSUE <- factor(all_tissue$TISSUE, levels = c("Blood", "Atrial Appendage","Left Ventricle","Muscle"))

png("depth_voilin_plot.png", width = 1000, height = 600)
ggplot(all_tissue, aes(x=TISSUE, y=DP, fill=TISSUE)) + 
  geom_violin(color = "grey") +
  geom_boxplot(width=0.1, color = "black", fill = "white") +
  ylab("Depth") +
  xlab("Tissue") +
  guides(fill = guide_legend(title="Tissue")) +
  scale_fill_manual(values = cbPalette) +
  ggtitle("Violin Plot of Average Depth Per Tissue Type") + 
  theme_bw() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))
dev.off()

############################
#transition vs. transversion
############################

unique.variants <- all_tissue[!duplicated(all_tissue$variant),]
unique.variants$transition <- ifelse(unique.variants$REF %in% c("C","T") & unique.variants$ALT %in% c("C","T") |
                                       unique.variants$REF %in% c("A", "G") & unique.variants$ALT %in% c("A", "G"), 
                                     "transition", 
                                     "transversion")

unique.variants$transition <- ifelse(nchar(unique.variants$REF) != 1 | nchar(unique.variants$ALT) != 1,
                                     "Multi-nucleotide/frameshift", unique.variants$transition)

summary.unique.variants <- unique.variants %>% group_by(transition) %>% summarise(perc = n()/565)
summary.unique.variants$labels <- paste(round(100*summary.unique.variants$perc), "%", sep="")

png("perc_transition.png", width = 600, height = 600)
ggplot(summary.unique.variants, aes(x = "", y = perc, fill = transition)) +
  geom_col() +
  scale_fill_brewer(palette = "Pastel1")+
  geom_text(aes(label = labels),
            position = position_stack(vjust = 0.5),
            size = 15) +
  coord_polar(theta = "y") +
  guides(fill = guide_legend(title = "")) +
  theme_void() +
  theme(legend.text = element_text(size = 17),
        plot.title = element_text(size = 20))
dev.off()

####################################
#frequency of heteroplasmic variants
####################################

all_result <- read.table("MitImpact_result_all.tsv", header = T, sep = "\t")

all_result <- select(all_result, c("combine", "APOGEE_score","TISSUE"))
all_tissue <- merge(all_tissue, all_result, by = c("combine","TISSUE"))
all_tissue$mutation <- ifelse(all_tissue$APOGEE_score == ".", "synonymous","nonsynonymous")
all_tissue$mutation[str_count(all_tissue$ALT) != str_count(all_tissue$REF)] <- "frameshift"
all_tissue$mutation[!(all_tissue$POS %in% coding.region)] <- "non-coding"

summary_frequency <- function(name){
  summary.count <- subset(all_tissue, TISSUE == name) %>% group_by(variant, mutation) %>% summarise(count = n())
  summary.count$group[summary.count$count == 1] <- "Singleton"
  summary.count$group[summary.count$count >= 2 & summary.count$count <= 10] <- "2-10 donors"
  summary.count$group[summary.count$count > 10 & summary.count$count <= 30] <- "10-30 donors"
  summary.count$group[summary.count$count > 30] <- ">30 donors"
  
  summary.freq <- summary.count %>% group_by(group, mutation) %>% summarise(freq = n())
  summary.freq$group <- factor(summary.freq$group, levels = c("Singleton", "2-10 donors","10-30 donors",">30 donors"))
  summary.freq$mutation <- factor(summary.freq$mutation, levels = c("nonsynonymous", "synonymous","frameshift","non-coding"))
  summary.freq$Tissue <- "Left Ventricle"
  return(summary.freq)
}

summary.freq.lv <- summary_frequency("Left Ventricle")
summary.freq.lv2 <- group_by(summary.freq.lv, group) %>% summarise(group, mutation, freq, percentage = freq/sum(freq))

png("variant_freq_lv.png", width = 600, height = 400)
ggplot(data=summary.freq.lv, aes(x=group, y=freq, fill = "Left Ventricle")) +
  geom_bar(stat="identity", width=0.5) +
  xlab("Number of donors") +
  ylab("Frequency") +
  scale_fill_manual(values = cbPalette[3])+
  guides(fill = guide_legend(title = "")) +
  #ggtitle("Left Ventricle") +
  theme_classic() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))
dev.off()

png("variant_freq_lv_stratified.png", width = 600, height = 400)
ggplot(data=summary.freq.lv2, aes(x=group, y=percentage, fill = mutation)) +
  geom_bar(stat="identity", width=0.5) +
  xlab("Number of donors") +
  ylab("Frequency") +
  scale_fill_brewer(palette = "Set2")+
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = "")) +
  ggtitle("Left Ventricle") +
  theme_classic() +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))
dev.off()

##################
#percentage of VAF
##################

vaf_percentage <- function(tissue, total){
  tissue$VAF[tissue$AF <= 0.25] <- "<= 25%"
  tissue$VAF[tissue$AF <= 0.5 & tissue$AF > 0.25] <- "25% - 50%"
  tissue$VAF[tissue$AF <= 0.75 & tissue$AF > 0.5] <- "50% - 75%"
  tissue$VAF[tissue$AF > 0.75] <- "> 75%"
  summary.tissue <- tissue %>% group_by(VAF) %>% summarise(count = n())
  summary.tissue$VAF <- factor(summary.tissue$VAF, levels = c("<= 25%", "25% - 50%","50% - 75%","> 75%"))
  summary.tissue$freq <- summary.tissue$count/total
  return(summary.tissue)
}

LV_unique_variant <- LV[!duplicated(LV$variant),]
summary.vaf.lv <- vaf_percentage(LV_unique_variant, 202)

png("vaf_lv.png", width = 600, height = 400)
ggplot(data=summary.vaf.lv, aes(x=VAF, y=freq, fill = "Left Ventricle")) +
  geom_bar(stat="identity", position = 'dodge', width = 0.5) +
  scale_fill_manual(values = cbPalette[3])+
  xlab("VAF") +
  ylab("Proportion of Variants, %") +
  guides(fill = guide_legend(title = "")) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic()+
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 17))
dev.off()



