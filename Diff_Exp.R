library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(FSA)  # Dunn test 
library(GmAMisc) # KwPlot 
library(rstatix)


##### Import Data ####
setwd("/Users/chensisi/Documents/RNAseq/")

# Transcript data import 
transcript <- read.table("TCGA_GTEx_RSK4_tpm.txt", header = TRUE,  row.names=1)
transcript1 <- as.data.frame(t(transcript)) # Transpose the data frame 
transcript1[1:3,]

# Metadata import 
metadata <- read.table("TCGA_GTEX_category.txt",  header = TRUE, row.names = 1, fill = TRUE,sep = "\t")

metadata <- metadata %>% 
  separate(TCGA_GTEX_main_category, c("Project", "Sample_type"), sep = "\\s", extra = "merge")

rownames(metadata)<-gsub(rownames(metadata), pattern="-", replace=".")
metadata[1:3,]


# Survival data import 
survival <- read.table("TCGA_survival_data", sep="\t", header=T, row.names=1) 
rownames(survival)<-gsub(rownames(survival), pattern="-", replace=".")
survival[1:3,]



# Merge two dataframes by matching sample IDs 
merge1 <- merge(transcript1, metadata,  by = 0)
rownames(merge1) <- merge1[,1]
merge1 <- merge1[,-1]
merge2 <- merge(merge1, survival, by = 0)

# Convert merged dataframe into tibble for downstream data manipulation 
tibble1 <- as_tibble(merge1)
unique(merge1$Project)
unique(merge1$Sample_type)

TCGA_by_sample <- tibble1 %>% 
  filter(Project == "TCGA") %>% 
  group_by(Sample_type) 

TCGA_samples <- unique(TCGA_by_sample$Sample_type)


GTEX_by_sample <- tibble1 %>% 
  filter(Project == "GTEX") %>% 
  group_by(Sample_type)

GTEX_samples <- unique(GTEX_by_sample$Sample_type)

###### Normality test ####
  # By project: TCGA vs GTEX 
qplot(sample = ENST00000620340.4, data = merge1 , color= Project) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression", 
       y = "log2(x+0.0001) transfromed TPM", x = "Threoretical") 

  # TCGA 
qplot(sample = ENST00000620340.4, data = TCGA_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression in TCGA samples", 
       y = "log2(x+0.0001) transfromed TPM", x = "Threoretical")

qplot(sample = ENST00000262752.4, data = TCGA_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 2 expression in TCGA samples", 
     y = "log2(x+0.0001) transfromed TPM", x = "Threoretical")

  # GTEX
qplot(sample = ENST00000620340.4, data = GTEX_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 1 expression in GTEX samples", 
       y = "log2(x+0.0001) transfromed TPM", x = "Threoretical")

qplot(sample = ENST00000262752.4, data = GTEX_by_sample, color = Sample_type) + 
  labs(title = "Q-Q plot for RSK4 isoform 2 expression in GTEX samples", 
       y = "log2(x+0.0001) transfromed TPM", x = "Threoretical")

# Shapiro test 

  #TCGA 
tcga_normality <- data.frame()

for (x in TCGA_samples){
  a <- TCGA_by_sample %>%
    select(ENST00000620340.4,ENST00000262752.4,Sample_type) %>%
    filter(Sample_type == x)
  shapiro_p_value <- shapiro.test(a$ENST00000620340.4)[[2]]
  result <- data.frame(x, shapiro_p_value)
  tcga_normality <- rbind(normality, result)
}

tcga_normality %>% filter(shapiro_p_value >= 0.05)

  # GTEX 
gtex_normality <- data.frame()

for (x in GTEX_samples){
  a <- GTEX_by_sample %>%
    select(ENST00000620340.4,ENST00000262752.4,Sample_type) %>%
    filter(Sample_type == x)
  shapiro_p_value <- shapiro.test(a$ENST00000620340.4)[[2]]
  result <- data.frame(x, shapiro_p_value)
  gtex_normality <- rbind(gtex_normality, result)
}

gtex_normality %>% filter(shapiro_p_value >= 0.05)

#### Histogram and density plot ####

  # TCGA 
ggplot(TCGA_by_sample, aes(x=ENST00000620340.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 1 expression in TCGA samples")

ggplot(TCGA_by_sample, aes(x=ENST00000262752.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 2 expression in TCGA samples")


  # GTEX
ggplot(GTEX_by_sample, aes(x=ENST00000620340.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 1 expression in GTEX samples")

ggplot(GTEX_by_sample, aes(x=ENST00000262752.4, fill = Sample_type)) + 
  geom_density(alpha=.3) + 
  labs(title = "Density plot for RSK4 isoform 2 expression in GTEX samples")

# Box plot
# Isoform I: ENST00000620340.4
# Isoform II: ENST00000262752.4
# TCGA samples 
ggplot(TCGA_by_sample, aes(x=Sample_type, y = ENST00000620340.4)) + 
  geom_boxplot() + 
  labs(title = "Comparing RSK4 isoform1 expression between different tumour samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1))

ggplot(TCGA_by_sample, aes(x=Sample_type, y = ENST00000262752.4)) + 
  geom_boxplot() + 
  labs(title = "Comparing RSK4 isoform2 expression between different tumour samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1))

# GTEX samples
ggplot(GTEX_by_sample, aes(x=Sample_type, y = ENST00000620340.4)) + 
  geom_boxplot() + 
  labs(title = "Comparing RSK4 isoform1 expression between normal tissue samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) + 
  theme(plot.title = element_text(hjust = 0.5)) #标题居中

ggplot(GTEX_by_sample, aes(x=Sample_type, y = ENST00000262752.4)) + 
  geom_boxplot() + 
  labs(title = "Comparing RSK4 isoform2 expression between normal tissue samples") + 
  theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust=1)) +
  theme(plot.title = element_text(hjust = 0.5))

##### Test for difference -- Kruskal-Wallis Test ####
# A collection of data samples are independent if they come from unrelated populations and the samples do not affect each other. 
# Using the Kruskal-Wallis Test, we can decide whether the population distributions are identical without assuming them to follow the normal distribution.
  ## TCGA 
tcga_iso1<- kruskal.test(ENST00000620340.4 ~ Sample_type, data = TCGA_by_sample)
tcga_iso1

tcga_iso2<- kruskal.test(ENST00000262752.4 ~ Sample_type, data = TCGA_by_sample)
tcga_iso2

  # Multiple Pairwise comparisons -- Dunn test 
tcga_iso1_pwc <- dunnTest(ENST00000620340.4 ~ as.factor(Sample_type), data = TCGA_by_sample, 
                          method="bonferroni")
df1 <- tcga_iso1_pwc$res

tcga_iso2_pwc <- dunnTest(ENST00000262752.4 ~ as.factor(Sample_type), data = TCGA_by_sample, 
                          method="bonferroni")
df2 <- tcga_iso2_pwc$res

df2 %>%
  filter(P.adj < 0.05 ) %>%
  arrange(P.adj)

write.csv(df2, "TCGA_by_sample_difference_iso2.csv")



  ## GTEX 
gtex_iso1<- kruskal.test(ENST00000620340.4 ~ Sample_type, data = GTEX_by_sample)
gtex_iso1

gtex_iso2<- kruskal.test(ENST00000262752.4 ~ Sample_type, data = GTEX_by_sample)
gtex_iso2

# Multiple Pairwise comparisons -- Dunn test 
gtex_iso1_pwc <- dunnTest(ENST00000620340.4 ~ as.factor(Sample_type), data = GTEX_by_sample, 
                          method="bonferroni")
df3 <- gtex_iso1_pwc$res

gtex_iso2_pwc <- dunnTest(ENST00000262752.4 ~ as.factor(Sample_type), data = GTEX_by_sample, 
                          method="bonferroni")
df4 <- gtex_iso2_pwc$res

df4 %>%
  filter(P.adj < 0.05 ) %>%
  arrange(P.adj)

write.csv(df3, "GTEX_by_sample_difference_iso1.csv")


#### One-way ANOVA #####
# Compute the analysis of variance
res.aov <- aov(ENST00000620340.4 ~ as.factor(Sample_type), data = TCGA_by_sample)
# Summary of the analysis
summary(res.aov)

# Tukey multiple pairwise-comparisons
TukeyHSD(res.aov)
# 1. Homogeneity of variances -- check the homogeneity of variance assumption 
plot(res.aov,1)

# 2. Normality
plot(res.aov, 2)
# Extract the residuals
aov_residuals <- residuals(object = res.aov )
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )




