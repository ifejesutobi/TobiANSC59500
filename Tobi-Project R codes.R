##############################################################
# title: "Alpha  and Beta diversity in R - qiime2 output"
# author: "Ogunribido Tobi -ANSC595"
# Working dir: C:\Users\17657\OneDrive\Documents\Project Data
# Insructor: Tim Johnson
# date: "March 26, 2021"
##############################################################

##Files Used (All created from qiime2 moving pictures tutorial)

# core-metrics-results/evenness_vector.qza (alpha diversity)
# core-metrics-results/faith_pd_vector.qza (alpha diversity)
# core-metrics-results/observed_features_vector.qza (alpha diversity)
# core-metrics-results/shannon_vector.qza (alpha diversity)
# TobiMetadata.txt (metadata)

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
# In the code, the text before = is what the file will be called in R. 
# Make this short but unique as this is how you will tell R to use this 
# file in later commands.

# header: tells R that the first row is column names, not data
# row.names: tells R that the first column is row names, not data
# sep: tells R that the data are tab-delimited. 
# If you had a comma-delimited file, you would us sep=","

# Load data

Tobimeta<-read_q2metadata("TobiMetadata.txt")

evenness = read_qza("evenness_vector.qza")
evenness <- evenness$data %>% rownames_to_column("SampleID") #this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged


simpson = read_qza("simpson_vector.qza")
simpson<-simpson$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\


## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(Tobimeta)
#observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
#observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)

###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = evenness
alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, simpson, by.x ="SampleID", by.y = "SampleID")


meta = merge(Tobimeta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) = meta$SampleID
#meta = meta[,-1]
str(Tobimeta)

#Alpha-diversity
# Alpha-diversity is within sample diversity. It is how many 
# different species (OTUs) are in each sample (richness) and how 
# evenly they are distributed (evenness), which together are diversity. 
# Each sample has one value for each metric.


##Explore alpha metrics
# Now we will start to look at our data. We will first start with 
# alpha-diversity and richness. 
#
# You want the data to be roughly normal so that you can run ANOVA 
# or t-tests. If it is not normally distributed, you will need to 
# consider if you should normalize the data or usenon-parametric 
# tests such as Kruskal-Wallis.

# Here, we see that none of the data are normally distributed, 
# with the exception of "Faith" and "Observed Features".


#Plots
hist(as.numeric(meta$evenness), main = "Evenness", xlab ="", breaks=10)
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)
hist(meta$simpson, main="Simpson diversity", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta$evenness, title= "Evenness")
ggqqplot(meta$shannon, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$observed_features, title = "Observed Features")
ggqqplot(meta$simpson, title = "Simpson")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$observed_features)
shapiro.test(meta$simpson)


# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant, the distribution is non-normal.

# We see that, as expected from the graphs, shannon and evenness 
# are normally distributed.


#Overall, for alpha-diversity:

# ANOVA, t-test, or general linear models with the normal distribution 
# are used when the data is roughly normal. Transforming the data to 
# achieve a normal distribution could also be completed.
#
# Kruskal-Wallis, Wilcoxon rank sum test, or general linear models 
# with another distribution are used when the data is not normal or if 
# the n is low, like less than 30.

# Our main variables of interest are

# body site: gut, tongue, right palm, left palm
# subject: 1 and 2
# month-year: 10-2008, 1-2009, 2-2009, 3-2009, 4-2009

## Categorical variables
# Now that we know which tests can be used, let's run them. 

## Normally distributed metrics

# Since it's the closest to normalcy, we will use **Evenness** as an 
#example. First, we will test body site, which is a categorical variable 
# with more than 2 levels. Thus, we run ANOVA. If age were only two 
# levels, we could run a t-test

# Does body site impact the Evenness of the microbiota?

#Run the ANOVA and save it as an object

colnames(meta)[3] <- "delivery_mode"
aov.evenness.delivery_mode = aov(pielou_evenness ~ delivery_mode, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.evenness.delivery_mode)

#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(aov.evenness.delivery_mode)

# We clearly see that the evenness between hands and gut are different. 
# When we plot the data, we see that evenness decreases in the gut 
# compared to palms.

levels(meta$delivery_mode)
#Re-order the groups because the default is 1yr-2w-8w
meta$delivery_mode.ord = factor(meta$delivery_mode, c("Post-labor CS", "Pre-labor CS", "Vaginal"))
levels(meta$delivery_mode.ord)

#Plot
boxplot(pielou_evenness ~ delivery_mode.ord, data=meta, ylab="Peilou evenness")

evenness <- ggplot(meta, aes(delivery_mode.ord, pielou_evenness)) + 
  geom_boxplot(aes(color = delivery_mode.ord)) + 
  ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/evenness.png", evenness, height = 3, width = 3)

# Now, the above graph is kind of not correct. Our test and our graphic do not exactly match. ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.

evenness_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(delivery_mode.ord) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

evenness_se <- ggplot(evenness_summary, aes(delivery_mode.ord, mean_evenness, fill = delivery_mode.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 

ggsave("output/evenness_se.png", evenness_se, height = 3, width = 3)


## **Non-normally distributed metrics**

# We will use **Faith's phylogenetic diversity** here. Since body site 
# is categorical, we use Kruskal-Wallis (non-parametric equivalent of 
# ANOVA). If we have only two levels, we would run Wilcoxon rank sum 
# test (non-parametric equivalent of t-test)

kruskal.test(faith_pd ~ delivery_mode.ord, data=meta)

# We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

pairwise.wilcox.test(meta$faith_pd, meta$delivery_mode.ord, p.adjust.method="BH")

# Like evenness, we see that pd also increases with age.

#Plot
boxplot(faith_pd ~ delivery_mode.ord, data=meta, ylab="Faith phylogenetic diversity")

# or with ggplot2

faith_pd <- ggplot(meta, aes(delivery_mode.ord, faith_pd)) + 
  geom_boxplot(aes(color = delivery_mode.ord)) + 
  #ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/pd.png", faith_pd, height = 3, width = 3)

##Continuous variables

shannon <- ggplot(meta, aes(delivery_mode.ord, shannon)) + 
  geom_boxplot(aes(color = delivery_mode.ord)) + 
  #ylim(c(0.5,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Diversity", x = "") 
ggsave("output/pd.png", Shannon, height = 3, width = 3)

# For continuous variables, we use general linear models, specifying 
# the distribution that best fits our data.

# **Normally distributed metrics**

# Since days.since.experiment.start is a continuous variable, we run a 
# general linear model. We will again use evenness as our roughly normal 
# metric. The default of `glm` and `lm` is the normal distribution so we 
# don't have to specify anything.

# Does days.since.experiment.start impact evenness of the microbiota?
colnames(meta)[9] <- "days.since.experiment.start"
glm.evenness.time = glm(pielou_evenness ~ delivery_mode, data=meta)
summary(glm.evenness.time)

#The output let's us know that the intercept of our model is significantly different from 0 but our slope (*e.g.* our variable of interest) is not. This makes sense when we look at the data.

plot(pielou_evenness ~ delivery_mode, data=meta)
#Add the glm best fit line
plot(pielou_evenness ~ days.since.experiment.start, data=meta) + abline(glm.evenness.time)

# **Non-normally distributed metrics**

# We will again use a *general linear model* for our non-normally 
# distributed metric Faith_pd. However, this time, we change the 
# distribution from normal to something that fits the data better. 

# But which distribution should we choose? In statistics, there is no 
# one "best" model. There are only good and better models. We will use 
# the plot() function to compare two models and pick the better one.

# First, the Gaussian (normal) distribution, which we already know is a bad fit.

gaussian.faith.time = glm(faith_pd ~ delivery_mode, data=meta, family="gaussian")
plot(gaussian.faith.time, which=c(1,2))

# Quasipoisson (log) distribution
qp.faith.time = glm(faith_pd ~ delivery_mode, data=meta, family="quasipoisson")
plot(qp.faith.time, which=c(1,2))

# What we're looking for is no pattern in the Residuals vs. Fitted graph 
# ("stars in the sky"), which shows that we picked a good distribution 
# family to fit our data. We also want our residuals to be normally 
# distributed, which is shown by most/all of the points falling on the 
# line in the Normal Q-Q plot.

# While it's still not perfect, the quasipoisson fits much better. 
# In the residuals vs fitted graph, the y axis is from -2 to 4  whereas 
# the axis with gaussian was from -5 to 10. So, we will use quasipoisson 
# and see that ADG does not to correlate to Chao richness.
summary(qp.faith.time)

# Plotting this we see that, indeed, there is a trend toward correlation between Faith_pd and days.since.experiment.start.

#Plot
plot(log(faith_pd) ~ delivery_mode, data=meta, ylab="ln(Faith Phylo. Diversity)")
plot(log(faith_pd) ~ delivery_mode, data=meta, ylab="ln(Faith Phylo. Diversity)") + abline(qp.faith.time)




##############################################################
# title: "Beta diversity in R - qiime2 output"
# author: "Ogunribido Tobi -ANSC595"
# Working dir: C:\Users\17657\OneDrive\Documents\Project Data
# Insructor: Tim Johnson
# date: "March 26, 2021"
##############################################################

#Load the packages

library(qiime2R)

library(tidyverse)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# To get these files you need to scp them from the cluster:
#
# first on  your laptop cd to the directory where you want to save them.
# Then use this code for our example dataset today:
# mkdir core-metrics-results/
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/* .
# scp john2185@bell.rcac.purdue.edu:/depot/microbiome/data/2021_ANSC595/john2185/qiime/moving_pictures_pipeline/core-metrics-results/* core-metrics-results/.
##############################################


setwd("C:\Users\17657\OneDrive\Documents\Downloads")
#setwd("C:\Users\17657\OneDrive\Documents\Downloads")
list.files()

if(!dir.exists("output"))
  dir.create("output")

metadata<-read_q2metadata("TobiMetadata.txt")
#metadata2 <- read.delim("TobiMetadata.txt", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#metadata2 <- metadata2[-1,]
str(metadata)
levels(metadata$`delivery_mode`)



row.names(metadata) <- metadata[ ,1]
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("bray_curtis_pcoa_results.qza")

body_colors <- c("Black", "Blue", "Green")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

my_column <- "delivery_mode"
#my_column <- "DietTreatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-basic_", my_column,".pdf"), height=2, width=3, device="pdf") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "delivery_mode"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= delivery_mode)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=corn_colors, name = my_column)
ggsave(paste0("output/BC-ellipse_", my_column,"-subject.pdf", main = "Bray-Curtis-PcoA"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

##SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Delivery Mode")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= delivery_mode), size = 2) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "delivery_mode")
ggsave(paste0("output/Wuni-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches


#Jaccard_pcoa_results

Jacc_PCoA<-read_qza("Jaccard_pcoa_results.qza")

Jacc_meta <- Jacc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Jacc_meta,mean)

ggplot(Jacc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jacc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jacc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Delivery Mode")
ggsave(paste0("output/Jacc-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Jacc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= delivery_mode), size = 2) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jacc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jacc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "delivery_mode")
ggsave(paste0("output/Jacc-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches



#Unweighted Unifrac PcoA Results

Uwuni_PCoA<-read_qza("unweighted_unifrac_pcoa_results.qza")

Uwuni_meta <- Uwuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Uwuni_meta,mean)

ggplot(Uwuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "Delivery Mode")
ggsave(paste0("output/Uwuni-ellipse_", my_column,".pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Uwuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= delivery_mode), size = 2) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Uwuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Uwuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=body_colors, name = "delivery_mode")
ggsave(paste0("output/Uwuni-ellipse_", my_column,"-subject.pdf"), height=3, width=4.5, device="pdf") # save a PDF 3 inches by 4 inches




####$$$ PART 2: DESeq2 codes
#### Get Ready!

## This script goes through making a taxa bar plot and then running
## DESeq2 to find differentially abundant ASVs

##################Working dir:  C:\Users\17657\OneDrive\Documents\Downloads

# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)

##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#TobiMetadata.tsv
#core-metrics-results/rarefied_table.qza
#Tobi-rooted-tree.qza
#Tobi-taxonomy.qza
##############################################


if(!dir.exists("output"))
  dir.create("output")

Tobimeta<-read_q2metadata("TobiMetadata.txt")
str(Tobimeta)
colnames(Tobimeta)[3] = "family_ID"
levels(Tobimeta$delivery_mode)

row.names(Tobimeta) <- Tobimeta[ ,1]
#Tobimeta <- Tobimeta[,-1]
row.names(Tobimeta)

##Qiime2r method of reading in the taxonomy files
taxonomy<-read_qza("newTobi-taxonomy.qza")
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)

#All this is OK except that in future use of the taxonomy table, 
#these ASVs will be ignored because they are not classified. Why 
#are ASVs not classified? Its because there is not a close enough 
#match in the database. Just because there is not a good match in 
#the database does not mean they don’t exist, so I wanted to make 
#sure this data was not lost. So in my new code, from lines 200 – 224 
#I make it so that ASVs that are unclassified at any level are 
#classified as the lowest taxonomic level for which there is a 
#classification.
#Next, all these `NA` classifications with the last level that was 
#classified

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}


#################################################################
##Taxa barplot
#################################################################

physeq <- qza_to_phyloseq(
  features="rarefied_table.qza",
  tree="rooted-tree.qza",
  "newTobi-taxonomy.qza",
  metadata = "TobiMetadata.txt"
)


#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = F)

tax.clean = tax.clean[row.names(tax.clean) %in% rownames(physeq_otu_table),]
metadata.filtered = Tobimeta[row.names(Tobimeta) %in% colnames(physeq_otu_table),]

#Assign as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(physeq_otu_table), taxa_are_rows=TRUE)

#our edited and formatted taxonomy table from the top of this script
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata.filtered)

#We then merge these into an object of class phyloseq.

physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)



# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black",
  "purple", "yellow", "green", "brown", "tan"
)

#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "delivery_mode"

rm(taxa.summary)

abund_filter <- 0.01 #means ASV should have at least, an abundance of 1% to be shown in plot. you can change this depending on your expectation.
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
   group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance), .groups="rowwise") 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.average <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.average <- as.data.frame(physeq.taxa.average)
  colnames(physeq.taxa.average)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.average)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  #physeq_meta_filtered$delivery_mode.ord = factor(physeq_meta_filtered$delivery_mode, c("left palm", "right palm", "gut", "tongue"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = get(my_column), y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~LitterTreatment) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    #theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 delivery group")) 
  ggsave(paste0("output/", ml, "BarPlot_", my_column, ".png"), height = 5)
}

#################################################################
###Differential Abundance with DESeq2
#################################################################


#Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html

#First load DESeq2.
#If you need help  with DESeq2 install, see this website
https://bioconductor.org/packages/release/bioc/html/DESeq2.html

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")

#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1

#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = FALSE)

OTU.clean2 <- physeq_otu_table + 1


#Now make the phyloseq object:


OTU.physeq = otu_table(as.matrix(OTU.clean2), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(tax.clean))
meta.physeq = sample_data(metadata.filtered)


#We then merge these into an object of class phyloseq.


physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)


#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.

diagdds = phyloseq_to_deseq2(physeq_deseq, ~ delivery_mode)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
#the test type of "Wald" tests for significance of coefficients in a Negative Binomial GLM. This is generally a pretty good assumption for sequencing experiments. This was designed with RNA-seq in mind, but also pretty good for 16S sequencing.


###Investigate test results table

#The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the diagdds object (see above). I then order by the adjusted p-value, removing the entries with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.

#Contrast: this argument specifies what comparison to extract from the object to build a results table. There are exactly three elements:

#  1. the name of a factor in the design formula, 
#  2. the name of the numerator level for the fold change, and 
#  3. the name of the denominator level for the fold change (simplest case)

alpha = 0.05
my_contrast = c("delivery_mode", "Post-labor CS", "Pre-labor CS") 
#my_contrast = c("delivery_mode", "Post-labor CS", "Vaginal") 
#my_contrast = c("delivery_mode", "Pre-labor CS", "Vaginal") 

res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)

sigtab = res[which(res$padj < alpha), ]
sigtab_test <- as(sigtab, "data.frame")
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
head(sigtab)


###Volcano Plot

with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
#x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
DESeq_fig = ggplot(sigtab, aes(x=Genus, y = log2FoldChange, color=Phylum)) + 
  geom_point(size=3) + 
  ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
  scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
  #ylim(0,8) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 5.5, width = 9.9)


#######################################################
####LEFSE ANALYSIS
#######################################################
library(qiime2R)
library(phyloseq)
library(tidyverse)

list.files()

genus_freq_meta <- read.csv("collapse.frequency.table.txt", sep = "\t", header = F)

genus_freq_meta2 <- separate(data = genus_freq_meta, col = "V1", into = c("kingdom", "phylum", "class", "order", "family", "genus"), sep = "[|]")
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "k__", replacement = "", fixed = TRUE)
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "p__", replacement = "", fixed = TRUE)
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "c__", replacement = "", fixed = TRUE)
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "o__", replacement = "", fixed = TRUE)
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "f__", replacement = "", fixed = TRUE)
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "g__", replacement = "", fixed = TRUE)
genus_freq_meta2[]  <- lapply(genus_freq_meta2, gsub, pattern = "__", replacement = "", fixed = TRUE)

n_metadata_rows <- 3

i = 4
for (i in (n_metadata_rows+1):nrow(genus_freq_meta2)){
  if (genus_freq_meta2[i,2] == ""){
    kingdom <- paste("unclassified_", genus_freq_meta2[i,1], sep = "")
    genus_freq_meta2[i, 2:6] <- kingdom
  } else if (genus_freq_meta2[i,3] == ""){
    phylum <- paste("unclassified_", genus_freq_meta2[i,2], sep = "")
    genus_freq_meta2[i, 3:6] <- phylum
  } else if (genus_freq_meta2[i,4] == ""){
    class <- paste("unclassified_", genus_freq_meta2[i,3], sep = "")
    genus_freq_meta2[i, 4:6] <- class
  } else if (genus_freq_meta2[i,5] == ""){
    order <- paste("unclassified_", genus_freq_meta2[i,4], sep = "")
    genus_freq_meta2[i, 5:6] <- order
  } else if (genus_freq_meta2[i,6] == ""){
    family <- paste("unclassified_", genus_freq_meta2[i,5], sep = "")
    genus_freq_meta2[i, 6:6] <- family
  } 
}
genus_freq_meta3 <- unite(genus_freq_meta2, col = "taxon", c("kingdom", "phylum", "class", "order", "family", "genus"), sep = "|")
genus_freq_meta3[]  <- lapply(genus_freq_meta3, gsub, pattern = "|NA|NA|NA|NA|NA", replacement = "", fixed = TRUE)

write_tsv(genus_freq_meta3, "collapse.frequency.lefse.tsv", col_names = F)




####################################################################
##Co-Occurrence Analysis
######################################################################
###
 

#make sure you have these libraries
library(Hmisc)
library(plyr)
library(reshape2)
library(qiime2R)
#library(igraph)
#library(fdrtool)
# this is the demo data, take a look at it. YOU MUST PROPERLY DIRECT THE FILE PATH

setwd(dir = "~/Desktop/ANSC595/moving-pictures/")

ASVs <- read_qza("newTobi-table.qza")
ASV_table <- as.data.frame(ASVs$data)

#####################################################################
##I added in this chuck. What does it do and why?

ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)]
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
######################################################################

dataset <- as.data.frame(t(ASV_table))


# we are going to create a network per treatment
head(dataset[,1:10])

metadata<-read_q2metadata("TobiMetadata.txt")
str(metadata)
#colnames(metadata)[3] = "family_ID"

dataset <- merge(metadata, dataset, by.x = "SampleID", by.y = 0)
treatments<-as.vector(unique(dataset$delivery_mode))
datasetn<-dataset
datasetn[datasetn==0]<-NA

#i = 1

summary(metadata$delivery_mode)

my_column <- "delivery_mode"
n1 <- 4
n2 <- 4
n3 <- 4
n4 <- 4
n5 <- 4

num_metadata_columns <- 24 # Here is the number of columns in the dataset dataframe before the ASVs

q_cutoff <- 0.10

final_results<-data.frame()

for(i in 1:length(treatments)){
  #subset the data for a particular treatment YOU MUST ENTER THE HEADER OF THE COLUMN THAT HAS THE DIFFERENT TREATMENTS IN THIS CASE “Foaming_Status”
  print(paste("reading ",treatments[i],sep=""))
  temp<-subset(dataset, get(my_column)==treatments[i])
  tempn<-subset(datasetn, get(my_column)==treatments[i])
  print(paste("finished reading ",treatments[i],sep=""))
  # making an object that has all the results in it (both rho and P values)
  results<-rcorr(as.matrix(temp[,-c(1:num_metadata_columns)]),type="spearman") ## use the "-c" parameter to remove metadata columns
  resultsn<-rcorr(as.matrix(tempn[,-c(1:num_metadata_columns)]),type="spearman")
  
  #make two seperate objects for p-value and correlation coefficients
  rhos<-results$r
  ps<-results$P
  ns<-resultsn$n
  # going to melt these objects to 'long form' where the first two columns make up the pairs of OTUs, I am also removing NA's as they are self-comparisons, not enough data, other bad stuff
  ps_melt<-na.omit(melt(ps))
  #creating a qvalue based on FDR
  ps_melt$qval<-p.adjust(ps_melt$value, method = "BH")
  #making column names more relevant
  
  names(ps_melt)[3]<-"pval"
  # if you are of the opinion that it is a good idea to subset your network based on adjusted P-values (qval in this case), you can then subset here
  ps_sub<-subset(ps_melt, qval < q_cutoff)
  
  # now melting the rhos, note the similarity between ps_melt and rhos_melt
  rhos_melt<-na.omit(melt(rhos))
  names(rhos_melt)[3]<-"rho"
  
  # now melting the ns
  ns_melt<-(melt(ns))
  names(ns_melt)[3]<-"n"
  
  #merging together and remove negative rhos
  merged<-merge(ps_sub,rhos_melt,by=c("Var1","Var2"))
  if (treatments[i]==treatments[1]) {
    merged<-merge(merged,subset(ns_melt, n > n1),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[2]) {
    merged<-merge(merged,subset(ns_melt, n > n2),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[3]) {
    merged<-merge(merged,subset(ns_melt, n > n3),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[4]) {
    merged<-merge(merged,subset(ns_melt, n > n4),by=c("Var1","Var2"))
  }   else if (treatments[i]==treatments[5]) {
    merged<-merge(merged,subset(ns_melt, n > n5),by=c("Var1","Var2"))
  }   else
    print("Somethings wrong with your treatment designations. Please Check!!")
  
  if (nrow(merged) > 0) {
    merged$trt<-treatments[i]
    final_results<-rbind(final_results, merged)
  }   else {
    print("no correlations for this variable")
  }
  
  print(paste("finished ",treatments[i],sep=""))
}

strong_results<-subset(final_results, rho >= 0.7)


###############################################################
# If you want to see the correlation scatterplot of 
# two significant ASVs
###############################################################


gut_ASVs<-subset(dataset, get(my_column)==treatments[1])
gut_ASVs<-subset(dataset, get(my_column)=="Post-labor CS")

colnames(gut_ASVs[1:10])
head(final_results)
ggplot(gut_ASVs, aes(x = ASV99, y = ASV158)) +
  geom_point()



###############################################################
# If you want to see the the taxonomic assignment of  
# these significant ASVs
###############################################################

taxonomy<-read_qza("newTobi-taxonomy.qza")
head(taxonomy$data)

tax.clean<-parse_taxonomy(taxonomy$data)
head(tax.clean)

#All this is OK except that in future use of the taxonomy table, 
#these ASVs will be ignored because they are not classified. Why 
#are ASVs not classified? Its because there is not a close enough 
#match in the database. Just because there is not a good match in 
#the database does not mean they don’t exist, so I wanted to make 
#sure this data was not lost. So in my new code, from lines 200 – 224 
#I make it so that ASVs that are unclassified at any level are 
#classified as the lowest taxonomic level for which there is a 
#classification.
#Next, all these `NA` classifications with the last level that was 
#classified

tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("unclassified_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("unclassified_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("unclassified_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("unclassified_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("unclassified_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("unclassified_",tax.clean$Genus[i], sep = "_")
  }
}

strong_results_taxa <- merge(strong_results, ASVkey, by.x = "Var1", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, ASVkey, by.x = "Var2", by.y = "ASVnos")
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.x", by.y = 0)
strong_results_taxa <- merge(strong_results_taxa, tax.clean, by.x = "ASVstring.y", by.y = 0)

write.csv(strong_results_taxa, "output/Projects-ASVs-strong-results-taxa.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="Vaginal"), "output/Projects-ASVs-strong-results-taxa-vaginal.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="Pre-labor CS"), "output/Projects-ASVs-strong-results-taxa-PrelaborCS.csv", row.names = F)
write.csv(subset(strong_results_taxa, trt=="Post-labor CS"), "output/Projects-ASVs-strong-results-taxa-PostlaborCS.csv", row.names = F)


############################################################
# Now you can import these into cytoscape to visualize the network
#Open Cytoscape
#Import Network from File system
#Anything with a .x indicate as a "source node attribute"
#Anything with a .y indicate as a "target node attribute"
#"Var2" indicate as "Target Node"
#"Var1" indicate as "Source Node"
#the rest leave as "edge attribute"
#Edit>Remove Duplicated Edges, click "Ignore edge direction"
#To manually assign aesthetics, use "Discrete Mapping"
#
#
#
#
#
######################################################################

















