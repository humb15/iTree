# iTree Data Analyses

library(dada2)

# Define the path

path<- ("/Users/humbe/Desktop/iTree/16S_Fasta/RawFasta")
list.files(path)

# Get matched lists of the forward and reverse fastq files.

# Filenames have format: SAMPLENAME_R1_001_trimmed.fq and SAMPLENAME_R2_001_trimmed.fq

fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Inspect the read numbers 

library(ShortRead)
nseq <- function(f) { length(readFastq(f)) }
nseqF <- sapply(fnFs, nseq)
nseqR <- sapply(fnRs, nseq)
head(cbind(nseqF, nseqR))
table(nseqF == nseqR)

# Files do have same number of sequences

# Filter and trim
# Assign the filenames for the filtered fastq.gz files and place them in filtered/ subdirectory

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# For truncLen and maxEE parameters I used figaro in the terminal
# figaro -i ~/Desktop/iTree/16S_Fasta/RawFasta -o ~/Desktop/iTree/16S_Fasta/figaro -a 261 -f 20 --r 20
# Reads cannot be longer than the sequenced amplicon, not even 1 nt.

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(288,33),
                     maxN=0, maxEE=c(3,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # Use matchIDs = TRUE when files have different number of sequences
head(out)

# Learn error rate

errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 15)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 15)

# Visualize the estimated error rates
# Check estimated error rates (black line) are a good fit to the observed rates 
# (points), and that error rates drop with increased quality.

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample Inference

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool = "pseudo")

# Inspect the returned dada-class object:

dadaFs[[1]]
dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

# Construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

# Remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline

# Number of reads that made it through each step in the pipeline:
  
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Save final seqtab object

save(seqtab.nochim, file="/Users/humbe/Desktop/iTree/16S_Fasta/SeqTabFinal.RData")

# load("/Users/humbe/Desktop/iTree/16S_Fasta/SeqTabFinal.RData") 

# Assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "/Users/humbe/Desktop/iTree/16S_Fasta/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread=TRUE, verbose = TRUE)

save(taxa, file="/Users/humbe/Desktop/iTree/TaxaThruGenus.RData")

#taxa <- addSpecies(taxa, "/Users/humbe/Desktop/iTree/16S_Fasta/silva_species_assignment_v138.1.fa.gz",
#                   verbose = TRUE)

#taxa.print <- taxa # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print)

# Handoff to phyloseq

library(phyloseq)
library(Biostrings)
library(ggplot2)

# Extract sample names from seqtab.nochim
samples.out <- rownames(seqtab.nochim)

sample_data <- read.csv("~/Desktop/iTree/itree_metadata2.csv")

# Extract the numeric part from the sample names in seqtab.nochim
numeric_samples_otu <- as.numeric(sub("i", "", rownames(seqtab.nochim)))

# Extract the numeric part from the sample names in sample_data
numeric_samples_sample_data <- sample_data$sample

# Sort the sample names in seqtab.nochim and sample_data
seqtab.nochim <- seqtab.nochim[order(numeric_samples_otu), ]
sample_data <- sample_data[order(numeric_samples_sample_data), ]

rownames(seqtab.nochim) <- sample_data$sample

identical(rownames(seqtab.nochim), sample_data$sample)

# Check the data and structure to ensure they match
str(taxa)
str(sample_data)
str(seqtab.nochim)

DataSet <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                    sample_data(sample_data), 
                    tax_table(taxa))

DataSet

dna_Set <- Biostrings::DNAStringSet(taxa_names(DataSet))
names(dna_Set) <- taxa_names(DataSet)
DataSet <- merge_phyloseq(DataSet, dna_Set)
taxa_names(DataSet) <- paste0("ASV", seq(ntaxa(DataSet)))
DataSet

#If there are ASVs with no counts, we should remove them 

DataSet <- prune_samples(sample_sums(DataSet) > 0, DataSet)
DataSet <- subset_samples(DataSet, sample_sums(DataSet) > 1e-5)

# Filter ASVs for which the variance across all samples is very low

varianceThreshold <- 1e-5
x <- taxa_sums(DataSet)
keepTaxa <- (x / sum(x)) > varianceThreshold
DataSet <- prune_taxa(keepTaxa, DataSet)

# Check if there are ASVs with no counts and how many there are.
# If there are ASVs with no counts, we should remove them.

any(taxa_sums(DataSet) == 0)
sum(taxa_sums(DataSet) == 0)

# Remove ASVs with less than a total of 50 counts.

DataSet <- prune_taxa(taxa_sums(DataSet) > 50, DataSet)

# Define prevalence threshold as 0.1% of total samples

prevalenceThreshold <- 0.001 * nsamples(DataSet)
prevalenceThreshold

# Compute prevalence of each feature, store as data.frame
# Counts the number of elements that are greater than zero

prevdf <- apply(X = otu_table(DataSet),
                MARGIN = ifelse(taxa_are_rows(DataSet), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame

prevdf <- data.frame(Prevalence = prevdf,
                     Abundance = taxa_sums(DataSet),
                     tax_table(DataSet))

# Execute prevalence filter, using `prune_taxa()` function

TaxatoKeep <- rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
DataSet<- prune_taxa(TaxatoKeep, DataSet)

# Number of taxa that are kept

ntaxa(DataSet)

# Number of samples that are kept

nsamples(DataSet)

head(sample_sums(DataSet))

# Check distribution

SeqDepth <- rowSums(otu_table(DataSet))

sort(SeqDepth)
min(SeqDepth)
max(SeqDepth)

# Number of unique ASVs per sample

rowSums(otu_table(DataSet) != 0)

# Total number of reads

sum(otu_table(DataSet))

# Extract abundance matrix from the phyloseq object
ASV_DataSet <- as(otu_table(DataSet), "matrix")
# transpose if necessary
if(taxa_are_rows(DataSet)){ASV_DataSet <- t(ASV_DataSet)}
# Coerce to data.frame
ASVdf_DataSet <- as.data.frame(ASV_DataSet)

# Extract taxa matrix from the phyloseq object
Taxa_DataSet <- as(tax_table(DataSet), "matrix")
# transpose if necessary
if(taxa_are_rows(DataSet)){Taxa_DataSet <- t(Taxa_DataSet)}
# Coerce to data.frame
TaxaDF_DataSet <- as.data.frame(Taxa_DataSet)

# Save to export to Excel

library("openxlsx")
write.xlsx(ASVdf_DataSet, "/Users/humbe/Desktop/iTree/ASV_DataSet.xlsx", rowNames = TRUE)
write.xlsx(TaxaDF_DataSet, "/Users/humbe/Desktop/iTree/TaxaDF_DataSet.xlsx", rowNames = TRUE)

# Save pruned DataSet

save(DataSet, file="~/Desktop/iTree/DataSet.RData")

# Summarize the contents

library(microbiome)
summarize_phyloseq(DataSet)

# Check the structure of your sample data to ensure 'quality' is a category
sample_data(DataSet)

# Calculate the average read number per quality type
avg_read_per_quality <- aggregate(sample_sums(DataSet), by = list(quality = sample_data(DataSet)$quality), FUN = mean)

# Print the result
print(avg_read_per_quality)

#Abundances for taxonomic groups ('OTU table') as a Taxa x Samples matrix:

# Absolute abundances
otu.absolute <- abundances(DataSet)

# Relative abundances
otu.relative <- abundances(DataSet, "compositional")

# List unique order-level groups:

sort(get_taxa_unique(DataSet, "Order"))

# Alpha Diversity

# Plot 

# Plot richness measures with boxplots

# Define colors for each quality level
my_colors <- c('#FEDC56', '#A1DAB9', '#3F51B5')

# Plot richness measures with boxplots
p <- plot_richness(DataSet, x = "quality", color = "quality", measures = c("Chao1", "Shannon", "Simpson")) + 
  geom_boxplot(outlier.shape = NA, outlier.colour = NA) +
  
  # Remove legend
  guides(color = FALSE) +
  
  # Assign custom colors
  scale_color_manual(values = my_colors) +
  
  # Adjust facet labels
  labs(x = "Quality", y = "Richness") +
  
  # Override facet labels with desired names
  scale_x_discrete(labels = c("1high" = "High", "2mid" = "Mid", "3low" = "Low")) +
  
  # Further customization as needed
  theme_minimal()

print(p)

rich_T <- estimate_richness(DataSet, measures = c("Shannon", "Simpson", "Chao1"))
rownames(rich_T) <- DataSet@sam_data[["sample"]]
rich_T

# Investigate if data is normal: plot each metric.

par(mfrow=c(2,2))

hist(rich_T$Hill, main="Shannon diversity", xlab="", breaks=10)
hist(rich_T$Simpson, main="Simpson diversity", xlab="", breaks=10)
hist(rich_T$InvSimpson, main="Inv Simpson diversity", xlab="", breaks=10)
hist(rich_T$Chao1, main="Chao richness", xlab="", breaks=15)

par(mfrow=c(1,1))

# Density plot

# library("ggpubr")

# ggdensity(RD_L$Shannon, 
#main = "Shannon diversity",
#xlab = "")

# Q-Q Plot

qq1<- ggqqplot(rich_T$Shannon, title = "Shannon")
qq2<- ggqqplot(rich_T$Hill, title = "Hill")
qq3<- ggqqplot(rich_T$InvSimpson, title = "Inv Simpson")
qq4<- ggqqplot(rich_T$Chao1, title = "Chao")

plot_grid(qq1, qq2, qq3, qq4,
          ncol=2, nrow=2)

# If all the points fall approximately along this reference line, we can assume normality.

# You want the data to be roughly normal so that you can run ANOVA or t-tests

# To test for normality statistically, we can run the Shapiro-Wilk test of normality.

shapiro.test(rich_T$Shannon)
shapiro.test(rich_T$Simpson) # not normal
shapiro.test(rich_T$Chao1)

rich_T <- merge(rich_T, DataSet@sam_data, by=0, all=TRUE) 

anova_shannon <- aov(Shannon ~ quality, data = rich_T)
anova_chao1 <- aov(Chao1 ~ quality, data = rich_T)
kw_test_simpson <- kruskal.test(Simpson ~ quality, data = rich_T)

# Print the results
print(summary(anova_shannon))
print(summary(anova_chao1))
print(kw_test_simpson)

#QSeq

library(devtools)
library(QSeq)

Q.ps<-QSeq(DataSet, "x16S")

# Analyses of data

library(vegan)

# Multivariate statistical analysis (PERMANOVA)

# Beta Diversity

# NMDS

# Com_mat <- as.matrix(otu_table(DataSet))

set.seed(1992)
nmdsB <- metaMDS(otu_table(DataSet), distance = "bray", trymax = 100)

# Extracting scores

data.scores <- scores(nmdsB, display="sites", tidy = TRUE)

# Add columns to data frame 

data.scores <- merge(data.scores, DataSet@sam_data, by=0, all=TRUE)  # by=0 equals by="row.names"

# Plotting with ggplot

gg <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = quality), size = 3, alpha = 0.5) + 
  #scale_color_gradient(low = "lightgreen", high = "#FAAED2") +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Quality")

# Change labels manually

gg <- ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = factor(quality, levels = c("1high", "2mid", "3low"))), size = 3, alpha = 0.5) + 
  scale_color_manual(values = c("1high" = "salmon", "2mid" = "lightgreen", "3low" = "steelblue"), 
                     labels = c("1high" = "High", "2mid" = "Mid", "3low" = "Low")) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Quality")

gg


#Biplot arrows

# Select the desired metrics

env <- data.scores[, c("quality_score", "pH", "CEC", "Mg", "Ca")]
colnames(env)[colnames(env) == "quality_score"] <- "Quality score"

en <- envfit(nmdsB, env, permutations = 999, na.rm = TRUE)

# Access the arrow coordinates
arrow_data <- as.data.frame(en[["vectors"]][["arrows"]])
arrow_labels <- rownames(arrow_data)

# Make the CEC arrow shorter (e.g., multiply by 0.4)
arrow_data["CEC", ] <- arrow_data["CEC", ] * 0.9

# Add the biplot arrows to the existing plot 'gg'

gg +
  geom_segment(data = arrow_data, 
               aes(x = 0, y = 0, xend = arrow_data[, 1], yend = arrow_data[, 2]), 
               arrow = arrow(type = "closed", length = unit(0.1, "inches")),
               color = c("red", "green", "blue", "purple", "orange"), size = 0.5) +
  geom_text(data = arrow_data, aes(x = arrow_data[, 1], y = arrow_data[, 2], label = arrow_labels), 
            size = 3, nudge_x = c(0, -0.05, 0.05, 0, 0)) +
  theme_minimal()

# Statistical analyses

# Transform data to proportions as appropriate for Bray-Curtis distances 

Q.ps_prop <- transform_sample_counts(Q.ps, function(otu) otu / sum(otu))

BC.dist <- phyloseq::distance(Q.ps_prop, method = "bray")

# Calculate dispersion (variances) within each group.

var.metr <- betadisper(BC.dist, as.factor(Q.ps@sam_data[["quality"]]), type=c("median"))

# Perform an ANOVA-like test to determine if the variances differ by groups.
# H0: no difference between the groups being compared.

permutest(var.metr, permutations=1000)

# within-group variances are equal or homogenous so I can follow up with a permanova
# partially true, permanova doesn't make the asumption but unequal variances can affect

plot(var.metr, hull = FALSE, ellipse = TRUE, main = NULL)

var.metr

# Run PERMANOVA on distances.
# H0: no difference in multivariate dispersion between groups

adonis2(BC.dist ~ DataSet@sam_data$quality, permutations = 1000, method = "bray")

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

pairwise_adonis <- pairwise.adonis2(BC.dist ~ quality, data = as(sample_data(DataSet), "data.frame"), nperm = 999)

# Print the results
print(pairwise_adonis)

# Prepare and subset data

library(plyr)

Ab_Fam <-tax_glom(DataSet,taxrank = "Family", NArm=FALSE)

AbFam_melt <- psmelt(Ab_Fam)  

sub_AbFam <- subset(AbFam_melt, select = c("Sample", "Abundance", "quality", "Family"))

# Merge family sums and the total numbers

Totals <- ddply(sub_AbFam, c("Sample"), summarise, total = sum(Abundance))  

fam_totals <- merge(sub_AbFam, Totals, by = "Sample")  

# Calculate the relative abundance (how is the total abundance of a specific 
# family distributed among samples)

fam_totals$RelAbundance <- fam_totals$Abundance/fam_totals$total  

#  Calculate the Percent Abundance

fam_totals$PercentAbund <- fam_totals$RelAbundance * 100  

#  Calculate the family abundance based on quality.

fam_stats <- ddply(fam_totals, c("quality","Family"), summarise, 
                   N = length(PercentAbund),
                   mean_abundance = mean(PercentAbund),
                   sd   = sd(PercentAbund),
                   se   = sd / sqrt(N))

abund_by_fam <- ddply(fam_stats, c("Family"), summarise, 
                      N = length(mean_abundance),
                      family_mean = mean(mean_abundance))

abund_by_fam <-arrange(abund_by_fam, desc(family_mean))

# Create a vector of the families 

fam_order <- as.character(abund_by_fam$Family) 

# Top 12 most abundant families (>1.5%)

fam12_order <- as.character(abund_by_fam$Family)[1:12]

top12fams <- fam_stats[fam_stats$Family %in% fam_order[1:12], ] 

top12fams$Family <- factor(top12fams$Family,
                           levels = c(fam12_order))

# Plot

# Define colors for bars
my_colors <- c('#FFFCC6', '#A1DAB4', '#3F51B5')  

# Modify facet labels
top12fams$quality <- factor(top12fams$quality, levels = c("1high", "2mid", "3low"), 
                            labels = c("High", "Mid", "Low"))

# Plotting code
ggplot(data = top12fams, aes(x = mean_abundance, y = reorder(Family, desc(Family)), fill = quality)) +
  geom_bar(stat = 'identity', color = "black") +
  geom_errorbar(aes(xmin = mean_abundance - se, xmax = mean_abundance + se), width = 0.25, color = "black") +
  scale_fill_manual(values = my_colors) +  # Apply custom colors
  xlab("Mean Percent Relative Abundance (%)") +
  facet_wrap(~ quality) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1), legend.position="none") +
  ylab("Family")

# ITS

# Modify facet labels
top10fams$quality <- factor(top10fams$quality, levels = c("1high", "2mid", "3low"), 
                            labels = c("High", "Mid", "Low"))

# Define colors for bars
my_colors <- c('#FFFCC6', '#A1DAB4', '#3F51B5')
#my_colors <- c("#2ca25f", "#99d8c9", "#e5f5f9")

# Plotting code with modified Family names
ggplot(data = top10fams, aes(x = mean_abundance, y = reorder(gsub("^f__", "", Family), desc(Family)), fill = quality)) +
  geom_bar(stat = 'identity', color = "black") +
  geom_errorbar(aes(xmin = mean_abundance - se, xmax = mean_abundance + se), width = 0.25, color = "black") +
  scale_fill_manual(values = my_colors) +  # Apply custom colors
  xlab("Mean Percent Relative Abundance (%)") +
  facet_wrap(~ quality) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1), legend.position = "none") +
  ylab("Family")

library(microeco)

# load 16S data
# load("/Users/humbe/Desktop/iTree/16S_Fasta/DataSet.RData")

t_table <- as.data.frame(DataSet@tax_table)
colnames(t_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus")

# create microtable object
meco_bact <- microtable$new(sample_table = DataSet@sam_data, otu_table = as.data.frame(t(DataSet@otu_table)), t_table)

# use tidy_dataset() to make OTUs and samples information consistent across files
meco_bact$tidy_dataset()

t1 <- trans_func$new(meco_bact)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)

# use list to prepare data
tmp <- list()

# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp$func <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)

# assign the list as taxa_abund in your microtable object
meco_bact$taxa_abund <- tmp

# use trans_diff class to perform differential test
t2 <- trans_diff$new(dataset = meco_bact, method = "anova", group = "quality", taxa_level = "all")
t2$plot_diff_abund(add_sig = T) + ggplot2::ylab("Relative abundance (%)")
