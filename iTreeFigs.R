library(vegan)
library(ggplot2)

load("/Users/humbe/Desktop/iTree/16S_Fasta/DataSet.RData")
load("/Users/humbe/Desktop/iTree/ITS_Fasta/DataSet.RData")

# Your data
data <- data.frame(
  sample = 1:26,
  ITS = c(1.21E+08, 5.72E+07, 2.38E+08, 1.96E+08, 2.92E+08, 2.47E+08, 4.65E+07, 8.16E+07, 1.17E+08, 2.46E+08, 2.39E+07, 1.54E+08, 3.35E+07, 2.07E+08, 3.70E+07, 7.08E+06, 6.95E+06, 9.17E+07, 1.66E+08, 1.19E+07, 1.54E+08, 1.21E+07, 2.73E+07, 2.70E+08, 5.62E+07, 3.63E+07),
  MBC = c(33, 54, 76, 67, 118, 64, 209, 188, 62, -11, 74, 60, 75, 22, 20, 155, 91, 64, 156, 92, 69, 20, 120, 21, 23, -164)
)

# Define the custom order
custom_order <- c(1, 10:19, 2, 20:26, 3:9)

# Reorder the rows
reordered_data <- data[match(custom_order, data$sample), ]

DataSet@sam_data$ITS <- reordered_data$ITS
DataSet@sam_data$MBC <- reordered_data$MBC

save(DataSet, file="~/Desktop/iTree/16S_Fasta/DataSet.RData")
save(DataSet, file="~/Desktop/iTree/ITS_Fasta/DataSet.RData")

set.seed(1992)
nmdsB <- metaMDS(otu_table(DataSet), distance = "bray", trymax = 100)

# Extracting scores

data.scores <- scores(nmdsB, display="sites", tidy = TRUE)
#data.scores <- scores(nmdsB, display="species", tidy = TRUE)

# Add columns to data frame 

data.scores <- merge(data.scores, DataSet@sam_data, by=0, all=TRUE)  # by=0 equals by="row.names"

ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, colour = quality_score)) + 
  geom_point(size = 3, alpha = 0.5) +
  scale_color_viridis_c(breaks = c(min(data.scores$quality_score), median(data.scores$quality_score), max(data.scores$quality_score)),
                        labels = c("Low", "Mid", "High")) +
  theme(
    axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
    panel.background = element_blank(), 
    panel.border = element_rect(fill = NA, colour = "grey30"), 
    axis.ticks = element_blank(), 
    axis.text = element_blank(), 
    legend.key = element_blank(), 
    legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  labs(colour = "Quality")

env <- data.scores[, c("quality_score", "x16S", "OM","pH","Ca", "MBC")]
env <- data.scores[, c("quality_score", "MBC", "ITS", "OM","pH","Ca")]

colnames(env)[colnames(env) == "quality_score"] <- "Quality score"
colnames(env)[colnames(env) == "x16S"] <- "Gene copy"
colnames(env)[colnames(env) == "ITS"] <- "Gene copy"

# Extract site coordinates from nmdsB for 16S
#site_coordinates <- nmdsB$species

# Ensure the site coordinates are in the correct format (data frame with two columns)
#site_coordinates <- as.data.frame(site_coordinates)

# Ensure row names match those in env
#rownames(site_coordinates) <- rownames(env)

# Now, you can use envfit with the site coordinates
#en <- envfit(site_coordinates, env, permutations = 999, na.rm = TRUE)

en <- envfit(nmdsB, env, permutations = 999, na.rm = TRUE)

en

# Access the arrow coordinates, they're already in the 'en' object
#arrow_data <- as.data.frame(en[["vectors"]][["arrows"]])
#arrow_data

# Create a list of factor names, their corresponding NMDS coordinates, and p-values
factors <- list(
  "Quality score" = list(coords = c(-0.9988521, -0.04790144), p_value = "0.005"),
  "MBC" = list(coords = c(-0.1821004, 0.98327994), p_value = "0.196"),
  "Gene copy" = list(coords = c(-0.1768432, -0.98423904), p_value = "0.844"),
  "OM" = list(coords = c(-0.8015206, 0.59796718), p_value = "0.036"),
  "pH" = list(coords = c(0.9970344, 0.07695780), p_value = "0.001"),
  "Ca" = list(coords = c(0.9634637, 0.26783881), p_value = "0.001")
)

#FOR ITS:

factors <- list(
  "Quality score" = list(coords = c(-0.9242046, -0.3818978), p_value = "0.001"),
  "MBC" = list(coords = c(-0.3374578, 0.9413407), p_value = "0.622"),
  "Gene copy" = list(coords = c(-0.6344845, 0.7729356), p_value = "0.004"),
  "OM" = list(coords = c(-0.9876030, -0.1569724), p_value = "0.247"),
  "pH" = list(coords = c(0.9949107, 0.1007602), p_value = "0.001"),
  "Ca" = list(coords = c(0.9836450, 0.1801182), p_value = "0.001")
)

# Set up the layout for multiple plots
par(mfrow = c(2, 3))  # 2 rows, 3 columns

# Loop through each factor and create a plot
for (factor_name in names(factors)) {
  plot(x = c(-1, 1), y = c(-1, 1), type = "n", xlab = "", ylab = "", main = factor_name)
  lines(c(-1, 1), c(0, 0), col = "gray")  # X-axis line
  lines(c(0, 0), c(-1, 1), col = "gray")  # Y-axis line
  arrows(0, 0, factors[[factor_name]]$coords[1], factors[[factor_name]]$coords[2], col = "red", length = 0.1)
  mtext(side = 1, text = paste("p-value:", factors[[factor_name]]$p_value), adj = 0.5, line = 2, cex = 0.8)  # Add p-value as subtitle centered
}

library(viridis)

# Extract three colors from the Viridis color palette
low_color <- viridis(3)[1]  # Low value color
mid_color <- viridis(3)[2]  # Mid value color
high_color <- viridis(3)[3] # High value color

"#FDE725FF"
"#21908CFF"
"#440154FF"


ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = factor(quality, levels = c("1high", "2mid", "3low"))), size = 3, alpha = 0.5) + 
  scale_color_manual(values = c("1high" = "#FDE725FF", "2mid" = "#43A047", "3low" = "#1A237E"), 
                     labels = c("1high" = "High", "2mid" = "Mid", "3low" = "Low")) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Quality")

#------------

# Extract sample data from phyloseq object
sam_data <- as.data.frame(sample_data(DataSet))

# Extract relevant variables
variables <- c("quality", "quality_score", "MBC", "ITS", "x16S", "OM", "pH", "M3_P", "K", "Mg", "Ca")
relevant_data <- sam_data[, variables]

# Rename 'quality' levels
relevant_data$quality <- factor(relevant_data$quality, levels = c("1high", "2mid", "3low"), labels = c("High", "Mid", "Low"))

library(dplyr)

# Initialize a list to store results
test_results <- list()

# Loop through each variable and perform Shapiro-Wilk test before Kruskal-Wallis test
for (variable in variables) {
  var_data <- relevant_data[[variable]]
  if (is.numeric(var_data)) {
    shapiro_result <- shapiro.test(var_data)
    if (shapiro_result$p.value > 0.05) {
      # If data is normally distributed, perform ANOVA
      test_result <- aov(var_data ~ relevant_data$quality)
      test_results[[variable]] <- summary(test_result)
    } else {
      # If data is not normally distributed, perform Kruskal-Wallis test
      test_result <- kruskal.test(var_data ~ relevant_data$quality)
      test_results[[variable]] <- test_result
    }
  }
}

# Print out the results
for (variable in names(test_results)) {
  cat("Variable:", variable, "\n")
  print(test_results[[variable]])
  cat("\n")
}

# Group by 'quality' and calculate mean and SD for each variable
stats <- relevant_data %>%
  group_by(quality) %>%
  summarise(
    quality_score_mean = mean(quality_score, na.rm = TRUE),
    quality_score_sd = sd(quality_score, na.rm = TRUE),
    MBC_mean = mean(MBC, na.rm = TRUE),
    MBC_sd = sd(MBC, na.rm = TRUE),
    x16S_mean = mean(x16S, na.rm = TRUE),
    x16S_sd = sd(x16S, na.rm = TRUE),
    ITS_mean = mean(ITS, na.rm = TRUE),
    ITS_sd = sd(ITS, na.rm = TRUE),
    OM_mean = mean(OM, na.rm = TRUE),
    OM_sd = sd(OM, na.rm = TRUE),
    pH_mean = mean(pH, na.rm = TRUE),
    pH_sd = sd(pH, na.rm = TRUE),
    M3_P_mean = mean(M3_P, na.rm = TRUE),
    M3_P_sd = sd(M3_P, na.rm = TRUE),
    K_mean = mean(K, na.rm = TRUE),
    K_sd = sd(K, na.rm = TRUE),
    Mg_mean = mean(Mg, na.rm = TRUE),
    Mg_sd = sd(Mg, na.rm = TRUE),
    Ca_mean = mean(Ca, na.rm = TRUE),
    Ca_sd = sd(Ca, na.rm = TRUE)
  )

# View the summary statistics
print(stats)

# Group by 'quality' and calculate mean and SE for each variable
stats_se <- relevant_data %>%
  group_by(quality) %>%
  summarise(
    quality_score_mean = mean(quality_score, na.rm = TRUE),
    quality_score_se = sd(quality_score, na.rm = TRUE) / sqrt(n()),
    MBC_mean = mean(MBC, na.rm = TRUE),
    MBC_se = sd(MBC, na.rm = TRUE) / sqrt(n()),
    x16S_mean = mean(x16S, na.rm = TRUE),
    x16S_se = sd(x16S, na.rm = TRUE) / sqrt(n()),
    ITS_mean = mean(ITS, na.rm = TRUE),
    ITS_se = sd(ITS, na.rm = TRUE) / sqrt(n()),
    OM_mean = mean(OM, na.rm = TRUE),
    OM_se = sd(OM, na.rm = TRUE) / sqrt(n()),
    pH_mean = mean(pH, na.rm = TRUE),
    pH_se = sd(pH, na.rm = TRUE) / sqrt(n()),
    M3_P_mean = mean(M3_P, na.rm = TRUE),
    M3_P_se = sd(M3_P, na.rm = TRUE) / sqrt(n()),
    K_mean = mean(K, na.rm = TRUE),
    K_se = sd(K, na.rm = TRUE) / sqrt(n()),
    Mg_mean = mean(Mg, na.rm = TRUE),
    Mg_se = sd(Mg, na.rm = TRUE) / sqrt(n()),
    Ca_mean = mean(Ca, na.rm = TRUE),
    Ca_se = sd(Ca, na.rm = TRUE) / sqrt(n())
  )

# Define a function to perform pairwise comparisons using Wilcoxon rank sum test
perform_pairwise_comparison <- function(data, variable) {
  # Extract data for the variable
  var_data <- data[[variable]]
  
  # Perform pairwise comparison using Wilcoxon rank sum test
  pairwise_test <- pairwise.wilcox.test(var_data, data$quality, p.adjust.method = "bonferroni")
  
  return(pairwise_test)
}

# List of variables to analyze
variables <- names(relevant_data)[-1]

# Initialize a list to store results
comparison_results <- list()

# Loop through each variable and perform pairwise comparisons
for (variable in variables) {
  comparison_results[[variable]] <- perform_pairwise_comparison(relevant_data, variable)
}

# View the results
comparison_results

library(dplyr)
library(tidyr)

# Create the plot data
plot_data <- reshape2::melt(cbind.data.frame(sampname = rownames(t2$res_spe_func_perc), 
                                             t2$res_spe_func_perc), 
                            id.vars = "sampname")

# Extract big categories and subcategories
plot_data <- plot_data %>%
  separate(variable, into = c("category", "subcategory"), sep = "\\|")

# Capitalize the first letter of the category column
plot_data$category <- tools::toTitleCase(plot_data$category)

# Replace underscores with spaces in the subcategory column
plot_data$subcategory <- gsub("_", " ", plot_data$subcategory)

# Define order and quality groups
order_x <- paste0("i", 1:26)
quality_groups <- cut(as.numeric(gsub("i", "", order_x)), 
                      breaks = c(0, 10, 18, 26), 
                      labels = c("High", "Mid", "Low"))

# Update quality groups for each sample
plot_data$quality <- ifelse(as.numeric(gsub("i", "", plot_data$sampname)) %in% 1:10, "High",
                            ifelse(as.numeric(gsub("i", "", plot_data$sampname)) %in% 11:18, "Mid", "Low"))

# Reorder levels of the quality variable
plot_data$quality <- factor(plot_data$quality, levels = c("High", "Mid", "Low"))

# Reorder levels of the 'category' variable
plot_data$category <- factor(plot_data$category, 
                             levels = c("Primary Lifestyle", 
                                        "Endophytic Interaction Capability", 
                                        "Plant Pathogenic Capacity", 
                                        "Ectomycorrhiza Exploration Type"))

# Create the plot
heatmap_plot <- ggplot(aes(x = factor(sampname, levels = order_x), 
                           y = factor(subcategory, levels = rev(unique(subcategory))), 
                           fill = value), 
                       data = plot_data) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#00695C") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(color = "black", fill = "white"),  # Add border around facets
        strip.placement = "outside",  # Move facet labels outside the plot area
        axis.text.y = element_text(size = 10), 
        axis.title.y = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),  # Adjust orientation of big categories
        strip.text.y = element_text(size = 10, angle = 0)) +  # Adjust orientation of big categories
  scale_y_discrete(position = "right") +  # Move y-axis ticks to the right side of the plot
  facet_grid(category ~ quality, scales = "free", space = "free", switch = "y") +  # Facet big 4 groups on the left
  theme(strip.text.y.left = element_text(angle = 0))  # Adjust orientation of facet labels

# Print the plot
print(heatmap_plot)

# Bacterial dataset
library(reshape2)
library(microeco)

t_table <- as.data.frame(DataSet@tax_table)
colnames(t_table) <- c("Kingdom","Phylum","Class","Order","Family","Genus")

# create microtable object
meco_bact <- microtable$new(sample_table = DataSet@sam_data, otu_table = as.data.frame(t(DataSet@otu_table)), t_table)

# use tidy_dataset() to make OTUs and samples information consistent across files
meco_bact$tidy_dataset()

t1 <- trans_func$new(meco_bact)
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = FALSE)
t1$plot_spe_func_perc(order_x = paste0(1:26))

# If you want to change the group list, reset the list t2$func_group_list
t1$func_group_list

# use show_prok_func to see the detailed information of prokaryotic traits
t1$show_prok_func("methanotrophy")


#---------

# Extract abundance data
abundance_data <- t1$res_spe_func_perc

# Create a column named "Sample" and copy row names
abundance_data$Sample <- rownames(abundance_data)

# Add 'i' before each sample ID
abundance_data$Sample <- paste0("i", abundance_data$Sample)

# Classify samples based on quality
abundance_data$Quality <- ifelse(as.numeric(gsub("[^0-9]", "", abundance_data$Sample)) <= 10, "High",
                                 ifelse(as.numeric(gsub("[^0-9]", "", abundance_data$Sample)) <= 18, "Mid", "Low"))

# Melt the data for plotting
plot_data <- melt(abundance_data, id.vars = c("Sample", "Quality"), variable.name = "Category", value.name = "Relative_Abundance")

# Calculate total relative abundance for each category
category_totals <- aggregate(Relative_Abundance ~ Category, data = plot_data, FUN = sum)

# Sort categories based on total relative abundance and select the top 13 (>10)
top_categories <- head(category_totals[order(category_totals$Relative_Abundance, decreasing = TRUE), ], 13)

# Filter plot_data to include only top categories
plot_data_filtered <- plot_data[plot_data$Category %in% top_categories$Category, ]

library(stringr)

# Capitalize first letter and replace underscores with spaces in Category column
plot_data_filtered$Category <- str_replace_all(str_to_title(plot_data_filtered$Category), "_", " ")

# Reorder the levels of the Quality column
plot_data_filtered$Quality <- factor(plot_data_filtered$Quality, levels = c("High", "Mid", "Low"))

# Convert Sample column to factor with correct order
plot_data_filtered$Sample <- factor(plot_data_filtered$Sample, levels = paste0("i", 1:26))

# Plot with the reordered Sample column
heatmap <- ggplot(plot_data_filtered, aes(x = Sample, y = Category, fill = Relative_Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#558B2F") +
  facet_grid(. ~ Quality, scales = "free_x", space = "free_x") +
  theme_minimal() +
  labs(x = "Sample", y = "Category", fill = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

# Print the heatmap
print(heatmap)


# By quality

# Extract relative abundance information per group in each category
abundance_info <- t2$res_abund

# Make changes to the Taxa column
abundance_info$Taxa <- gsub("_", " ", abundance_info$Taxa)
abundance_info$Taxa <- tools::toTitleCase(abundance_info$Taxa)

# Make changes to the Group column
abundance_info$Group <- gsub("1high", "High", abundance_info$Group)
abundance_info$Group <- gsub("2mid", "Mid", abundance_info$Group)
abundance_info$Group <- gsub("3low", "Low", abundance_info$Group)

# Create the heatmap
heatmap <- ggplot(abundance_info, aes(x = Group, y = Taxa, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#558B2F") +
  labs(x = "Quality Group", y = "Category", fill = "Mean Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "right")

# Display the heatmap
print(heatmap)

# Calculate total abundance for each category
library(dplyr)

total_abundance <- abundance_info %>%
  group_by(Taxa) %>%
  summarise(Total_Abundance = sum(Mean)) %>%
  arrange(desc(Total_Abundance))  # Arrange in descending order of total abundance

# Select the top 15 
top_categories <- total_abundance$Taxa[1:15]  # Adjust the number as needed

# Filter abundance_info to keep only the top categories
filtered_abundance_info <- abundance_info %>%
  filter(Taxa %in% top_categories)

# Create the heatmap
heatmap <- ggplot(filtered_abundance_info, aes(x = Group, y = Taxa, fill = Mean)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#558B2F") +
  labs(x = "Quality Group", y = "Category", fill = "Mean Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "right")

# Display the heatmap
print(heatmap)


# Define a function to group categories
group_categories <- function(category) {
  if (category %in% c("Photoheterotrophy", "Phototrophy", "Chemoheterotrophy", "Aerobic chemoheterotrophy")) {
    return("Energy source")
  } else if (category %in% c("Aerobic ammonia oxidation", "Nitrification", "Nitrogen fixation", "Nitrate respiration", 
                             "Nitrate reduction", "Nitrogen respiration")) {
    return("N - cycle")
  } else if (category == "Cellulolysis") {
    return("C - cycle")
  } else {
    return("Other")
  }
}

# Apply the function to create the new grouping variable
plot_data_filtered$Category_Group <- sapply(plot_data_filtered$Category, group_categories)

# Capitalize first letter and replace underscores with spaces in Category column
plot_data_filtered$Category <- str_replace_all(str_to_title(plot_data_filtered$Category), "_", " ")

# Reorder the levels of the Quality column
plot_data_filtered$Quality <- factor(plot_data_filtered$Quality, levels = c("High", "Mid", "Low"))

# Convert Sample column to factor with correct order
plot_data_filtered$Sample <- factor(plot_data_filtered$Sample, levels = paste0("i", 1:26))

# Plot with the reordered Sample column and faceted by the new Category_Group
heatmap <- ggplot(plot_data_filtered, aes(x = Sample, y = factor(Category, levels = rev(unique(Category))), fill = Relative_Abundance)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#558B2F") +
  facet_grid(Category_Group ~ Quality, scales = "free", space = "free") +  # Facet by Category_Group
  theme_minimal() +
  labs(x = "Sample", y = "Category", fill = "Relative Abundance (%)") +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))

# Print the heatmap

print(heatmap)

