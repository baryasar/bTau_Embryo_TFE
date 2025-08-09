library(tidyr)
library(tibble)
library(ggplot2)
library(viridis)

# NOTE: Some of the objects were pre-loaded in the environment. These objects were loaded or created in other scripts within this project.

# Extract scaled data for RNA assay
scaled_data <- as.data.frame(GetAssayData(strt.seurat.obj, assay = "RNA", layer = "scale.data"))

set.seed(20)
kClust <- kmeans(scaled_data, centers=6, nstart = 1000, iter.max = 20)
#warnings() 
kClusters <- kClust$cluster

# Check column names to ensure they correspond to the sample identifiers
sample_names <- colnames(scaled_data)

# Get developmental stage information from the Seurat object
stage_info <- strt.seurat.obj$orig.ident[sample_names]

# Convert stage info to a vector
stage_info_vector <- as.vector(stage_info)

# Rename columns with the stage information
colnames(scaled_data) <- stage_info_vector
# Create unique column names (if necessary) to avoid duplicate column names
colnames(scaled_data) <- make.unique(colnames(scaled_data), sep = "_")

# Convert the data to long format
scaled_data_long_df <- scaled_data %>%
  rownames_to_column(var = "gene") %>%  # Keep gene names as a column
  pivot_longer(-gene, names_to = "stage", values_to = "expression")

# Add cluster information (assuming you already have it)
scaled_data_long_df$cluster <- rep(kClusters, each = nrow(scaled_data_long_df) / length(kClusters))

# Remove suffixes from stage names
scaled_data_long_df <- scaled_data_long_df %>%
  mutate(stage = sub("_\\d+$", "", stage))  # Remove the suffix

# Calculate average expression per gene and stage
average_scaled_data_long_df <- scaled_data_long_df %>%
  group_by(gene, stage,cluster) %>%
  summarize(average_expression = mean(expression, na.rm = TRUE)) %>%
  ungroup()

average_scaled_data_long_df$stage <- recode(
  average_scaled_data_long_df$stage,
  "16c" = "16-cell",
  "2c" = "2-cell",
  "4c" = "4-cell",
  "8c" = "8-cell",
  "GV" = "GV oocyte",
  "MII" = "MII oocyte"
)

average_scaled_data_long_df$stage <- factor(
  average_scaled_data_long_df$stage,
  levels = c("GV oocyte", "MII oocyte", "2-cell", "4-cell", "8-cell", "16-cell", "Blastocyst")
)

# Calculate the average expression for each cluster and stage
cluster_average <- average_scaled_data_long_df %>%
  group_by(stage, cluster) %>%
  summarize(cluster_avg_expression = mean(average_expression, na.rm = TRUE)) %>%
  ungroup()


tiff(filename="Cluster_line6_plot.tiff", width=2000, height=2250, res=300)

# Custom labeller function to add 'Cluster' before the cluster number
cluster_labeller <- function(variable, value) {
  paste("Cluster", value)
}

# Plot average expression per gene and a thick line for the cluster average
ggplot(average_scaled_data_long_df, aes(x = stage, y = average_expression, color = factor(cluster), group = gene)) +
  geom_line(alpha = 0.3) +  # Line plot for each gene (lighter lines)
  geom_line(data = cluster_average, aes(x = stage, y = cluster_avg_expression, group = cluster), 
            linewidth = 2, color = "black") +  # Thick black line for cluster average
  facet_wrap(~ cluster, scales = "free_y", labeller = cluster_labeller, , ncol = 2) +  # Custom labels for facets
  scale_color_viridis_d(option = "viridis") +  # Use turbo color palette
  theme_minimal() +
  labs(x = NULL, y = "Average Expression", title = NULL) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Close TIFF device
dev.off()
