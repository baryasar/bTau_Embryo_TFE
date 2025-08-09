library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubfigs)
library(RColorBrewer)
library(svglite)
library(viridis)
library(scales)

################################################################################
### 1. Stacked Bar Plot of Annotation Percentages (Averaged per Stage)
################################################################################

# Load annotation counts and remove unwanted columns and samples
Annoations_stages_counts <- read.delim("BOV_LIB_1.3_byTFE-counts_annotation.txt")

Annoations_stages_counts <- Annoations_stages_counts[,-c(1,2,4,5,6,7,8,9)]

# Remove specific samples that are low quality or unwanted
Annoations_stages_counts <- Annoations_stages_counts %>%
  select(-c("X16c.BOV_LIB_1.3_37","X8c.BOV_LIB_1.3_26", "GV.BOV_LIB_1.3_11", "GV.BOV_LIB_1.3_13","X2c.BOV_LIB_1.3_23"))

# Sum counts across all samples for each annotation category
Annoations_stages_counts_summary <- Annoations_stages_counts %>%
  group_by(Annotation) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))

# Calculate percentage contribution of each annotation category per stage/sample
Annoations_percentage <- Annoations_stages_counts_summary %>%
  mutate(across(-Annotation, ~ .x / sum(.x) * 100))

# Create new columns for each prefix by calculating the row-wise mean
Annoations_percentage_avg <- Annoations_percentage %>%
  mutate(
    X16c_avg = rowMeans(select(., starts_with("X16c")), na.rm = TRUE),
    X2c_avg = rowMeans(select(., starts_with("X2c")), na.rm = TRUE),
    X4c_avg = rowMeans(select(., starts_with("X4c")), na.rm = TRUE),
    X8c_avg = rowMeans(select(., starts_with("X8c")), na.rm = TRUE),
    Blc_avg = rowMeans(select(., starts_with("Blc")), na.rm = TRUE),
    GV_avg = rowMeans(select(., starts_with("GV")), na.rm = TRUE),
    MII_avg = rowMeans(select(., starts_with("MII")), na.rm = TRUE)
  ) %>%
  select(Annotation, X16c_avg, X2c_avg, X4c_avg, X8c_avg, Blc_avg, GV_avg, MII_avg)

# Group similar annotations together for simplified categories for plotting
Annoations_percentage_avg_6_feature <- Annoations_percentage_avg %>%
  mutate(Annotation = case_when(
    Annotation %in% c("Coding_5UTR", "Coding_upstream") ~ "Upstream and 5′-UTR of coding",
    Annotation %in% c("Coding_3UTR", "Coding_CDS") ~ "CDS and 3′-UTR of coding",
    Annotation %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Upstream and 1st exon of noncoding",
    TRUE ~ Annotation
  )) %>%
  group_by(Annotation) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

# Rename columns for better readability in plots
colnames(Annoations_percentage_avg_6_feature) <- c("Annotation", "16-cell","2-cell","4-cell","8-cell","Blastocyst","GV oocyte","MII oocyte")

# Convert the wide-format data to long-format
Annoations_long <- Annoations_percentage_avg_6_feature %>%
  pivot_longer(cols = -Annotation, names_to = "Stage", values_to = "Percentage")

# Specify the desired order of stages
desired_order <- c("GV oocyte", "MII oocyte", "2-cell", "4-cell", "8-cell","16-cell", "Blastocyst")

# Convert the Stage column to a factor with specified levels
Annoations_long <- Annoations_long %>%
  mutate(Stage = factor(Stage, levels = desired_order))

# Set order and rename factor levels for Annotation
Annoations_long$Annotation <- factor(Annoations_long$Annotation, levels=c("Unannotated","Intron" ,"Noncoding_other-exon","Upstream and 1st exon of noncoding","CDS and 3′-UTR of coding" ,"Upstream and 5′-UTR of coding" ))
Annoations_long$Annotation <- as.factor(Annoations_long$Annotation)
Annoations_long$Annotation <- recode(Annoations_long$Annotation, 
                                      "Noncoding_other-exon" = "Noncoding other exon")

svglite(filename = "Annotation_percentages_average.svg", width = 15, height = 9)

ggplot(Annoations_long, aes(x = Stage, y = Percentage, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +
  scale_y_continuous(labels = scales::percent) +  # Format y-axis as percentage
  labs(
    x = "Stage",
    y = "Aligned reads",
    fill = "Annotation",
    title = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 22),  # Increase font size of x-axis labels
    axis.text.y = element_text(size = 18),  # Increase y-axis label size
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_text(size = 22),  # Increase y-axis title size
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend elements font size
    legend.position = "bottom"
  )

dev.off()

################################################################################
### 2. Boxplot + Dotplot Showing Per-sample Variability of Annotation Percentages
################################################################################


Annoations_percentage_6feature <- Annoations_percentage %>%
  mutate(Annotation = case_when(
    Annotation %in% c("Coding_5UTR", "Coding_upstream") ~ "Upstream and 5′-UTR of coding",
    Annotation %in% c("Coding_3UTR", "Coding_CDS") ~ "CDS and 3′-UTR of coding",
    Annotation %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Upstream and 1st exon of noncoding",
    TRUE ~ Annotation
  )) %>%
  group_by(Annotation) %>%
  summarise(across(everything(), sum, na.rm = TRUE), .groups = "drop")

Annoations_percentage_6feature$Annotation <- factor(Annoations_percentage_6feature$Annotation, levels=c("Unannotated","Intron" ,"Noncoding_other-exon","Upstream and 1st exon of noncoding","CDS and 3′-UTR of coding" ,"Upstream and 5′-UTR of coding" ))

# Reshape the data from wide to long format
Annoations_long2 <- Annoations_percentage_6feature %>%
  pivot_longer(
    cols = starts_with(c("X16c", "X2c", "X4c", "X8c", "Blc", "GV", "MII")),
    names_to = "Sample",
    values_to = "Percentage"
  ) %>%
  mutate(
    # Extract the stage from the sample column
    Stage = case_when(
      grepl("X16c", Sample) ~ "16-cell",
      grepl("X2c", Sample) ~ "2-cell",
      grepl("X4c", Sample) ~ "4-cell",
      grepl("X8c", Sample) ~ "8-cell",
      grepl("Blc", Sample) ~ "Blastocyst",
      grepl("GV", Sample) ~ "GV oocyte",
      grepl("MII", Sample) ~ "MII oocyte",
      TRUE ~ "Unknown"
    )
  )

Annoations_long2 <- Annoations_long2 %>%
  mutate(Stage = factor(Stage, levels = desired_order))

Annoations_long2 <- Annoations_long2 %>%
  mutate(
    Percentage = Percentage / 100  # Convert percentage values from 0-100 to 0-1
  )

# Fix annotation naming
Annoations_long2$Annotation <- as.factor(Annoations_long2$Annotation)
Annoations_long2$Annotation <- recode(Annoations_long2$Annotation, 
                                      "Noncoding_other-exon" = "Noncoding other exon")

svglite(filename = "Annotation_percentages.svg", width = 13, height = 10) 

# Create a combined boxplot and dot plot
ggplot(Annoations_long2, aes(x = Stage, y = Percentage)) +
  geom_jitter(aes(color = Annotation), shape = 19, width = 0.2, size = 4) +  # Dot plot
  geom_boxplot(aes(fill = Annotation), outlier.shape = NA, color = "black", width = 0.5, alpha = 0.3) +  # Boxplot
  labs(x = "Stage", y = "Aligned reads") +
  scale_y_continuous(labels = scales::percent) +  # Format y-axis as percentage
  scale_color_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +  # Apply custom colors to dots
  scale_fill_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +  # Apply custom colors to boxplots
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
    axis.text.y = element_text(size = 18),
    legend.position = "bottom",  # Adjust legend position
    #strip.text = element_text(size = 18),  # Adjust facet label size if needed
    strip.text = element_blank(),
    legend.title = element_text(size = 18),  # Increase legend title size
    legend.text = element_text(size = 18),  # Increase legend elements font size
    panel.spacing = unit(1, "lines"),  # Increase space between facets
    panel.border = element_rect(color = "grey", fill = NA, size = 0.5)  # Add grey border around each facet
  ) +
  facet_wrap(~ Annotation, scales = "free_y", ncol = 2)  # Wrap annotations into facets

dev.off()



################################################################################
### 3. RNA Content Plot: Normalized Poly(A)-tailed RNA Content Across Stages
################################################################################


TFEclass <- read.delim("TFEclass.txt", header=FALSE)
Map_and_Spike_in_5_end  <- read.delim("BOV_LIB_1.3-QC.txt", header=TRUE)

# Select and transform the data
Map_and_Spike_in_5_end  <- Map_and_Spike_in_5_end %>%
  select(column1 = 1, column5 = 5, column8 = 8) %>% # Select the 1st, 5th, and 8th columns
  mutate(column1 = substr(column1, 13, nchar(column1)),
         column1 = as.numeric(column1)) # Modify column1 to start from the 13th character

colnames(TFEclass) <- c("Sample","Stage")
colnames(Map_and_Spike_in_5_end) <- c("Sample","Mapped reads","Mapped 5-end spike-in reads")

Samples_and_Mapped_Reads<- left_join(Map_and_Spike_in_5_end, TFEclass, by = "Sample")

#Remove samples 11, 13, 26, 37, 23 and controls
Samples_and_Mapped_Reads <-Samples_and_Mapped_Reads[-c(1,3,5,19,36,31,16),]

Samples_and_Mapped_Reads <- Samples_and_Mapped_Reads %>%
  mutate(Relative_Content = `Mapped reads` / `Mapped 5-end spike-in reads`)

data <- Samples_and_Mapped_Reads


#GV_MII
# Compute averages for GV and MII
averages_GV_MII <- data %>%
  filter(Stage %in% c("GV", "MII")) %>%
  group_by(Stage) %>%
  summarize(Average_Relative_Content = mean(Relative_Content, na.rm = TRUE))

# Normalize data such that GV is 1 on the Y-axis
normalized_data_GV_MII <- data %>%
  filter(Stage %in% c("GV", "MII")) %>%
  left_join(averages_GV_MII, by = c("Stage" = "Stage")) %>%
  mutate(normalized=Average_Relative_Content/as.data.frame(averages_GV_MII)[1,2],normalized_relative_content=Relative_Content/as.data.frame(averages_GV_MII)[1,2])


#MII_2c
averages_MII_2c <- data %>%
  filter(Stage %in% c("2c", "MII")) %>%
  group_by(Stage) %>%
  summarize(Average_Relative_Content = mean(Relative_Content, na.rm = TRUE))

normalized_data_MII_2c <- data %>%
  filter(Stage %in% c("2c", "MII")) %>%
  left_join(averages_MII_2c, by = c("Stage" = "Stage")) %>%
  mutate(normalized=Average_Relative_Content/as.data.frame(averages_MII_2c)[2,2],normalized_relative_content=Relative_Content/as.data.frame(averages_MII_2c)[2,2])


#2c_4c
averages_2c_4c <- data %>%
  filter(Stage %in% c("2c", "4c")) %>%
  group_by(Stage) %>%
  summarize(Average_Relative_Content = mean(Relative_Content, na.rm = TRUE))

normalized_data_2c_4c <- data %>%
  filter(Stage %in% c("2c", "4c")) %>%
  left_join(averages_2c_4c, by = c("Stage" = "Stage")) %>%
  mutate(normalized=Average_Relative_Content/as.data.frame(averages_2c_4c)[1,2],normalized_relative_content=Relative_Content/as.data.frame(averages_2c_4c)[1,2])


#4c_8c
averages_4c_8c <- data %>%
  filter(Stage %in% c("8c", "4c")) %>%
  group_by(Stage) %>%
  summarize(Average_Relative_Content = mean(Relative_Content, na.rm = TRUE))

normalized_data_4c_8c <- data %>%
  filter(Stage %in% c("8c", "4c")) %>%
  left_join(averages_4c_8c, by = c("Stage" = "Stage")) %>%
  mutate(normalized=Average_Relative_Content/as.data.frame(averages_4c_8c)[1,2],normalized_relative_content=Relative_Content/as.data.frame(averages_4c_8c)[1,2])


#8c_16c
averages_8c_16c <- data %>%
  filter(Stage %in% c("8c", "16c")) %>%
  group_by(Stage) %>%
  summarize(Average_Relative_Content = mean(Relative_Content, na.rm = TRUE))

normalized_data_8c_16c <- data %>%
  filter(Stage %in% c("8c", "16c")) %>%
  left_join(averages_8c_16c, by = c("Stage" = "Stage")) %>%
  mutate(normalized=Average_Relative_Content/as.data.frame(averages_8c_16c)[2,2],normalized_relative_content=Relative_Content/as.data.frame(averages_8c_16c)[2,2])


#16c_Blc
averages_16c_Blc <- data %>%
  filter(Stage %in% c("Blc", "16c")) %>%
  group_by(Stage) %>%
  summarize(Average_Relative_Content = mean(Relative_Content, na.rm = TRUE))

normalized_data_16c_Blc <- data %>%
  filter(Stage %in% c("Blc", "16c")) %>%
  left_join(averages_16c_Blc, by = c("Stage" = "Stage")) %>%
  mutate(normalized=Average_Relative_Content/as.data.frame(averages_16c_Blc)[1,2],normalized_relative_content=Relative_Content/as.data.frame(averages_16c_Blc)[1,2])


library(ggpubfigs)
library(viridis)

# Define or extract the color palette
colors <- friendly_pal("wong_eight")

# Define the color mapping for each condition
condition_colors <- c(
  "GV" = colors[1],
  "MII" = colors[2],
  "2c" = colors[3],
  "4c" = colors[4],
  "8c" = colors[5],
  "16c" = colors[6],
  "Blastocyst" = colors[7]
)

library(dplyr)
Samples_and_Mapped_Reads_ref <- Samples_and_Mapped_Reads %>%
  mutate(GV_avg_ref = mean(Relative_Content[Stage == "GV"], na.rm = TRUE)) %>%
  mutate(as_comp_Gv_avg = Relative_Content/GV_avg_ref)

Samples_and_Mapped_Reads_avg <- Samples_and_Mapped_Reads %>%
  mutate(GV_avg_ref = mean(Relative_Content[Stage == "GV"], na.rm = TRUE)) %>%
  select(Stage,Relative_Content,GV_avg_ref)%>%
  group_by(Stage) %>%
  mutate(Stage_avg_ref = mean(Relative_Content, na.rm = TRUE)) %>%
  ungroup()%>%
  distinct(Stage,Stage_avg_ref,GV_avg_ref)%>%
  mutate(Relative_avg = Stage_avg_ref/GV_avg_ref)


# Update the stage levels and labels
Samples_and_Mapped_Reads_ref$Stage <- factor(
  Samples_and_Mapped_Reads_ref$Stage,
  levels = c("GV", "MII", "2c", "4c", "8c", "16c", "Blc"), # Correct levels
  labels = c("GV oocyte", "MII oocyte", "2-cell", "4-cell", "8-cell", "16-cell", "Blastocyst")
)

Samples_and_Mapped_Reads_avg$Stage <- factor(
  Samples_and_Mapped_Reads_avg$Stage,
  levels = c("GV", "MII", "2c", "4c", "8c", "16c", "Blc"), # Correct levels
  labels = c("GV oocyte", "MII oocyte", "2-cell", "4-cell", "8-cell", "16-cell", "Blastocyst")
)

svg(filename = "RNA_content_ref_GV.svg", width = 14, height = 7)

# Plotting
ggplot(Samples_and_Mapped_Reads_ref, aes(x = Stage)) +
  # Add jitter for individual data points
  geom_jitter(aes(y = as_comp_Gv_avg, color = Stage), 
              shape = 19, width = 0.05, size = 6, alpha = 0.8) +
  # Add boxplot for averages from Samples_and_Mapped_Reads_avg
  geom_boxplot(data = Samples_and_Mapped_Reads_avg, 
               aes(x = Stage, y = Relative_avg), 
               outlier.shape = NA, fill = NA, color = "red2", 
               width = 0.2, size = 1) +
  # Labels and scaling
  labs(x = "Stage", y = "Normalized Poly(A)-tailed RNA Content") +
  scale_y_continuous(limits = c(1 / 16, 1.5)) + # Adjust y-axis limits
  theme_minimal() +
  # Customize the x-axis and other themes
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 23, face = "bold", margin = margin(r = 20)),
    axis.text.x = element_text(size = 23, face = "bold"),
    axis.text.y = element_text(size = 20),
    legend.position = "none",
    panel.border = element_rect(color = "gray", fill = NA, size = 1)
  ) +
  scale_color_manual(values = condition_colors)

dev.off()
