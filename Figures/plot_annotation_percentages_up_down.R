library(RColorBrewer)
library(dplyr)
library(scales)
BOV_LIB_1.3_annotation <- read.delim("BOV_LIB_1.3_annotation.txt", header=FALSE)

# Convert row names of res1 to a column
res1_df <- as.data.frame(res1) %>%
  rownames_to_column(var = "TFE")

# Perform the left join with BOV_LIB_1.3_annotation
res1_df_annot <- res1_df %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1"))

res1_df_annot_perc_DOWN <- res1_df_annot %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

res1_df_annot_perc_UP <- res1_df_annot %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

#-----------------------MII vs GV-----------------------------

# Add a new column to indicate "UP" or "DOWN"
res1_df_annot_perc_DOWN <- res1_df_annot_perc_DOWN %>%
  mutate(direction = "Down-")

res1_df_annot_perc_UP <- res1_df_annot_perc_UP %>%
  mutate(direction = "Up-")

# Combine the two dataframes
combined_res1 <- bind_rows(res1_df_annot_perc_DOWN, res1_df_annot_perc_UP)

# Define a named vector for renaming
rename_annotations <- c(
  "Coding_3UTR+Coding_CDS" = "CDS and 3′-UTR of coding",
  "Coding_5UTR+Coding_upstream" = "Upstream and 5′-UTR of coding",
  "Intron" = "Intron",
  "Noncoding_1st-exon+Noncoding_upstream" = "Upstream and 1st exon of noncoding",
  "Noncoding_other-exon" = "Noncoding other exon",
  "Unannotated" = "Unannotated"
)

# Update the annotations in combined_res1
combined_res1 <- combined_res1 %>%
  mutate(V3 = recode(V3, !!!rename_annotations))

# Define the desired order of levels
desired_levels <- c(
  "Unannotated",
  "Intron",
  "Noncoding other exon",
  "Upstream and 1st exon of noncoding",
  "CDS and 3′-UTR of coding",
  "Upstream and 5′-UTR of coding"
)

# Apply the levels to the factor
combined_res1$V3 <- factor(combined_res1$V3, levels = desired_levels)
combined_res1$direction <- factor(combined_res1$direction, levels = c("Up-", "Down-"))


svglite(filename = "TFE_annotations_res1.svg", width = 15, height = 8)  

ggplot(combined_res1, aes(x = "", y = percentage, fill = V3)) +
  geom_bar(stat = "identity", width = 0.7, position = "stack") +
  coord_flip() +  # Flip to horizontal bars
  facet_wrap(~direction, ncol = 1) +  # Create separate plots for UP and DOWN
  scale_fill_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +  # Format y-axis as percentage with % sign
  labs(
    x = "TFE",
    y = "Percentage",
    fill = "Annotation",
    title = "MII oocyte vs GV oocyte"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 22, hjust = 0, face = "bold"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.text = element_text(size = 22)  # Increase legend elements font size
  )

dev.off()

#-----------------------2-cell vs MII-----------------------------
# Convert row names of res2 to a column
res2_df <- as.data.frame(res2) %>%
  rownames_to_column(var = "TFE")

# Perform the left join with BOV_LIB_1.3_annotation
res2_df_annot <- res2_df %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1"))

res2_df_annot_perc_DOWN <- res2_df_annot %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

res2_df_annot_perc_UP <- res2_df_annot %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

# Add a new column to indicate "UP" or "DOWN"
res2_df_annot_perc_DOWN <- res2_df_annot_perc_DOWN %>%
  mutate(direction = "Down-")

res2_df_annot_perc_UP <- res2_df_annot_perc_UP %>%
  mutate(direction = "Up-")

# Combine the two dataframes
combined_res2 <- bind_rows(res2_df_annot_perc_DOWN, res2_df_annot_perc_UP)

# Define a named vector for renaming
rename_annotations <- c(
  "Coding_3UTR+Coding_CDS" = "CDS and 3′-UTR of coding",
  "Coding_5UTR+Coding_upstream" = "Upstream and 5′-UTR of coding",
  "Intron" = "Intron",
  "Noncoding_1st-exon+Noncoding_upstream" = "Upstream and 1st exon of noncoding",
  "Noncoding_other-exon" = "Noncoding other exon",
  "Unannotated" = "Unannotated"
)

# Update the annotations in combined_res2
combined_res2 <- combined_res2 %>%
  mutate(V3 = recode(V3, !!!rename_annotations))

# Define the desired order of levels
desired_levels <- c(
  "Unannotated",
  "Intron",
  "Noncoding other exon",
  "Upstream and 1st exon of noncoding",
  "CDS and 3′-UTR of coding",
  "Upstream and 5′-UTR of coding"
)

# Apply the levels to the factor
combined_res2$V3 <- factor(combined_res2$V3, levels = desired_levels)

combined_res2$direction <- factor(combined_res2$direction, levels = c("Up-", "Down-"))


svglite(filename = "TFE_annotations_res2.svg", width = 15, height = 8)

ggplot(combined_res2, aes(x = "", y = percentage, fill = V3)) +
  geom_bar(stat = "identity", width = 0.7, position = "stack") +
  coord_flip() +  # Flip to horizontal bars
  facet_wrap(~direction, ncol = 1) +  # Create separate plots for UP and DOWN
  scale_fill_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +  # Format y-axis as percentage with % sign
  labs(
    x = "TFE",
    y = "Percentage",
    fill = "Annotation",
    title = "2-cell vs MII oocyte"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 22, hjust = 0, face = "bold"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.text = element_text(size = 22)  # Increase legend elements font size
  )

dev.off()

#-----------------------16-cell vs 8-cell-----------------------------
# Convert row names of res5 to a column
res5_df <- as.data.frame(res5) %>%
  rownames_to_column(var = "TFE")

# Perform the left join with BOV_LIB_1.3_annotation
res5_df_annot <- res5_df %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1"))

res5_df_annot_perc_DOWN <- res5_df_annot %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

res5_df_annot_perc_UP <- res5_df_annot %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

# Add a new column to indicate "UP" or "DOWN"
res5_df_annot_perc_DOWN <- res5_df_annot_perc_DOWN %>%
  mutate(direction = "Down-")

res5_df_annot_perc_UP <- res5_df_annot_perc_UP %>%
  mutate(direction = "Up-")

# Combine the two dataframes
combined_res5 <- bind_rows(res5_df_annot_perc_DOWN, res5_df_annot_perc_UP)

# Define a named vector for renaming
rename_annotations <- c(
  "Coding_3UTR+Coding_CDS" = "CDS and 3′-UTR of coding",
  "Coding_5UTR+Coding_upstream" = "Upstream and 5′-UTR of coding",
  "Intron" = "Intron",
  "Noncoding_1st-exon+Noncoding_upstream" = "Upstream and 1st exon of noncoding",
  "Noncoding_other-exon" = "Noncoding other exon",
  "Unannotated" = "Unannotated"
)

# Update the annotations in combined_res5
combined_res5 <- combined_res5 %>%
  mutate(V3 = recode(V3, !!!rename_annotations))

# Define the desired order of levels
desired_levels <- c(
  "Unannotated",
  "Intron",
  "Noncoding other exon",
  "Upstream and 1st exon of noncoding",
  "CDS and 3′-UTR of coding",
  "Upstream and 5′-UTR of coding"
)

# Apply the levels to the factor
combined_res5$V3 <- factor(combined_res5$V3, levels = desired_levels)
combined_res5$direction <- factor(combined_res5$direction, levels = c("Up-", "Down-"))


svglite(filename = "TFE_annotations_res5.svg", width = 15, height = 8)

ggplot(combined_res5, aes(x = "", y = percentage, fill = V3)) +
  geom_bar(stat = "identity", width = 0.7, position = "stack") +
  coord_flip() +  # Flip to horizontal bars
  facet_wrap(~direction, ncol = 1) +  # Create separate plots for UP and DOWN
  scale_fill_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +  # Format y-axis as percentage with % sign
  labs(
    x = "TFE",
    y = "Percentage",
    fill = "Annotation",
    title = "16-cell vs 8-cell"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 22, hjust = 0, face = "bold"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.text = element_text(size = 22)  # Increase legend elements font size
  )

dev.off()

#-----------------------Blastocyst vs 16-cell-----------------------------
# Convert row names of res6 to a column
res6_df <- as.data.frame(res6) %>%
  rownames_to_column(var = "TFE")

# Perform the left join with BOV_LIB_1.3_annotation
res6_df_annot <- res6_df %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1"))

res6_df_annot_perc_DOWN <- res6_df_annot %>%
  filter(padj < 0.05 & log2FoldChange < 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

res6_df_annot_perc_UP <- res6_df_annot %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  count(V3) %>%
  mutate(V3 = case_when(
    V3 %in% c("Coding_5UTR", "Coding_upstream") ~ "Coding_5UTR+Coding_upstream",
    V3 %in% c("Coding_3UTR", "Coding_CDS") ~ "Coding_3UTR+Coding_CDS",
    V3 %in% c("Noncoding_1st-exon", "Noncoding_upstream") ~ "Noncoding_1st-exon+Noncoding_upstream",
    TRUE ~ V3
  )) %>%
  group_by(V3) %>%
  summarize(n = sum(n))%>%
  mutate(percentage = (n / sum(n)) * 100)

# Add a new column to indicate "UP" or "DOWN"
res6_df_annot_perc_DOWN <- res6_df_annot_perc_DOWN %>%
  mutate(direction = "Down-")

res6_df_annot_perc_UP <- res6_df_annot_perc_UP %>%
  mutate(direction = "Up-")

# Combine the two dataframes
combined_res6 <- bind_rows(res6_df_annot_perc_DOWN, res6_df_annot_perc_UP)

# Define a named vector for renaming
rename_annotations <- c(
  "Coding_3UTR+Coding_CDS" = "CDS and 3′-UTR of coding",
  "Coding_5UTR+Coding_upstream" = "Upstream and 5′-UTR of coding",
  "Intron" = "Intron",
  "Noncoding_1st-exon+Noncoding_upstream" = "Upstream and 1st exon of noncoding",
  "Noncoding_other-exon" = "Noncoding other exon",
  "Unannotated" = "Unannotated"
)

# Update the annotations in combined_res6
combined_res6 <- combined_res6 %>%
  mutate(V3 = recode(V3, !!!rename_annotations))

# Define the desired order of levels
desired_levels <- c(
  "Unannotated",
  "Intron",
  "Noncoding other exon",
  "Upstream and 1st exon of noncoding",
  "CDS and 3′-UTR of coding",
  "Upstream and 5′-UTR of coding"
)

# Apply the levels to the factor
combined_res6$V3 <- factor(combined_res6$V3, levels = desired_levels)
combined_res6$direction <- factor(combined_res6$direction, levels = c("Up-", "Down-"))

svglite(filename = "TFE_annotations_res6.svg", width = 15, height = 8)

ggplot(combined_res6, aes(x = "", y = percentage, fill = V3)) +
  geom_bar(stat = "identity", width = 0.7, position = "stack") +
  coord_flip() +  # Flip to horizontal bars
  facet_wrap(~direction, ncol = 1) +  # Create separate plots for UP and DOWN
  scale_fill_manual(values = brewer.pal(n = 6, name = "YlGnBu")) +
  scale_y_continuous(labels = label_percent(scale = 1)) +  # Format y-axis as percentage with % sign
  labs(
    x = "TFE",
    y = "Percentage",
    fill = "Annotation",
    title = "Blastocyst vs 16-cell"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 22),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 22, hjust = 0, face = "bold"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 22),  # Increase legend title size
    legend.text = element_text(size = 22)  # Increase legend elements font size
  )

dev.off()
