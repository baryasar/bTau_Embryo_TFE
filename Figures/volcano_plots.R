library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggpubfigs)
library(viridis)
library(svglite)
library(tibble)


TFEcounts <- read.delim("../BOV_LIB_1.3_byTFE-counts_annotation.txt")
rownames(TFEcounts) <- TFEcounts[,1]

TFEcounts <- TFEcounts[,-c(1,2,3,4,5,6,7,8,9)]

#Remove sample 37 (16-cell), 11, 13, 26 (8-cell), 23 (2-cell)
TFEcounts <- TFEcounts[,-c(5,35,37,20,12)]


condition <- c("16c","16c","16c","16c","16c","16c","16c","2c","2c","2c","2c","4c","4c","4c","4c","4c","8c","8c","8c","8c","8c","8c","8c","Blastocyst","Blastocyst","Blastocyst","Blastocyst","Blastocyst","Blastocyst","Blastocyst","GV","GV","GV","GV","MII","MII","MII","MII","MII","MII","MII")

length(condition) == length(colnames(TFEcounts))

colnames(TFEcounts) <- c("16-cell rep1", "16-cell rep2", "16-cell rep3",
                         "16-cell rep4", "16-cell rep5",
                         "16-cell rep6", "16-cell rep7", "2-cell rep1", 
                         "2-cell rep2",  "2-cell rep3",  "2-cell rep4",
                         "4-cell rep1",  "4-cell rep2", 
                         "4-cell rep3",  "4-cell rep4",  "4-cell rep5", 
                         "8-cell rep1",  "8-cell rep2",  "8-cell rep3", 
                         "8-cell rep4",  "8-cell rep5",  "8-cell rep6", 
                         "8-cell rep7", "Blastocyst rep1", 
                         "Blastocyst rep2",  "Blastocyst rep3",  "Blastocyst rep4", 
                         "Blastocyst rep5",  "Blastocyst rep6",  "Blastocyst rep7", 
                         "GV oocyte rep1",   "GV oocyte rep2",   "GV oocyte rep3",  
                         "GV oocyte rep4",  
                         "MII oocyte rep1",   "MII oocyte rep2",   "MII oocyte rep3",  
                         "MII oocyte rep4",   "MII oocyte rep5",   "MII oocyte rep6",  
                         "MII oocyte rep7")


sample_names <- colnames(TFEcounts)

meta <- data.frame(condition,row.names = sample_names)

#checking if column and row names are in the same order
all(colnames(TFEcounts) == rownames(meta))

# Convert condition  to factors
meta$condition <- factor(meta$condition)

#Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = TFEcounts,
                              colData = meta,
                              design = ~ condition)

dds$condition <- factor(dds$condition, levels = c("GV","MII","2c","4c","8c","16c","Blastocyst"))

levels(dds$condition)

# Filter out genes with zero counts
dds <- dds[ rowSums(counts(dds)) > 0, ]

#extract the row names from the dds object to find where your control genes are located.
all_genes <- rownames(counts(dds))

# Find indices of control genes
ctrlGeneIndices <- grep("^RNA_SPIKE_", all_genes)

# Determine the size factors to use for normalization
dds <- estimateSizeFactors(dds, controlGenes = ctrlGeneIndices)

# See size factors for each column
sizeFactors(dds)

# Transform the normalized counts, variance stabilizing transformation (VST)
vst <- vst(dds, blind = TRUE)

# Extract the vst matrix from the object
vst_matrix <- assay(vst)

# Define or extract the color palette
colors <- friendly_pal("wong_eight")

# Define the color mapping for each condition
condition_colors <- c(
  "GV oocyte" = colors[1],
  "MII oocyte" = colors[2],
  "2-cell" = colors[3],
  "4-cell" = colors[4],
  "8-cell" = colors[5],
  "16-cell" = colors[6],
  "Blastocyst" = colors[7]
)

colnames(meta) <- c("Stage")
meta$Stage <- c("16-cell", "16-cell", "16-cell",
                "16-cell", "16-cell", "16-cell",
                "16-cell", "2-cell", 
                "2-cell",  "2-cell",
                "2-cell",   "4-cell",  "4-cell", 
                "4-cell",  "4-cell",  "4-cell", 
                "8-cell",  "8-cell",  "8-cell", 
                "8-cell",  "8-cell",  "8-cell", 
                "8-cell", "Blastocyst", 
                "Blastocyst",  "Blastocyst",  "Blastocyst", 
                "Blastocyst",  "Blastocyst",  "Blastocyst", 
                "GV oocyte",   "GV oocyte",   "GV oocyte",  
                "GV oocyte",  
                "MII oocyte",   "MII oocyte",   "MII oocyte",  
                "MII oocyte",   "MII oocyte",   "MII oocyte",  
                "MII oocyte")
# Create a list for annotation colors
annotation_colors <- list(
  Stage = condition_colors
)

# Ensure `meta` is properly formatted
meta$Stage <- factor(meta$Stage, levels = names(condition_colors))

#----------------------------DESEQ---------------------------------------

# Run the DESeq2 analysis
dds <- DESeq(dds)

resultsNames(dds)

plotDispEsts(dds)


res1 <- results(dds, contrast = c("condition", "MII", "GV"), alpha=0.05)
res2 <- results(dds, contrast = c("condition", "2c", "MII"), alpha=0.05)
res3 <- results(dds, contrast = c("condition", "4c", "2c"), alpha=0.05)
res4 <- results(dds, contrast = c("condition", "8c", "4c"), alpha=0.05)
res5 <- results(dds, contrast = c("condition", "16c", "8c"), alpha=0.05)
res6 <- results(dds, contrast = c("condition", "Blastocyst", "16c"), alpha=0.05)


#------------Convert TFE numbers to gene numbers-----------------

TFE_corresponding_genes <- read.delim("../TFE_corresponding_genes.txt", header=FALSE)
BOV_LIB_1.3_annotation <- read.delim("../BOV_LIB_1.3_annotation.txt", header=FALSE)

#RES1 UP
as.data.frame(res1) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#2733

#RES1 DOWN
as.data.frame(res1) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange < 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#3618

#RES2 UP
as.data.frame(res2) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#149

#RES2 DOWN
as.data.frame(res2) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange < 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#7400


#RES3 UP
as.data.frame(res3) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#21

#RES3 DOWN
as.data.frame(res3) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange < 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#87


#RES4 UP
as.data.frame(res4) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#1

#RES4 DOWN
as.data.frame(res4) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange < 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#70


#RES5 UP
as.data.frame(res5) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#731


#RES5 DOWN
as.data.frame(res5) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange < 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#8536


#RES6 UP
as.data.frame(res6) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#6501

#RES6 DOWN
as.data.frame(res6) %>%
  rownames_to_column(var = "TFE") %>%
  left_join(BOV_LIB_1.3_annotation, by = c("TFE" = "V1")) %>%
  filter(V3 != "Unannotated" & padj < 0.05 & log2FoldChange < 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#1301


# Get names for legend
down <- unlist(strsplit('Downregulated,Not Significant,Upregulated', split = ","))[1]
notsig <- unlist(strsplit('Downregulated,Not Significant,Upregulated', split = ","))[2]
up <- unlist(strsplit('Downregulated,Not Significant,Upregulated', split = ","))[3]


# Set colours
colours <- setNames(c("#3E4A89FF", "#707173FF", "#D8576BFF"), c(down, notsig, up))

# Convert the results table into a df

#MII/GV
res_df <- res1 %>%
  data.frame() %>%
  rownames_to_column(var="TFE")
#2C/MII
res_df <- res2 %>%
  data.frame() %>%
  rownames_to_column(var="TFE")
#4C/2C
res_df <- res3 %>%
  data.frame() %>%
  rownames_to_column(var="TFE")
#8C/4C
res_df <- res4 %>%
  data.frame() %>%
  rownames_to_column(var="TFE")
#16C/8C
res_df <- res5 %>%
  data.frame() %>%
  rownames_to_column(var="TFE")
#Blastocyst/16C
res_df <- res6 %>%
  data.frame() %>%
  rownames_to_column(var="TFE")

# Create columns from the column numbers specified
res_df <- res_df %>% mutate(fdr = .[[7]],
                            pvalue = .[[6]],
                            logfc = .[[3]])

# Create significant (sig) column
res_df <- mutate(res_df, sig = case_when(
  fdr < 0.05 & logfc > 0.0 ~ up,
  fdr < 0.05 & logfc < -0.0 ~ down,
  TRUE ~ notsig))

# NOTE: The object `res_df` is overwritten multiple times below to create volcano plots for different comparisons.
# To generate plots for different contrasts, assign the relevant results object (e.g., res1, res2, etc.) to `res_df` before plotting.

# Create the plot
p <- ggplot(data = res_df, aes(x = logfc, y = -log10(pvalue))) +
  geom_point(aes(colour = sig), size = 2, alpha = 0.6) +  # Adjust point size and transparency
  scale_color_manual(values = colours) +  # Apply your custom colors
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "MII oocyte vs GV oocyte"  ) +  # Clear axis labels
  theme_minimal(base_size = 14) +  # Use a clean, minimal theme and adjust base font size
  theme(
    panel.grid.major = element_blank(),  # Remove grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),  # Make axis lines bold
    axis.title = element_text(size = 18),  # Bold axis titles
    axis.text = element_text(size = 18),  # Larger axis text
    legend.position = "bottom",  # Move legend to top for better readability
    legend.key = element_blank(),  # Clean legend background
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 18),  # Adjust legend text size
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22)#,  # Optional title formatting
    #plot.margin = margin(1, 1, 1, 1, "cm")  # Adjust margins
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4)))# Legend points size


######### MII vs GV oocyte #########
tiff(filename="MII_vs_GV.tiff",width=3000, height=3000,res=300)
p + 
  annotate("text", x = Inf, y = Inf, label = "▲ 8902 TFEs/2736 genes", color = "#D8576BFF", size = 4.8, alpha = 0.8, hjust = 1, vjust = 1) + 
  annotate("text", x = -Inf, y = Inf, label = "▼ 6508 TFEs/3623 genes", color = "#3E4A89FF", size = 4.8, alpha = 0.8, hjust = 0, vjust = 1) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 60)) +
  theme(plot.margin = margin(1, 3, 1, 3, "cm")) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "MII oocyte vs GV oocyte")
dev.off()

######### 2-cell vs MII oocyte #########
tiff(filename="2c_vs_MII.tiff",width=3000, height=3000,res=300)
p + 
  annotate("text", x = Inf, y = Inf, label = "▲ 211 TFEs/150 genes", color = "#D8576BFF", size = 4.8, alpha = 0.8, hjust = 1, vjust = 1) +
  annotate("text", x = -Inf, y = Inf, label = "▼ 21918 TFEs/7426 genes", color = "#3E4A89FF", size = 4.8, alpha = 0.8, hjust = 0, vjust = 1) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 60)) +
  theme(plot.margin = margin(1, 3, 1, 3, "cm")) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "2-cell vs MII oocyte")
dev.off()

######### 4-cell vs 2-cell #########
tiff(filename="4c_vs_2c.tiff",width=3000, height=3000,res=300)
p + 
  annotate("text", x = Inf, y = Inf, label = "▲ 29 TFEs/21 genes", color = "#D8576BFF", size = 4.8, alpha = 0.8, hjust = 1, vjust = 1) +
  annotate("text", x = -Inf, y = Inf, label = "▼ 122 TFEs/88 genes", color = "#3E4A89FF", size = 4.8, alpha = 0.8, hjust = 0, vjust = 1) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 60)) +
  theme(plot.margin = margin(1, 3, 1, 3, "cm")) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "4-cell vs 2-cell")
dev.off()

######### 8-cell vs 4-cell #########
tiff(filename="8c_vs_4c.tiff",width=3000, height=3000,res=300)
p + 
  annotate("text", x = Inf, y = Inf, label = "▲ 4 TFEs/1 genes", color = "#D8576BFF", size = 4.8, alpha = 0.8, hjust = 1, vjust = 1) +
  annotate("text", x = -Inf, y = Inf, label = "▼ 92 TFEs/70 genes", color = "#3E4A89FF", size = 4.8, alpha = 0.8, hjust = 0, vjust = 1) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 60)) +
  theme(plot.margin = margin(1, 3, 1, 3, "cm")) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "8-cell vs 4-cell")
dev.off()

######### 16-cell vs 8-cell #########
tiff(filename="16c_vs_8c.tiff",width=3000, height=3000,res=300)
p + 
  annotate("text", x = Inf, y = Inf, label = "▲ 1621 TFEs/733 genes", color = "#D8576BFF", size = 4.8, alpha = 0.8, hjust = 1, vjust = 1) +
  annotate("text", x = -Inf, y = Inf, label = "▼ 20461 TFEs/8564 genes", color = "#3E4A89FF", size = 4.8, alpha = 0.8, hjust = 0, vjust = 1) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 60)) +
  theme(plot.margin = margin(1, 3, 1, 3, "cm")) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "16-cell vs 8-cell")
dev.off()

######### Blastocyst vs 16-cell #########
tiff(filename="Blc_vs_16c.tiff",width=3000, height=3000,res=300)
p + 
  annotate("text", x = Inf, y = Inf, label = "▲ 12989 TFEs/6510 genes", color = "#D8576BFF", size = 4.8, alpha = 0.8, hjust = 1, vjust = 1) +
  annotate("text", x = -Inf, y = Inf, label = "▼ 2057 TFEs/1302 genes", color = "#3E4A89FF", size = 4.8, alpha = 0.8, hjust = 0, vjust = 1) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 200)) +
  theme(plot.margin = margin(1, 3, 1, 3, "cm")) +
  labs(x = "Log2 Fold Change", y = "-Log10(P-value)",title = "Blastocyst vs 16-cell")
dev.off()
