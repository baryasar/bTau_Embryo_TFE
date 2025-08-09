library(DESeq2)
library(tidyverse)
library(pheatmap)
library(MASS)
library(ggpubfigs)
library(viridis)

TFEcounts <- read.delim("BOV_LIB_1.3_byTFE-counts_annotation.txt")
rownames(TFEcounts) <- TFEcounts[,1]

TFEcounts <- TFEcounts[,-c(1,2,3,4,5,6,7,8,9)]

#Remove sample 37, 11, 13, 26, 23
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

#Top Variable Genes by CV² filtered by adjusted p-value
nreads <- as.matrix(TFEcounts)

# Calculate row-wise coefficient of variation squared (CV²)
reads.rowVars <- apply(nreads, 1, var)
reads.rowMeans <- rowMeans(nreads)
y <- reads.rowVars / (reads.rowMeans^2)

# Normalize mean expression
x <- 1 / rowMeans(nreads)

# Identify spike-in controls
is.guide <- (substr(rownames(nreads), 1, 10) == "RNA_SPIKE_") == TRUE
x.guide <- x[is.guide]
y.guide <- y[is.guide]

# Fit a Gamma GLM to the guide data
fit <- glm(y.guide ~ x.guide, family = Gamma(link = 'identity'))
shape <- gamma.shape(fit)$alpha
scale <- (fit$coefficients[1] + fit$coefficients[2] / rowMeans(nreads)) / shape

# Adjust p-values
p.adj <- p.adjust(pgamma(q = y, lower.tail = FALSE, shape = shape, scale = scale), method = 'fdr')

# Create a data frame with gene names and their variability
gene_variability <- data.frame(
  gene = rownames(nreads),
  CV2 = y,
  p_adj = p.adj
)

vst_matrix <- as.data.frame(vst_matrix) %>%
  rownames_to_column(var= "gene")

vst_matrix_top <- gene_variability %>%
  arrange(desc(CV2)) %>%
  filter(p_adj < 0.05) %>%  # Filter for padj < 0.05
  slice_head(n=20000)%>%
  left_join(vst_matrix, by = c("gene" = "gene")) %>%
  dplyr::select(-CV2, -p_adj) %>%
  column_to_rownames(var = "gene")

vst_cor <- cor(vst_matrix_top)

svg(filename = "Heatmap_CV2_padj_top20000TFEs.svg", width = 8, height = 8) 

# Plot with pheatmap
pheatmap(
  vst_cor,
  annotation_col = meta,
  annotation_colors = annotation_colors,color = mako(50),
  fontsize_row = 9,  # Adjust font size for row labels
  fontsize_col = 9   # Adjust font size for column labels
)

dev.off()

#Number of TFEs p< 0.05
gene_variability %>%
  filter(p_adj < 0.05) %>%
  nrow()
#58737

#Count the number of annotated and unannotated
BOV_LIB_1.3_annotation <- read.delim("BOV_LIB_1.3_annotation.txt", header=FALSE)
gene_variability %>%
  arrange(desc(CV2)) %>%
  filter(p_adj < 0.05) %>%  # Filter for padj < 0.05
  slice_head(n = 20000) %>%
  left_join(BOV_LIB_1.3_annotation, by = c("gene" = "V1")) %>%
  summarise(
    Unannotated_count = sum(V3 == "Unannotated"),
    Annotated_count = sum(V3 != "Unannotated")
  )

TFE_corresponding_genes <- read.delim("TFE_corresponding_genes.txt", header=FALSE)

gene_variability %>%
  arrange(desc(CV2)) %>%
  filter(p_adj < 0.05) %>%  # Filter for padj < 0.05
  slice_head(n = 20000) %>%
  left_join(BOV_LIB_1.3_annotation, by = c("gene" = "V1")) %>%
  filter(V3 != "Unannotated") %>%
  left_join(TFE_corresponding_genes, by = c("gene" = "V1")) %>%
  filter(!is.na(V2.y)) %>%  # Remove rows where V2.y is NA
  mutate(
    V2_part1 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 1], V2.y),
    V2_part2 = ifelse(str_detect(V2.y, ":"), str_split(V2.y, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene)%>%
  nrow()
#6695

#Number of Known genes
TFE_corresponding_genes %>%
  filter(!is.na(V2)) %>%  # Ensure V2 is not NA
  mutate(
    V2_part1 = ifelse(str_detect(V2, ":"), str_split(V2, ":", simplify = TRUE)[, 1], V2),
    V2_part2 = ifelse(str_detect(V2, ":"), str_split(V2, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%
  nrow()
#18501
