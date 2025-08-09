library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(viridis)
library(ggpubfigs)
library(tibble)
library(svglite)

# NOTE: Some of the objects were pre-loaded in the environment. These objects were loaded or created in other scripts within this project.

# Define or extract the color palette
colors <- friendly_pal("wong_eight")

TFE_corresponding_genes <- read.delim("TFE_corresponding_genes.txt", header=FALSE)

#Convert the extracted data to a data frame for easier manipulation
normalized_data_df <- as.data.frame(as.matrix(normalized_data_filtered))

# Check the column names to ensure they correspond to the sample identifiers
sample_names <- colnames(normalized_data_df)

# Get the developmental stage information from the Seurat object
stage_info <- strt.seurat.obj$orig.ident[sample_names]

stage_info_vector <- as.vector(stage_info)

colnames(normalized_data_df) <- stage_info_vector

# Create unique column names by appending suffixes
colnames(normalized_data_df) <- make.unique(colnames(normalized_data_df), sep = "_")

# Convert to long format
normalized_long_df <- normalized_data_df %>%
  rownames_to_column(var = "gene") %>%  # Keep gene names as a column
  pivot_longer(-gene, names_to = "stage", values_to = "expression") # Gather data

# Remove suffixes from stage names
normalized_long_df <- normalized_long_df %>%
  mutate(stage = sub("_\\d+$", "", stage))  # Remove the suffix

median_expression_df <- normalized_long_df %>%
  group_by(gene, stage) %>%
  summarize(average_expression = median(expression, na.rm = TRUE)) %>%
  ungroup()

# Pivot the dataframe to wide format
wide_median_expression_df <- median_expression_df %>%
  pivot_wider(names_from = stage, values_from = average_expression)

colnames(wide_median_expression_df) <- c("gene","16-cell","2-cell","4-cell","8-cell","Blastocyst","GV oocyte","MII oocyte")

#-------------Only MII expressed-----------------------
res1_df <- as.data.frame(res1) %>%
  rownames_to_column(var = "TFE")

res1_df_UPgenes <- res1_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE" = "V1")) %>%
  replace_na(list(V2 = "Unannot.")) %>%
  mutate(TFE_V2 = paste(TFE, V2, sep = "/"))

only_MII_expressed <- wide_median_expression_df %>%
  filter(
    `GV oocyte` < -10 & 
      `MII oocyte` > -10 & 
      `2-cell` < -10 & 
      `4-cell` < -10 & 
      `8-cell` < -10 &
      `16-cell` < -10 &
      `Blastocyst` < -10 
  ) %>%
  dplyr::select(gene) %>%
  inner_join(res1_df_UPgenes, by = c("gene" = "TFE")) %>%
  filter(log2FoldChange > 3)%>%
  dplyr::select(gene,V2)

#------------Only 16-cell expressed----------------------
res5_df <- as.data.frame(res5) %>%
  rownames_to_column(var = "TFE")

res5_df_UPgenes <- res5_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE"= "V1"))%>%
  replace_na(list(V2 = "Unannot.")) %>%
  mutate(TFE_V2 = paste(TFE, V2, sep = "/"))


#Only 16-cell expressed
only_16c_expressed <- wide_median_expression_df %>%
  filter(
    `GV oocyte` < -10 & 
      `MII oocyte` < -10 & 
      `2-cell` < -10 & 
      `4-cell` < -10 & 
      `8-cell` < -10 &
      `16-cell` > -10 &
      `Blastocyst` < -10 
  ) %>%
  dplyr::select(gene) %>%
  inner_join(res5_df_UPgenes, by = c("gene" = "TFE")) %>%
  filter(log2FoldChange > 3)%>%
  dplyr::select(gene,V2)

#------------Only Blastocyst expressed-----------------------
res6_df <- as.data.frame(res6) %>%
  rownames_to_column(var = "TFE")

res6_df_UPgenes <- res6_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE"= "V1"))%>%
  replace_na(list(V2 = "Unannot.")) %>%
  mutate(TFE_V2 = paste(TFE, V2, sep = "/"))

#Only Blastocyst expressed
only_blc_expressed <- wide_median_expression_df %>%
  filter(
    `GV oocyte` < -10 & 
      `MII oocyte` < -10 & 
      `2-cell` < -10 & 
      `4-cell` < -10 & 
      `8-cell` < -10 &
      `16-cell` < -10 &
      `Blastocyst` > -10 
  ) %>%
  dplyr::select(gene) %>%
  inner_join(res6_df_UPgenes, by = c("gene" = "TFE")) %>%
  filter(log2FoldChange > 3)%>%
  dplyr::select(gene,V2)

#------------Only GV oocyte expressed-----------------------
res0 <- results(dds, contrast = c("condition", "GV", "MII"), alpha=0.05)

res0_df <- as.data.frame(res0) %>%
  rownames_to_column(var = "TFE")

res0_df_UPgenes <- res0_df %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>%
  left_join(TFE_corresponding_genes, by = c("TFE"= "V1"))%>%
  replace_na(list(V2 = "Unannot.")) %>%
  mutate(TFE_V2 = paste(TFE, V2, sep = "/"))

#Only GV expressed foldchange 4.5
only_GV_expressed <- wide_median_expression_df %>%
  filter(
    `GV oocyte` > -10 & 
      `MII oocyte` < -10 & 
      `2-cell` < -10 & 
      `4-cell` < -10 & 
      `8-cell` < -10 &
      `16-cell` < -10 &
      `Blastocyst` < -10 
  ) %>%
  dplyr::select(gene) %>%
  inner_join(res0_df_UPgenes, by = c("gene" = "TFE")) %>%
  filter(log2FoldChange > 3)%>%
  dplyr::select(gene,V2)


#COUNTING STAGE TOTAL NUMBER OF DEFINED GENES AND STAGE-SPECIFIC GENE NUMBERS AND STAGE_SPECIFIC UNANNOTATED NUMBERS

# Split V2 only for rows with a colon and gather all strings into a single column
TFE_corresponding_genes %>%
  mutate(
    V2_part1 = ifelse(str_detect(V2, ":"), str_split(V2, ":", simplify = TRUE)[, 1], V2),
    V2_part2 = ifelse(str_detect(V2, ":"), str_split(V2, ":", simplify = TRUE)[, 2], NA)
  ) %>%
  # Combine V2, V2_part1, and V2_part2 into one long vector of strings
  select(V2_part1, V2_part2) %>%
  pivot_longer(cols = everything(), values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(gene) %>%   # Keep only distinct gene names
  nrow()
#18501

only_GV_expressed %>%
  filter(grepl("Unannot\\.", V2)) %>%
  nrow()
#221

only_MII_expressed %>%
  filter(grepl("Unannot\\.", V2)) %>%
  nrow()
#572

only_16c_expressed %>%
  filter(grepl("Unannot\\.", V2)) %>%
  nrow()
#117

only_blc_expressed %>%
  filter(grepl("Unannot\\.", V2)) %>%
  nrow()
#609

#-1 is for removing Unannotated
only_GV_expressed %>%
  select(V2) %>% 
  arrange(V2) %>% 
  distinct(V2) %>% 
  nrow()-1
#339

only_MII_expressed %>%
  select(V2) %>% 
  arrange(V2) %>% 
  distinct(V2) %>% 
  nrow()-1
#155

only_16c_expressed %>%
  select(V2) %>% 
  arrange(V2) %>% 
  distinct(V2) %>% 
  nrow()-1
#75

only_blc_expressed %>%
  select(V2) %>% 
  arrange(V2) %>% 
  distinct(V2) %>% 
  nrow()-1
#1400

#----------------------------

Bos_taurus_TF <- read.delim("~/ARS-UCD1.3/Bos_taurus_TF.txt")
#replaces empty strings ("") in Symbol with the corresponding value from Ensembl.
Bos_taurus_TF <- Bos_taurus_TF %>%
  mutate(Symbol = ifelse(Symbol == "", Ensembl, Symbol))

# Define genes to annotate on the right
genes_to_show1 <- Bos_taurus_TF %>%
  inner_join(only_GV_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene) %>% 
  pull(gene)

genes_to_show2 <- Bos_taurus_TF %>%
  inner_join(only_MII_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)%>% 
  pull(gene)

genes_to_show3 <- Bos_taurus_TF %>%
  inner_join(only_16c_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)%>% 
  pull(gene)

genes_to_show4 <- Bos_taurus_TF %>%
  inner_join(only_blc_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)%>% 
  pull(gene)

Symbolgene1 <- Bos_taurus_TF %>%
  inner_join(only_GV_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)

Symbolgene2 <- Bos_taurus_TF %>%
  inner_join(only_MII_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)

Symbolgene3 <- Bos_taurus_TF %>%
  inner_join(only_16c_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)

Symbolgene4 <- Bos_taurus_TF %>%
  inner_join(only_blc_expressed, by = c("Symbol" = "V2")) %>% 
  dplyr::select(Symbol,gene)

Symbolgene <- rbind(Symbolgene1,Symbolgene2,Symbolgene3,Symbolgene4)

#------#-----#
#For plotting all stages
genes_to_show <- c(genes_to_show1,genes_to_show2,genes_to_show3,genes_to_show4)
#For plotting the 16-cell stage
genes_to_show <- genes_to_show3
#------#------#

genes_to_show_df <- data.frame(genes = genes_to_show)

genes_to_show_name <- genes_to_show_df %>% 
  inner_join(TFE_corresponding_genes, by = c("genes" = "V1")) %>% 
  pull(V2)

genes_to_show <- c(genes_to_show,"TFE43635","TFE45905","TFE1057")
genes_to_show_name <- c(genes_to_show_name,"LEUTX","DUXA","DPPA2")

##################################################
#HEATMAP WITH COMPLEXHEATMAP WITH LABELS
#TFEs
selected_genes <- c(only_GV_expressed$gene,only_MII_expressed$gene,only_16c_expressed$gene,only_blc_expressed$gene)
#genes
gene_labels <- c(only_GV_expressed$V2,only_MII_expressed$V2,only_16c_expressed$V2,only_blc_expressed$V2)


# Unique is used to eliminate duplicated genes
selected_genes <- unique(selected_genes)

# Extract the expression data for the selected genes from the scaled layer
expr_matrix <- as.matrix(
  GetAssayData(
    strt.seurat.obj, 
    layer = "scale.data", 
    assay = "RNA"
  )[
    selected_genes,
  ]
)

# Prepare cell type annotations
stages <- Idents(strt.seurat.obj)
# Create vector with unique and sorted cell types
unique_stages <- sort(unique(stages))

stages_colors <- setNames(
  friendly_pal("wong_eight")[1:length(unique_stages)],  # Index to get the right number of colors
  unique_stages
)

# Values greater than 2.5 become 2.5.
# Values less than -2.5 become -2.5.
# Values between -2.5 and 2.5 remain unchanged.
expr_matrix_clipped <- pmax(pmin(expr_matrix, 2.5), -2.5)

palette_expression_level <- colorRamp2(
  breaks = seq(min(expr_matrix_clipped), max(expr_matrix_clipped), length.out = 256),  # Create a sequence of breakpoints covering your data range
  colors = viridis(256, option = "magma")  # Use the full viridis magma palette
)

mark_at <- which(rownames(expr_matrix_clipped) %in% genes_to_show)
ha <- rowAnnotation(
  foo = anno_mark(
    at = mark_at,
    labels = genes_to_show_name
  )
)

# Find the indices of the genes in the heatmap matrix
mark_at <- which(rownames(expr_matrix_clipped) %in% genes_to_show)

# Match the correct order of names in expr_matrix
mark_labels <- genes_to_show_name[match(rownames(expr_matrix_clipped)[mark_at], genes_to_show)]

# Create the annotation with correct labels
ha <- rowAnnotation(
  foo = anno_mark(
    at = mark_at,
    labels = mark_labels,
    labels_gp = gpar(fontface = "italic") 
  )
)

library(svglite)
# Plotting the heatmap with all stages marked with TFs
svglite(filename = "heatmap_labeled.svg", width = 20, height = 24) 
# Plotting the heatmap with the 16-cel stage marked with TFs
svglite(filename = "heatmap_labeled16c.svg", width = 8, height = 10)

Heatmap(
  matrix = expr_matrix_clipped,
  row_order = selected_genes,
  name = "Scaled expression",
  column_split = factor(stages, levels = unique_stages),
  # Do not cluster rows or columns otherwise row_order/column_order are ignored
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title_rot = 45,
  show_column_names = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  cluster_column_slices = TRUE,
  column_gap = unit(1.5, "mm"),
  row_names_gp = gpar(fontsize = 5),
  column_title = unique_stages,
  column_title_gp = gpar(fontsize = 15),
  # Add the annotation bars on top
  top_annotation = HeatmapAnnotation(
    cluster = stages,
    # Color palettes
    col = list(cluster = stages_colors),
    show_legend = FALSE,
    gap = unit(1, "mm"),
    simple_anno_size = unit(2, "mm"),
    annotation_name_gp =  gpar(fontsize = 1, col = "white")
  ),
  right_annotation = ha,
  # Color scale for expression
  col = palette_expression_level,
  use_raster = FALSE,
  heatmap_legend_param = list(
    title = "Scaled\nexpression",
    at = c(-2, -1, 0, 1, 2),  # Explicitly setting the legend breaks to match your data range
    labels = c("-2","-1", "0","1", "2")  # Labels for the legend at the set breaks
  )
)
dev.off()



#-----------Heatmap for selected PRDL genes-------------------

# Step 1: Define your target transcription factor aliases
target_genes <- toupper(c("Alx1","Alx3","Alx4", "Argfx", "Arx", "Dmbx", "Dprx", "Drgx", "Duxa", "Esx1", 
                          "Gsc", "Hesx1", "Hopx", "Isx", "Leutx", "Mixl1", "Nobox", "Otp", 
                          "Otx1","Otx2", "Phox2a","Phox2b", "Pitx1","Pitx2","Pitx3", "Prop", "Prrx2", "Rax", "Rhox", "Sebox", 
                          "Shox","Shox2", "Tprx", "Uncx", "Vsx1","Vsx2","Tprx1","Tprx2","Tprx3"))

# Step 2: Match those aliases to actual gene names via your mapping table
filtered_genes <- TFE_corresponding_genes %>%
  filter(V2 %in% target_genes) %>%
  distinct(V1, V2) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = " / "))

# Step 3: Extract the actual gene names to be plotted
selected_genes <- filtered_genes$V1

# Step 4: Build expression-based axis labels (italic for V2 part)
label_expressions <- mapply(function(gene, tf) {
  bquote(.(gene) ~ "/" ~ italic(.(tf)))
}, filtered_genes$V1, filtered_genes$V2, SIMPLIFY = FALSE)

# Step 5: Create a full named label vector for all selected genes (fallback to gene name if no match)
final_labels <- setNames(as.list(selected_genes), selected_genes)
final_labels[names(label_expressions)] <- label_expressions

# Step 6: Generate the heatmap
pdf("PRDL_heatmap.pdf", width = 9, height = 12)  

DoHeatmap(strt.seurat.obj, 
          features = selected_genes, 
          slot = "scale.data", 
          group.by = "ident", 
          label = TRUE, 
          draw.lines = TRUE, 
          group.colors = colors) + 
  scale_color_manual(values = colors) +
  guides(color = "none") +
  scale_y_discrete(labels = final_labels) +  # ← Custom label vector with italic TF names
  scale_fill_gradientn(name = "Scaled\nexpression",  # Set the legend title here
    colors = viridis::viridis(256, option = "magma"), na.value = "white") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  )

dev.off()
