library(patchwork)
library(dplyr)

# NOTE: Some of the objects were pre-loaded in the environment. These objects were loaded or created in other scripts within this project.

TFEpeaks <- read.delim("BOV_LIB_1.3_peaks.bed", header=FALSE)
TFEs_DeepTSS <- read.delim("TFEs_DeepTSS.txt", header=FALSE)

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
  select(gene) %>%
  inner_join(res5_df_UPgenes, by = c("gene" = "TFE")) %>%
  filter(log2FoldChange > 0)%>%
  dplyr::select(gene,V2)

only_16c <- TFEpeaks %>%
  semi_join(only_16c_expressed, by = c("V4" = "gene")) %>%
  semi_join(TFEs_DeepTSS,by = c("V4" = "V1"))

adjust_promoters <- function(df) {
  df <- df %>%
    mutate(
      promoter_start = ifelse(V6 == "+", V2 - 400, V3 - 100),
      promoter_end = ifelse(V6 == "+", V2 + 100, V3 + 400)
    ) %>%
    dplyr::select(V1, promoter_start, promoter_end, V4, V5, V6) %>%
    mutate(promoter_start = ifelse(promoter_start < 0, 0, promoter_start)) # Ensure no negative start positions
  return(df)
}

# extract promoters
only_16c_promoters <- adjust_promoters(only_16c)

TFE_corresponding_genes <- read.delim("TFE_corresponding_genes.txt", header=FALSE)

only_16c_promoters <- only_16c_promoters %>%
left_join(TFE_corresponding_genes, by = c("V4" = "V1")) %>%
  replace_na(list(V2 = "Unannot.")) %>%
  mutate(TFE_V2 = paste(V4, V2, sep = "_")) %>%
  dplyr::select(V1,promoter_start,promoter_end,TFE_V2,V5,V6)

write.table(only_16c_promoters, "promoter_only_16c.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# After running XSTREME I selected only motifs in the first 20 rank and with known binders

#MEME-14 rank 1
binders_1 <- c("DUXA", "PHOX2A", "PRRX1", "PHOX2B", 
       "PROP1", "CART1", "DRGX", "ALX4", "ISX", "UNCX", 
       "ALX3", "ALX1", "ARX")
binders_1 <- sort(binders_1)

#MEME-18 rank 3
binders_2 <- c("EHF", "ELF1", "GABPA", "SPI1", "IKZF3", "ELF5", "SPIB", 
       "ETV1", "ERG", "ELF3", "ZNF675", "ZIC2", "ZIC3")
binders_2 <- sort(binders_2)

#MEME-4 rank 12
binders_3 <- c("ZNF528", "EBF2", "EBF1")
binders_3 <- sort(binders_3)

#MEME-6 rank 16
binders_4 <- c("RELA")

#MEME-8 rank 17
binders_5 <- c("KLF1", "KLF2", "KLF3", "KLF4", "KLF5", "KLF6", "KLF7", 
       "KLF9", "KLF10", "KLF11", "KLF12", "KLF13", "KLF14", "KLF15", 
       "KLF16", "KLF17", "ZFP281", "ZFP410", "ZFP740", "SP1", "SP2", 
       "SP3", "SP4", "SP5", "SP8", "SP9", "ZNF148", "ZNF263", "ZNF281", 
       "ZNF320", "ZNF530", "ZNF610", "ZNF682", "ZNF740", "SOX13", "SMAD3", 
       "PATZ1", "MAZ", "BCL6B", "E2F6", "ASCL2", "WT1", "EGR1", 
       "PRDM9", "ZBTB14")
binders_5 <- sort(binders_5)

#MEME-17 rank 19
binders_6 <- c("OBOX6", "OTX1", "OTX2", "PITX1", "PITX2")
binders_6 <- sort(binders_6)

# Extract available features (genes) from the Seurat object
available_features <- rownames(strt.seurat.obj)

# Filter the dataframe to find matching rows and extract V1
binders_1_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_1) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

binders_2_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_2) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

binders_3_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_3) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

binders_4_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_4) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

binders_5_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_5) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

binders_6_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_6) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

binders_X_TFEs <- TFE_corresponding_genes %>%
  filter(V2 %in% binders_X) %>%
  filter(V1 %in% available_features) %>%
  mutate(TFE_V2 = paste(V1, V2, sep = "/"))

# NOTE: Before plotting, assign `features` and `labels` to the desired binders set (e.g., binders_1_TFEs, binders_2_TFEs, etc.)
# Only one set should be active at a time for correct plotting.

features <- binders_1_TFEs$V1
labels <- binders_1_TFEs$TFE_V2

features <- binders_2_TFEs$V1
labels <- binders_2_TFEs$TFE_V2

features <- binders_3_TFEs$V1
labels <- binders_3_TFEs$TFE_V2

features <- binders_4_TFEs$V1
labels <- binders_4_TFEs$TFE_V2

features <- binders_5_TFEs$V1
labels <- binders_5_TFEs$TFE_V2

features <- binders_6_TFEs$V1
labels <- binders_6_TFEs$TFE_V2

features <- binders_X_TFEs$V1
labels <- binders_X_TFEs$TFE_V2

plots <- VlnPlot(
  object = strt.seurat.obj, 
  features = features, 
  group.by = "orig.ident", 
  ncol = 3, 
  combine = FALSE,
  layer = "data"
)

# Create individual plots with titles
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]]
  p + 
    scale_fill_manual(values = friendly_pal("wong_eight")) +  # Color scale
    labs(y = "Log2(Normalized Count)", title = labels[i]) +  # Y-axis label and title
    theme_classic() +  # Use classic theme
    theme(
      text = element_text(family = "sans", color = "black"),  # Font settings
      plot.title = element_text(size = 14, face = "bold"),
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.title.y = element_text(size = 18),  # Y-axis title
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      axis.text.y = element_text(size = 16, face = "plain"),  # Y-axis text
      axis.ticks = element_line(size = 0.5),  # Axis ticks
      legend.text = element_text(size = 18)  # Legend text size
    )  +
    guides(fill = guide_legend(override.aes = list(size = 6)))  # Adjust legend key size
})

# set y-axis limits
plots <- lapply(plots, function(p) {
  p + coord_cartesian(ylim = c(-13, -5))  # Set y-axis limits for each plot
})


library(svglite)
svglite(filename = "Binders_1.svg", width = 24, height = 6) #-13, 0, ncol = 7
svglite(filename = "Binders_2.svg", width = 30, height = 12) #-13, -4, ncol = 9
svglite(filename = "Binders_3.svg", width = 21, height = 6) #-13, -4, ncol = 6
svglite(filename = "Binders_4.svg", width = 18, height = 3) #-13, -4, ncol = 5
svglite(filename = "Binders_5.svg", width = 39, height = 27) #-13, 0, ncol =12
svglite(filename = "Binders_6.svg", width = 27, height = 6) #-13, -5, ncol = 8

# Combine plots
wrap_plots(plots) + plot_layout(ncol = 8, guides = 'collect')
dev.off()
