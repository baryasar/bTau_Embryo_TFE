targetPackages <- c("stringr", "dplyr","ggplot2", "cowplot", "ggbeeswarm", "forcats", "Seurat")
for(package in targetPackages) library(package, character.only = T)
library(ggpubfigs)

# Note: Data objects such as 'TFEcounts' need to be loaded or created before running this script.

# Create Seurat object
strt.seurat.obj <- CreateSeuratObject(counts = TFEcounts, project = "STRT", min.cells = 3, min.features = 200)

condition <- c("16c","16c","16c","16c","16c","16c","16c","2c","2c","2c","2c","4c","4c","4c","4c","4c","8c","8c","8c","8c","8c","8c","8c","Blastocyst","Blastocyst","Blastocyst","Blastocyst","Blastocyst","Blastocyst","Blastocyst","GV","GV","GV","GV","MII","MII","MII","MII","MII","MII","MII")
#strt.seurat.obj$orig.ident <- condition

strt.seurat.obj$orig.ident <- factor(condition, 
                                     levels = c('GV', 'MII', '2c', '4c', '8c', '16c', 'Blastocyst'),
                                     labels = c('GV oocyte', 'MII oocyte', '2-cell', '4-cell', '8-cell', '16-cell', 'Blastocyst'))

#Custom Spike-in normalization
reads.p1 <- TFEcounts + 1
all_genes <- rownames(TFEcounts)
ctrlGeneIndices <- grep("^RNA_SPIKE_", all_genes)
reads.p1.spikes <- colSums(reads.p1[ctrlGeneIndices, ])

normalized_data <- as(as.matrix(log2(reads.p1 / rep(reads.p1.spikes, each = nrow(reads.p1)))), "dgCMatrix")

retained_genes <- rownames(strt.seurat.obj)
rownames(normalized_data) <- gsub("_", "-", rownames(normalized_data))
normalized_data_filtered <- normalized_data[intersect(rownames(normalized_data), retained_genes), ]
strt.seurat.obj <- SetAssayData(object = strt.seurat.obj, assay = "RNA", layer = "data", new.data = normalized_data_filtered)

# Check if the feature names match
all.equal(rownames(normalized_data_filtered), rownames(GetAssayData(object = strt.seurat.obj, assay = "RNA", layer = "data")))

Idents(strt.seurat.obj) <- strt.seurat.obj$orig.ident

strt.seurat.obj <- FindVariableFeatures(strt.seurat.obj, selection.method = "vst", nfeatures = 2000)
#The vst method in FindVariableFeatures() uses the data stored in the “data” slot of the Seurat object

all.genes <- row.names(strt.seurat.obj)
strt.seurat.obj <- ScaleData(strt.seurat.obj, features = all.genes)

# PCA plotting
strt.seurat.obj <- RunPCA(strt.seurat.obj, features = VariableFeatures(object = strt.seurat.obj), npcs = 32)
strt.seurat.obj <- JackStraw(strt.seurat.obj, num.replicate = 100) #default number of dimensions to compute the jackstraw for is 20
strt.seurat.obj <- ScoreJackStraw(strt.seurat.obj, dims = 1:20)
JackStrawPlot(strt.seurat.obj, dims = 1:20)
strt.seurat.obj[["pca"]]@jackstraw@overall.p.values
ElbowPlot(strt.seurat.obj)

#PCA PLOT
tiff(filename="PCA_plot.tiff",width=2300, height=2000,res=300)
DimPlot(strt.seurat.obj, 
        reduction = "pca", 
        group.by = "orig.ident", 
        label = FALSE, 
        label.size = 5, 
        repel = TRUE, 
        pt.size = 4, 
        alpha = 0.7) +
  scale_color_manual(values = friendly_pal("wong_eight")) +  # Apply custom color palette
  theme_classic() +
  theme(
    text = element_text(family = "sans", color = "black"),  # Apply built-in sans font to all text elements
    plot.title = element_blank(),                          # Remove plot title
    axis.title.x = element_text(size = 16, face = "bold"), # X-axis title
    axis.title.y = element_text(size = 16, face = "bold"), # Y-axis title
    axis.text = element_text(size = 14, face = "plain"),   # Axis text
    axis.ticks = element_line(size = 0.5),                 # Axis ticks
    legend.text = element_text(size = 16)                 # Legend text
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) + # Adjust legend key size
  labs(x = "PC1", y = "PC2")  # Customize axis labels
dev.off()

# UMAP plotting
strt.seurat.obj <- RunUMAP(strt.seurat.obj, dims = 1:8, verbose = F) #This is UMAP

svg(filename = "UMAP_8dim.svg", width = 8, height = 8) 

DimPlot(strt.seurat.obj, 
        label.size = 5, 
        repel = TRUE, 
        label = FALSE, 
        group.by = "orig.ident", 
        pt.size = 8,
        alpha = 0.8) +
  scale_color_manual(values = friendly_pal("wong_eight")) +
  theme_classic() +
  theme(
    text = element_text(family = "sans", color = "black"),  # Apply built-in sans font to all text elements
    plot.title = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),  # X-axis title
    axis.title.y = element_text(size = 16, face = "bold"),  # Y-axis title
    axis.text = element_text(size = 14, face = "plain"),  # Axis text
    axis.ticks = element_line(size = 0.5),  # Axis ticks
    legend.text = element_text(size = 16),  # Legend text
    #legend.key.size = unit(1, "cm")  # Increase the size of the legend key
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) + # Adjust legend key size
  labs(x = "UMAP 1", y = "UMAP 2") # Customize axis label
dev.off()
