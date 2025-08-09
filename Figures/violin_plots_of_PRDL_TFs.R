library(ggpubfigs)
library(patchwork)
library(svglite)

# NOTE: Some of the objects were pre-loaded in the environment. These objects were loaded or created in other scripts within this project.

TFE_corresponding_genes <- read.delim("TFE_corresponding_genes.txt", header=FALSE)

# Manually assign gene names to each TFE ID
gene_names <- c(
  "TFE1467" = "ARGFX", "TFE1468" = "ARGFX",
  "TFE45904" = "DUXA", "TFE45905" = "DUXA", "TFE45906" = "DUXA",
  "TFE43634" = "LEUTX", "TFE43635" = "LEUTX", "TFE43636" = "LEUTX", "TFE43637" = "LEUTX",
  "TFE97580" = "NOBOX", "TFE97582" = "NOBOX", "TFE97583" = "NOBOX", "TFE97584" = "NOBOX", "TFE97585" = "NOBOX",
  "TFE44646" = "TPRX1", "TFE44647" = "TPRX1",
  "TFE44648" = "TPRX2", "TFE44649" = "TPRX2",
  "TFE45706" = "TPRX3"
)


# NOTE:
# Adjust the `features` vector to select the desired TFEs before plotting.
# Update the `svglite()` filename and `width`/`height` parameters according to the number of features plotted.
# Set the `ncol` argument in `wrap_plots()` to match the number of columns you want in the final layout.


#ARGFX
features <- c("TFE1467","TFE1468")

#DUXA
features <- c("TFE45904","TFE45905","TFE45906")

#LEUTX
features <- c("TFE43634","TFE43635", "TFE43636", "TFE43637")

#NOBOX
features <-  c("TFE97580", "TFE97582", "TFE97583", "TFE97584", "TFE97585")

#TPRX1
features <- c("TFE44646", "TFE44647")

#TPRX2
features <- c("TFE44648", "TFE44649")

#TPRX3
features <-  c("TFE45706")

plots <- VlnPlot(
  object = strt.seurat.obj, 
  features = features, 
  group.by = "orig.ident", 
  ncol = 3, 
  combine = FALSE
)

plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]]
  tfe_id <- features[i]
  gene <- gene_names[[tfe_id]]
  title_expr <- bquote(.(tfe_id) ~ "(" * italic(.(gene)) * ")")
  
  p + 
    scale_fill_manual(values = friendly_pal("wong_eight")) +  
    labs(y = "Log2(Normalized Count)", title = title_expr) +
    theme_classic() +  
    theme(
      text = element_text(family = "sans", color = "black"),
      plot.title = element_text(size = 26, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 24, face = "plain"),
      axis.ticks = element_line(size = 0.5),
      legend.text = element_text(size = 24)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 8.5)))
})

# set y-axis limits
plots <- lapply(plots, function(p) {
  p + coord_cartesian(ylim = c(-13, 0))  # Set y-axis limits for each plot
})

svglite(filename = "ARGFX_TFEs_v3.svg", width = 9, height = 3) #ncol =2 (3.5*2 + 2)
svglite(filename = "DUXA_TFEs_v3.svg", width = 12.5, height = 3) #ncol =3 (3.5*3 + 2)
svglite(filename = "LEUTX_TFEs_v3.svg", width = 16, height = 3) #ncol =4 (3.5*4 + 2)
svglite(filename = "NOBOX_TFEs_v3.svg", width = 19.5, height = 3) #ncol =5 (3.5*5 + 2)
svglite(filename = "TPRX1_TFEs_v3.svg", width = 9, height = 3) #ncol =2 (3.5*2 + 2)
svglite(filename = "TPRX2_TFEs_v3.svg", width = 9, height = 3) #ncol =2 (3.5*2 + 2)
svglite(filename = "TPRX3_TFEs_v3.svg", width = 5.5, height = 3) #ncol =1 (3.5 +2)

wrap_plots(plots) + plot_layout(ncol = 5, guides = 'collect')
dev.off()
