library(Gviz)
library(rtracklayer)
library(GenomicRanges)


options(ucscChromosomeNames = FALSE)

stages <- list(
  "GV"     = list(dir = "~/ARS-UCD1.3/bw_files/GV",     file = "BOV_LIB_1.3_12_"),
  "MII"    = list(dir = "~/ARS-UCD1.3/bw_files/MII",    file = "BOV_LIB_1.3_8_"),
  "2-c"    = list(dir = "~/ARS-UCD1.3/bw_files/2c",     file = "BOV_LIB_1.3_24_"),
  "4-c"    = list(dir = "~/ARS-UCD1.3/bw_files/4c",     file = "BOV_LIB_1.3_17_"),
  "8-c"    = list(dir = "~/ARS-UCD1.3/bw_files/8c",     file = "BOV_LIB_1.3_28_"),
  "16-c"   = list(dir = "~/ARS-UCD1.3/bw_files/16c",    file = "BOV_LIB_1.3_33_"),
  "Blc"    = list(dir = "~/ARS-UCD1.3/bw_files/Blc",    file = "BOV_LIB_1.3_42_")
)

colors <- c(
  "GV" = "#E69F00", "MII" = "#56B4E9", "2-c" = "#009E73",
  "4-c" = "#F0E442", "8-c" = "#0072B2", "16-c" = "#D55E00", "Blc" = "#CC79A7"
)

peak_file <- "~/ARS-UCD1.3/BOV_LIB_1.3_peaks.bed"
ensgene_file <- "~/ARS-UCD1.3/ensGene.txt"
gff_dir <- "~/ARS-UCD1.3/HD_annotation"
output_dir <- "~/ARS-UCD1.3"
dir.create(output_dir, showWarnings = FALSE)

gene_info <- list(
  ARGFX  = list(chr="1",  start=66051000, end=66088200, strand="+", gff="HD_ARGFX.gff3", ylim = c(0, 0.02)),
  DUXA   = list(chr="18", start=64390000, end=64407000, strand="-", gff="HD_DUXA.gff3", ylim = c(0, 0.2)),
  LEUTX  = list(chr="18", start=49322000, end=49339500, strand="+", gff="HD_LEUTX.gff3", ylim = c(0, 0.05)),
  NOBOX  = list(chr="4",  start=107659000, end=107674700, strand="-", gff="HD_NOBOX.gff3", ylim = c(0, 0.08)),
  TPRX1  = list(chr="18", start=54678000, end=54692000, strand="-", gff="HD_TPRX1.gff3", ylim = c(0, 0.08)),
  TPRX2  = list(chr="18", start=54724400, end=54728000, strand="+", gff="HD_TPRX2.gff3", ylim = c(0, 0.06)),
  TPRX3  = list(chr="18", start=62770500, end=62777800, strand="-", gff="HD_TPRX3.gff3", ylim = c(0, 0.2)),
  DPPA3  = list(chr="5", start=101325000, end=101341700, strand="+", ylim = c(0, 0.9))
)

plot_gene_tracks <- function(gene_name, info) {
  message("Plotting: ", gene_name)
  chr <- info$chr; start <- info$start; end <- info$end; strand <- info$strand
  gff_path <- file.path(gff_dir, info$gff)
  tracks <- list(GenomeAxisTrack(fontsize = 58,col = "black",          # axis line
                                 col.title = "black",    # title color
                                 col.axis = "black",     # tick labels
                                 fontcolor = "black"
                                 
  ))
  
  for (stage in names(stages)) {
    dir_path <- stages[[stage]]$dir
    file_prefix <- stages[[stage]]$file
    strand_tag <- ifelse(strand == "+", "plus", "minus")
    bw_file <- file.path(dir_path, paste0(file_prefix, strand_tag, ".bw"))
    
    if (file.exists(bw_file)) {
      
      ylim_val <- if (!is.null(info$ylim)) info$ylim else c(0, 0.08)
      
      bw_track <- DataTrack(
        range = bw_file, genome = "bosTau9", chromosome = chr,
        name = stage, type = "h", col = colors[[stage]],
        fill = colors[[stage]], lwd = 10, ylim = ylim_val,
        yTicksAt = ylim_val, cex.axis = 1.2, size = 18,
        col.axis = "black",
        col.title = "black",
        fontcolor = "black",
        fontface.title = 1,
        fontface.group = 1        
      )
      
      tracks <- append(tracks, list(bw_track))
    }
  }
  
  if (file.exists(peak_file)) {
    peak_df <- read.table(peak_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(peak_df)[1:6] <- c("chrom", "start", "end", "name", "score", "strand")
    peak_df <- subset(peak_df, chrom == chr & end >= (start - 10000) & start <= (end + 10000))
    peak_gr <- GRanges(seqnames = peak_df$chrom, ranges = IRanges(start = peak_df$start + 1, end = peak_df$end),
                       strand = peak_df$strand, transcript = peak_df$name, symbol = peak_df$name)
    
    peak_track <- GeneRegionTrack(peak_gr, genome = "bosTau9", chromosome = chr, name = "TFE",
                                  fill = "#77C679", col = NA, shape = "box", stacking = "squish",
                                  transcriptAnnotation = "transcript", fontcolor.group = "black",
                                  col.title = "#77C679", fontsize.group = 77, size = 12,
                                  fontface.title = 1,
                                  fontface.group = 1
                                  )
    
    tracks <- append(tracks, list(peak_track))
  }
  
  gene_table <- read.table(ensgene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")
  colnames(gene_table)[1:16] <- c("bin", "name", "chrom", "strand", "txStart", "txEnd",
                                  "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
                                  "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
  ens_gene_subset <- subset(gene_table, chrom == chr & txEnd >= (start - 10000) & txStart <= (end + 10000))
  
  ens_gr_list <- GRangesList()
  for (i in seq_len(nrow(ens_gene_subset))) {
    row <- ens_gene_subset[i, ]
    starts <- as.numeric(unlist(strsplit(row$exonStarts, ",")))
    ends <- as.numeric(unlist(strsplit(row$exonEnds, ",")))
    strand_i <- row$strand
    cds_start <- row$cdsStart
    cds_end <- row$cdsEnd
    
    for (j in seq_along(starts)) {
      exon_start <- starts[j] + 1
      exon_end <- ends[j]
      if (exon_start >= exon_end) next
      
      transcript_clean <- sub("\\.path1$", "", row$name)
      symbol_clean <- sub("\\.path1$", "", row$name2)
      
      if (exon_end < cds_start || exon_start > cds_end) {
        feature <- ifelse(strand_i == "+", ifelse(exon_end < cds_start, "5UTR", "3UTR"), ifelse(exon_end < cds_start, "3UTR", "5UTR"))
        region <- GRanges(chr, IRanges(exon_start, exon_end), strand = strand_i, feature = feature,
                          transcript = transcript_clean, symbol = symbol_clean)
        ens_gr_list[[paste0(row$name, "_", j)]] <- region
      } else if (exon_start >= cds_start && exon_end <= cds_end) {
        region <- GRanges(chr, IRanges(exon_start, exon_end), strand = strand_i, feature = "CDS",
                          transcript = transcript_clean, symbol = symbol_clean)
        ens_gr_list[[paste0(row$name, "_", j)]] <- region
      } else {
        if (exon_start < cds_start) {
          region <- GRanges(chr, IRanges(exon_start, cds_start - 1), strand = strand_i,
                            feature = ifelse(strand_i == "+", "5UTR", "3UTR"),
                            transcript = transcript_clean, symbol = symbol_clean)
          ens_gr_list[[paste0(row$name, "_", j, "_utrL")]] <- region
        }
        region <- GRanges(chr, IRanges(max(exon_start, cds_start), min(exon_end, cds_end)),
                          strand = strand_i, feature = "CDS",
                          transcript = transcript_clean, symbol = symbol_clean)
        ens_gr_list[[paste0(row$name, "_", j, "_cds")]] <- region
        if (exon_end > cds_end) {
          region <- GRanges(chr, IRanges(cds_end + 1, exon_end), strand = strand_i,
                            feature = ifelse(strand_i == "+", "3UTR", "5UTR"),
                            transcript = transcript_clean, symbol = symbol_clean)
          ens_gr_list[[paste0(row$name, "_", j, "_utrR")]] <- region
        }
      }
    }
  }
  
  if (length(ens_gr_list) > 0) {
    ensTrack <- GeneRegionTrack(unlist(ens_gr_list, use.names = FALSE), genome = "bosTau9", chromosome = chr,
                                transcriptAnnotation = "symbol", name = "Gene",
                                fill = "#F46D6F", col = NA, shape = "box",
                                stacking = "full", col.title = "#F46D6F", fontsize.group = 77, size = 12,
                                thinBoxFeature = c("5UTR", "3UTR"),
                                col.line = "black", fontcolor.group = "black",
                                fontface.title = 1,
                                fontface.group = 1)
    tracks <- append(tracks, list(ensTrack))
  }
  
  if (gene_name != "DPPA3" && file.exists(gff_path)) {
    gff_raw <- import.gff3(gff_path)
    exons <- gff_raw[gff_raw$type == "exon"]
    exons$transcript <- sub("\\.mrna1$", "", exons$Parent)
    
    gffTrack <- GeneRegionTrack(exons, genome = "bosTau9", chromosome = chr,
                                name = "HD", transcriptAnnotation = "transcript",
                                col = NA, fill = "#6BAED6", shape = "box",
                                stacking = "full", fontcolor.group = "black",
                                col.title = "#6BAED6", fontsize.group = 77, size = 12,
                                col.line = "black",
                                fontface.title = 1,
                                fontface.group = 1)
    tracks <- append(tracks, list(gffTrack))
  }
  
  pdf(file.path(output_dir, paste0(gene_name, "_Gviz_plot.pdf")), width = 16, height = 12)
  plotTracks(tracks, from = start, to = end, chromosome = chr,
             background.title = "white",
             rotate.title = TRUE, margin = 24, cex.title = 4,
             col.title = "black",       # track title text color
             col.axis = "black",        # y-axis tick color
             fontface.title = "plain"
  )
  dev.off()
}

for (gene in names(gene_info)) {
  plot_gene_tracks(gene, gene_info[[gene]])
}
