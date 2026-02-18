library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library(dplyr)
library(ggplot2)


countdata<-read.table("Stammera_raw_counts_single.txt", header=T, row.names=1)
head(countdata)
ncol(countdata)
coldata=read.table("Stammera_metadata.txt", header=T,row.names=1)
head(coldata)
nrow(coldata)
coldata$titers <- scale(as.numeric(coldata$titers))  # standardizes titers

# Yellow

# Subset metadata YELLOW 

coldata_yell <- coldata[coldata$Treatment %in% c("YELL_UNT", "YELL_HUM"), ]
countdata_yell <- countdata[, rownames(coldata_yell)]
dds_yell_titers <- DESeqDataSetFromMatrix(
  countData = countdata_yell,
  colData = coldata_yell,
  design =  ~ titers + Replicate+Treatment
)
keep <- rowSums(counts(dds_yell_titers)) > 10
dds_yell_titers <- dds_yell_titers[keep,]
dds_yell_titers <- DESeq(dds_yell_titers)
res_yell_titers <- results(dds_yell_titers,contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))
summary (res_yell_titers)


Figure 7D. Volcano plot comparing gene expression under control vs. low humidity conditions across untreated dsRNA yellow treatments. 

# Volcano plot

p_yell <- EnhancedVolcano(res_yell_titers,
  lab = rownames(res_yell_titers),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  col = c("gray", "gray", "gray", "red3"),
  colAlpha = 0.5,
  pointSize = 5,
  xlim = c(-4, 4),
  ylim = c(0, 4),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  title = NULL,
  subtitle = NULL,
  caption = NULL
)

p_yell + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_blank()
)

ggsave("Fig7D_yellow_volcano_plot.pdf", p_yell)
ggsave("Fig7D_yellow_volcano_plot.svg", p_yell)

# GFP

coldata_GFP <- coldata[coldata$Treatment %in% c("GFP_UNT", "GFP_HUM"), ]
countdata_GFP <- countdata[, rownames(coldata_GFP)]
dds_GFP_titers <- DESeqDataSetFromMatrix(
  countData = countdata_GFP,
  colData = coldata_GFP,
  design =  ~ titers + Replicate+Treatment
)
keep <- rowSums(counts(dds_GFP_titers)) > 10
dds_GFP_titers <- dds_GFP_titers[keep,]
dds_GFP_titers <- DESeq(dds_GFP_titers)
res_GFP_titers <- results(dds_GFP_titers,contrast = c("Treatment", "GFP_HUM", "GFP_UNT"))
summary (res_GFP_titers)


Figure 7D. Volcano plot comparing gene expression under control vs. low humidity conditions across untreated dsRNA GFP treatments. 

# Volcano plot

p_gfp <- EnhancedVolcano(res_GFP_titers,
  lab = rownames(res_GFP_titers),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  col = c("gray", "gray", "gray", "red3"),
  colAlpha = 0.5,
  pointSize = 5,
  xlim = c(-4, 4),
  ylim = c(0, 4),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  title = NULL,
  subtitle = NULL,
  caption = NULL
)

p_gfp + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_blank()
)

ggsave("Fig7D_GFP_volcano_plot.pdf", p_gfp)
ggsave("Fig7D_GFP_volcano_plot.svg", p_gfp)


##############CTL UNT and  CTL HUM ########################

coldata_CTL <- coldata[coldata$Treatment %in% c("CTL_UNT", "CTL_HUM"), ]
countdata_CTL <- countdata[, rownames(coldata_CTL)]
dds_CTL_titers <- DESeqDataSetFromMatrix(
  countData = countdata_CTL,
  colData = coldata_CTL,
  design =  ~ titers + Replicate+Treatment
)
keep <- rowSums(counts(dds_CTL_titers)) > 10
dds_CTL_titers <- dds_CTL_titers[keep,]
dds_CTL_titers <- DESeq(dds_CTL_titers)
res_CTL_titers <- results(dds_CTL_titers,contrast = c("Treatment", "CTL_HUM", "CTL_UNT"))
summary (res_CTL_titers)

Figure 7D. Volcano plot comparing gene expression under control vs. low humidity conditions across untreated control. 

# Volcano plot

p_ctl <- EnhancedVolcano(res_CTL_titers,
  lab = rownames(res_CTL_titers),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  col = c("gray", "gray", "gray", "red3"),
  colAlpha = 0.5,
  pointSize = 5,
  xlim = c(-4, 4),
  ylim = c(0, 4),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  title = NULL,
  subtitle = NULL,
  caption = NULL
)

p_ctl + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.background = element_blank()
)

ggsave("Fig7D_CTL_volcano_plot.pdf", p_gfp)
ggsave("Fig7D_CTL_volcano_plot.svg", p_gfp)


write.csv( as.data.frame(res_yell_titers), file="YELL_DEGs_titers.csv")
write.csv( as.data.frame(res_GFP_titers), file="GFP_DEGs_titers.csv")
write.csv( as.data.frame(res_CTL_titers), file="CTL_DEGs_titers.csv")

# Figure 7E. Heatmap of Stammera gene expression in eggs from dsRNA-yellow females under control and low humidity conditions. 

ntd_yell_titers <- normTransform(dds_yell_titers)
data_ntd_yell_titers<-assay(ntd_yell_titers)

top_genes_yell_titers <- head(order(rowMeans(data_ntd_yell_titers), decreasing = TRUE), 300)
top_data_yell_titers <- data_ntd_yell_titers[top_genes_yell_titers, ]

# Annotation file

annotations <- read.delim("annotation_Stammera.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

annotations$gene <- trimws(annotations$gene)

common_genes <- intersect(rownames(top_data_yell_titers), annotations$gene)

top_data_yell_titers_sub <- top_data_yell_titers[common_genes, ]
annot_sub <- annotations[match(common_genes, annotations$gene), ]

annotation_row <- data.frame(Category = annot_sub$Categories)
rownames(annotation_row) <- annot_sub$gene
ordering <- order(annotation_row$Category)
top_data_yell_titers_sub_ordered <- top_data_yell_titers_sub[ordering, ]
annotation_ordered <- annotation_row[ordering, , drop = FALSE]

annotation_ordered$Category <- factor(annotation_ordered$Category, 
                                      levels = sort(unique(annotation_ordered$Category)))
sorted_categories <- levels(annotation_ordered$Category)
cat_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(sorted_categories))
names(cat_colors) <- sorted_categories
ann_colors <- list(Category = cat_colors)

my_palette <- colorRampPalette(c("#a6dba0","#01665e","#8c510a"))(100)

# heatmap
heatmap_yell<-pheatmap(top_data_yell_titers_sub_ordered,
         scale = "row",
         annotation_row = annotation_ordered,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         fontsize_row = 5,
         color = my_palette)


ggsave("Fig7E_yellow_heatmap.pdf", heatmap_yell)
ggsave("Fig7E_yellow_heatmap.svg", heatmap_yell)

