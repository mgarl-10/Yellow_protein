
library("DESeq2")
library(ltc)
library("RColorBrewer")
library("EnhancedVolcano")
library(ggplot2)


#Transcriptome normalization
countdata<-read.table("life_cycle_host_count_table.txt", header=T, row.names=1)
head(countdata)
ncol(countdata)
coldata=read.table("life_cycle_host_metadata.txt", header=T,row.names=1)
head(coldata)
nrow(coldata)
ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = countdata,
colData = coldata,
design = ~ Treatment)
ddsFullCountTable
dds <- DESeq(ddsFullCountTable)
res <- results( dds )
summary (res)
res <- res[order(res$padj),]
head(res)
write.csv( as.data.frame(res), file="Adult_FG_Adult_OV.csv")

Adult_OV_Adult_FG <- results(dds, contrast = c("Treatment","Adult_OV","Adult_FG"))
Adult_OV_Adult_FG <- na.omit(Adult_OV_Adult_FG)


# Volcano plot


Figure 1J. Volcano plot showing differential gene expression between foregut symbiotic organs and ovary-associated glands, analysed via DESeq2.


green_color  <- rgb(5, 143, 77,  maxColorValue = 255)
yellow_color <- rgb(247, 182, 39, maxColorValue = 255)
darker_yellow <- rgb(191, 144, 0,  maxColorValue = 255)


gene <- "FUN_023368-T1"

p <- EnhancedVolcano(
  Adult_OV_Adult_FG,
  lab = rownames(Adult_OV_Adult_FG),
  x = "log2FoldChange",
  y = "padj",
  selectLab = gene,
  pCutoff = 5e-2,
  FCcutoff = 2.0,
  pointSize = 5,
  labSize = 2,
  labCol = "black",
  labFace = "bold",
  boxedLabels = TRUE,
  legendPosition = "none",
  drawConnectors = TRUE,
  widthConnectors = 1,
  colConnectors = "black",
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  cutoffLineType = 'dashed',
  cutoffLineWidth = 0.6,
  col = c("grey70", green_color, "grey70", green_color),
  title = NULL,
  subtitle = NULL,
  caption = NULL
)

df <- as.data.frame(Adult_OV_Adult_FG)
df$gene <- rownames(df)

p +
  geom_point(
    data = df[df$gene == gene, ],
    aes(x = log2FoldChange, y = -log10(padj)),
    inherit.aes = FALSE,
    shape = 21,                
    fill = yellow_color,       
    color = darker_yellow,     
    size = 6,
    stroke = 1.2               
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(color = "black")
  )


ggsave("Fig1J_volcano_plot.pdf", p)
ggsave("Fig1J_volcano_plot.svg", p)
