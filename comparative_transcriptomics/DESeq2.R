
library("DESeq2")
library(vegan)
library(pairwiseAdonis)
library(RVAideMemoire)
library(bbcRNA)
library("pheatmap")
library(ltc)
library("RColorBrewer")
library("EnhancedVolcano")


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


Adult_OV_Adult_FG<-results(dds, contrast=c("Treatment","Adult_OV","Adult_FG"))


#########green#############
# Load required libraries
library(DESeq2)        
library(EnhancedVolcano)

Adult_OV_Adult_FG <- results(dds, contrast = c("Treatment", "Adult_OV", "Adult_FG"))
Adult_OV_Adult_FG <- na.omit(Adult_OV_Adult_FG)

green_color <- rgb(5, 143, 77, maxColorValue = 255)
yellow_color <- rgb(247, 182, 39, maxColorValue = 255)

color_vector <- ifelse(rownames(Adult_OV_Adult_FG) == 'FUN_023368-T1', yellow_color, green_color)

EnhancedVolcano(Adult_OV_Adult_FG,
   lab = rownames(Adult_OV_Adult_FG),
   x = 'log2FoldChange',
   y = 'padj',
   selectLab = c('FUN_023368-T1'), 
   pCutoff = 5e-2,
   FCcutoff = 2.0,
   pointSize = 4,
   labSize = 2,
   labCol = 'black',
   labFace = 'bold',
   boxedLabels = TRUE,
   colAlpha = 4/5,
   legendPosition = 'none',  
   drawConnectors = TRUE,
   widthConnectors = 1.0,
   colConnectors = 'black',
   gridlines.major = FALSE, 
   gridlines.minor = FALSE,  
   col = color_vector,       
   xlab = expression(Log[2]~ 'fold change'), 
   ylab = expression(-Log[10]~ 'p-value'),  
   title = NULL            
)

