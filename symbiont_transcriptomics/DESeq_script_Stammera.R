library("DESeq2")
library(vegan)
library(pairwiseAdonis)
library(RVAideMemoire)
library(bbcRNA)
library("pheatmap")
library(ltc)
library("RColorBrewer")
library("EnhancedVolcano")


countdata<-read.table("Stammera_raw_counts_single.txt", header=T, row.names=1)
head(countdata)
ncol(countdata)
coldata=read.table("Stammera_metadata.txt", header=T,row.names=1)
head(coldata)
nrow(coldata)
coldata$titers <- scale(as.numeric(coldata$titers))  # standardizes titers
ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = countdata,
colData = coldata,
design = ~ Replicate+Treatment)
ddsFullCountTable
dds <- DESeq(ddsFullCountTable,fitType="mean")
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
res <- results(dds)
summary (res)
res <- res[order(res$padj),]
head(res)
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup="Treatment")
plotPCA(rld, intgroup=c("Treatment", "Replicate"))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE,fitType="mean")
plotPCA(vsd, intgroup="Treatment")
plotPCA(vsd, intgroup=c("Treatment", "Replicate"))

YELL<-results(dds, contrast=c("Treatment","YELL_UNT","YELL_HUM"))
summary(YELL)
YELL <- YELL[order(res$padj),]
as.data.frame(YELL)

GFP<-results(dds, contrast=c("Treatment","GFP_UNT","GFP_HUM"))
summary(GFP)
GFP <- GFP[order(res$padj),]
head(GFP)
CTL<-results(dds, contrast=c("Treatment","CTL_UNT","CTL_HUM"))
summary(CTL)
CTL <- CTL[order(res$padj),]
head(CTL)



ntd <- normTransform(dds)
data_ntd<-assay(ntd)
pheatmap(data_ntd, cluster_cols=FALSE)
pheatmap(data_ntd,cluster_rows=FALSE,cluster_cols=FALSE,show_rownames=FALSE,color = hcl.colors(100,"YlOrBr",rev=TRUE))

top_genes <- head(order(rowMeans(data_ntd), decreasing = TRUE), 60)
top_data <- data_ntd[top_genes, ]
annotation_col <- as.data.frame(colData(dds)[, c("Treatment", "Replicate")])
pheatmap(top_data,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardize each row (gene)
         fontsize_row = 6)

pheatmap(top_data,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         fontsize_row = 6)


barplot(colSums(counts(dds)), las=2, main="Total Reads Per Sample")

colSums(counts(dds))

dds <- dds[rowSums(counts(dds)) >= 10, ]

####################Analysis including paired counts and caplets###################
countdata<-read.table("Stammera_raw_counts_caplets.txt", header=T, row.names=1)
head(countdata)
ncol(countdata)
coldata=read.table("Stammera_metadata_caplets.txt", header=T,row.names=1)
head(coldata)
nrow(coldata)
ddsFullCountTable <- DESeqDataSetFromMatrix(
countData = countdata,
colData = coldata,
design = ~ Replicate+Treatment)
ddsFullCountTable
dds <- DESeq(ddsFullCountTable)
res <- results( dds )
summary (res)
res <- res[order(res$padj),]
head(res)
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup="Treatment")
plotPCA(rld, intgroup=c("Treatment", "Replicate"))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE,fitType="mean")
plotPCA(vsd, intgroup="Treatment")
plotPCA(vsd, intgroup=c("Treatment", "Replicate"))
pcaData <- plotPCA(vsd, intgroup=c("Treatment", "Replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


###################edgeR####################################

library(edgeR)
coldata_yell <- coldata[coldata$Treatment %in% c("YELL_UNT", "YELL_HUM"), ]
countdata_yell <- countdata[, rownames(coldata_yell)]

group <- coldata_yell$Treatment  
dge <- DGEList(counts = countdata_yell, group = group)

######coldata$titers <- scale(as.numeric(coldata$titers))  # same as before

design <- model.matrix(~ Replicate+Treatment, data = coldata_yell)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = "TreatmentYELL_HUM")  # change as needed
topTags(qlf)
results <- topTags(qlf, n = Inf)
summary(decideTests(qlf))


# Load edgeR
library(edgeR)


group <- coldata$Treatment  # Make sure it's a factor
dge <- DGEList(counts = countdata, group = group)
keep <- filterByExpr(dge, group = group)
dge <- dge[keep,, keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + group)  # No intercept; lets us compare any group directly
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
contrast <- makeContrasts(YELL_HUM - YELL_UNT, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)
topTags(qlf)  # shows top differentially expressed genes

# Optional: save full results to CSV
res <- topTags(qlf, n = Inf)
write.csv(as.data.frame(res), file = "YELL_HUM_vs_YELL_UNT_edgeR.csv")


#######DESEq with bacterial fraction#####################totalreads/Stammera reads#################

# Load required packages
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

countdata <- read.table("Stammera_raw_counts_single.txt", header=TRUE, row.names=1)
coldata <- read.table("Stammera_metadata.txt", header=TRUE, row.names=1)

total_reads <- c(15016769,23441654,19692030,23844717,20029013,20707343,
                 16076515,15966243,20828941,16173123,19655053,19717939,
                 18832436,19962178,18531906,18047044,14441033,15402674)

bacterial_reads <- c(125397,358937,43595,553418,1057516,276091,
                     91236,45867,612619,259789,46053,723446,
                     568544,105180,303448,60008,75145,3606)

samples <- colnames(countdata)
coldata$Sample <- samples
coldata$Total_Reads <- total_reads
coldata$Bacterial_Reads <- bacterial_reads
coldata$Bacterial_Fraction <- bacterial_reads / total_reads
coldata$Bacterial_Fraction_Scaled <- scale(coldata$Bacterial_Fraction)

coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Replicate <- as.factor(coldata$Replicate)

# Model 1: Treatment and Replicate only
dds_noBact <- DESeqDataSetFromMatrix(countData = countdata,
                                     colData = coldata,
                                     design = ~ Replicate+Treatment)
dds_noBact <- DESeq(dds_noBact)
keep <- rowSums(counts(dds_noBact)) > 10
dds_noBact <- dds_noBact[keep,]

# Model 2: Bacterial Fraction + Replicate + Treatment
dds_withBact <- DESeqDataSetFromMatrix(countData = countdata,
                                       colData = coldata,
                                       design = ~ Bacterial_Fraction_Scaled + Replicate+Treatment)
dds_withBact <- DESeq(dds_withBact)
keep <- rowSums(counts(dds_withBact)) > 10
dds_withBact <- dds_withBact[keep,]

# Model 3: Bacterial Titers + Replicate + Treatment
dds_titers <- DESeqDataSetFromMatrix(countData = countdata,
                                       colData = coldata,
                                       design = ~ titers + Replicate+Treatment)
dds_titers <- DESeq(dds_titers)
keep <- rowSums(counts(dds_titers)) > 10
dds_titers <- dds_titers[keep,]

# rlog transforms and PCA
rld_noBact <- rlog(dds_noBact, blind = TRUE)
rld_withBact <- rlog(dds_withBact, blind = TRUE)
plotPCA(rld_noBact, intgroup = "Treatment")
plotPCA(rld_withBact, intgroup = "Treatment")

# PCA colored by Bacterial RNA
pcaData <- plotPCA(rld_withBact, intgroup = "Treatment", returnData = TRUE)
pcaData$BacterialFraction <- coldata$Bacterial_Fraction

ggplot(pcaData, aes(x = PC1, y = PC2, color = BacterialFraction, shape = Treatment)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "PCA: Colored by Bacterial Fraction")


######Overall heatmap with no bacterial fraction###############
ntd_noBact <- normTransform(dds_noBact)
data_ntd_noBact<-assay(ntd_noBact)

top_genes_noBact <- head(order(rowMeans(data_ntd_noBact), decreasing = TRUE), 60)
top_data_noBact <- data_ntd_noBact[top_genes_noBact, ]
annotation_col_noBact <- as.data.frame(colData(dds_noBact)[, c("Treatment", "Replicate")])
pheatmap(top_data_noBact,
         annotation_col_noBact = annotation_col_noBact,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardise each row (gene)
         fontsize_row = 6)

pheatmap(top_data_yell_bact,
         annotation_col_yell_bact = annotation_col_yell_bact,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6)

######Overall heatmap with bacterial fraction###############
ntd_withBact <- normTransform(dds_withBact)
data_ntd_withBact<-assay(ntd_withBact)

top_genes_withBact <- head(order(rowMeans(data_ntd_withBact), decreasing = TRUE), 88)
top_data_withBact <- data_ntd_withBact[top_genes_withBact, ]
annotation_col_withBact <- as.data.frame(colData(dds_withBact)[, c("Treatment", "Replicate")])
pheatmap(top_data_withBact,
         annotation_col_withBact = annotation_col_withBact,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardise each row (gene)
         fontsize_row = 6)

######Overall heatmap with titers###############
ntd_titers <- normTransform(dds_titers)
data_ntd_titers<-assay(ntd_titers)

top_genes_titers <- head(order(rowMeans(data_ntd_titers), decreasing = TRUE), 60)
top_data_titers <- data_ntd_withBact[top_genes_titers, ]
annotation_col_titers <- as.data.frame(colData(dds_titers)[, c("Treatment", "Replicate")])
pheatmap(top_data_titers,
         annotation_col_titers = annotation_col_titers,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardise each row (gene)
         fontsize_row = 6)


****ALL HEATMAPS LOOOK THE SAME**********

#DE gene comparison for YELL_HUM vs YELL_UNT

# DE results without bacterial fraction
res_yell_noBact <- results(dds_noBact, contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))

# DE results with bacterial fraction
res_yell_withBact <- results(dds_withBact, contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))

# Compare number of significant genes
cat("Significant DE genes (FDR < 0.05):\n")
cat("Without bacterial fraction:", sum(res_noBact$padj < 0.05, na.rm = TRUE), "\n")
cat("With bacterial fraction   :", sum(res_withBact$padj < 0.05, na.rm = TRUE), "\n")
Without bacterial fraction: 11 
With bacterial fraction   : 5 

# Scatter plot of log2 fold changes
plot(res_noBact$log2FoldChange, res_withBact$log2FoldChange,
     xlab = "Without Bacterial Fraction",
     ylab = "With Bacterial Fraction",
     main = "Log2FC Comparison")
abline(0, 1, col = "red")

# RESULTS WITH BACTERIAL FRACTION#
res_yell_withBact <- results(dds_withBact, contrast = c("Treatment", "YELL_HUM", "YELL_UNT"), name = "Bacterial_Fraction_Scaled")

EnhancedVolcano(res_yell_withBact,
                lab = rownames(res_yell_withBact),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial RNA Fraction",
                subtitle = "Yellow low humidity vs Yellow control condition",
                pCutoff = 0.05,
                FCcutoff = 1)
res_GFP_withBact <- results(dds_withBact, contrast = c("Treatment", "GFP_HUM", "GFP_UNT"), name = "Bacterial_Fraction_Scaled")

EnhancedVolcano(res_GFP_withBact,
                lab = rownames(res_GFP_withBact),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial RNA Fraction",
                subtitle = "GFP low humidity vs GFP control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

res_CTL_withBact <- results(dds_withBact, contrast = c("Treatment", "CTL_HUM", "CTL_UNT"), name = "Bacterial_Fraction_Scaled")

EnhancedVolcano(res_CTL_withBact,
                lab = rownames(res_CTL_withBact),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial RNA Fraction",
                subtitle = "CTL low humidity vs CTL control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

#RESULTS NO bacterial fraction effect#
res_yell_noBact <- results(dds_noBact, contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))

EnhancedVolcano(res_yell_noBact,
                lab = rownames(res_yell_noBact),
                x = "log2FoldChange",
                y = "padj",
                title = "No Bacterial RNA Fraction",
                subtitle = "Yellow low humidity vs Yellow control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

res_GFP_noBact <- results(dds_noBact, contrast = c("Treatment", "GFP_HUM", "GFP_UNT"))

EnhancedVolcano(res_GFP_noBact,
                lab = rownames(res_GFP_noBact),
                x = "log2FoldChange",
                y = "padj",
                title = "No Bacterial RNA Fraction",
                subtitle = "GFP low humidity vs GFP control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

res_CTL_noBact <- results(dds_noBact, contrast = c("Treatment", "CTL_HUM", "CTL_UNT"))

EnhancedVolcano(res_CTL_noBact,
                lab = rownames(res_CTL_noBact),
                x = "log2FoldChange",
                y = "padj",
                title = "No Bacterial RNA Fraction",
                subtitle = "CTL low humidity vs CTL control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

#RESULTS BACTERIAL TITER EFFECTS#
res_yell_titers <- results(dds_titers, contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))

EnhancedVolcano(res_yell_titers,
                lab = rownames(res_yell_titers),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial Titers",
                subtitle = "Yellow low humidity vs Yellow control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

res_GFP_titers <- results(dds_titers, contrast = c("Treatment", "GFP_HUM", "GFP_UNT"))

EnhancedVolcano(res_GFP_titers,
                lab = rownames(res_GFP_titers),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial Titers",
                subtitle = "GFP low humidity vs GFP control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

res_CTL_titers <- results(dds_titers, contrast = c("Treatment", "CTL_HUM", "CTL_UNT"))

EnhancedVolcano(res_CTL_titers,
                lab = rownames(res_CTL_titers),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial Titers",
                subtitle = "CTL low humidity vs CTL control condition",
                pCutoff = 0.05,
                FCcutoff = 1)

# 
write.csv(as.data.frame(res_noBact), "DESeq2_YELL_noBacterialFraction.csv")
write.csv(as.data.frame(res_withBact), "DESeq2_YELL_withBacterialFraction.csv")
write.csv(as.data.frame(res_bact), "DESeq2_BacterialFraction_Effect.csv")


##############extracting YELL UNT and YELL HUM 
# Subset metadata YELLOW with Bacterial fractions
coldata_yell <- coldata[coldata$Treatment %in% c("YELL_UNT", "YELL_HUM"), ]
countdata_yell <- countdata[, rownames(coldata_yell)]
dds_yell_bact <- DESeqDataSetFromMatrix(
  countData = countdata_yell,
  colData = coldata_yell,
  design =  ~ Bacterial_Fraction_Scaled + Replicate+Treatment
)
keep <- rowSums(counts(dds_yell_bact)) > 10
dds_yell_bact <- dds_yell_bact[keep,]
dds_yell_bact <- DESeq(dds_yell_bact)
res_yell_bact <- results(dds_yell_bact,contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))
summary (res_yell_bact)
EnhancedVolcano(res_yell_bact,
                lab = rownames(res_yell_bact),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial RNA Fraction",
                subtitle = "Yellow low humidity vs Yellow control",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.7,
                pointSize = 3.0)

write.csv( as.data.frame(res_yell_bact), file="yell_DEGs_bacterial_fraction.csv")
normalized_counts_yell<-counts(dds_yell, normalized=TRUE)
write.csv(as.data.frame(normalized_counts_yell),file="normalized_counts_DESeq2_yell.csv")


# Subset metadata YELLOW with no  Bacterial fractions
coldata_yell <- coldata[coldata$Treatment %in% c("YELL_UNT", "YELL_HUM"), ]
countdata_yell <- countdata[, rownames(coldata_yell)]
dds_yell<- DESeqDataSetFromMatrix(
  countData = countdata_yell,
  colData = coldata_yell,
  design =  ~ Replicate+Treatment)
keep <- rowSums(counts(dds_yell)) > 10
dds_yell <- dds_yell[keep,]
dds_yell <- DESeq(dds_yell)
res_yell <- results(dds_yell,contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))
summary (res_yell)
EnhancedVolcano(res_yell,
                lab = rownames(res_yell),
                x = "log2FoldChange",
                y = "padj",
                title = "No Bacterial RNA Fraction",
                subtitle = "Yellow low humidity vs Yellow control",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.7,
                pointSize = 3.0)

##############extracting GFP UNT and GFP HUM 
# Subset metadata GFP with Bacterial fractions
coldata_GFP <- coldata[coldata$Treatment %in% c("GFP_UNT", "GFP_HUM"), ]
countdata_GFP <- countdata[, rownames(coldata_GFP)]
dds_GFP_bact <- DESeqDataSetFromMatrix(
  countData = countdata_GFP,
  colData = coldata_GFP,
  design =  ~ Bacterial_Fraction_Scaled + Replicate+Treatment
)
keep <- rowSums(counts(dds_GFP_bact)) > 10
dds_GFP_bact <- dds_GFP_bact[keep,]
dds_GFP_bact <- DESeq(dds_GFP_bact)
res_GFP_bact <- results(dds_GFP_bact,contrast = c("Treatment", "GFP_HUM", "GFP_UNT"))
summary (res_GFP_bact)
EnhancedVolcano(res_GFP_bact,
                lab = rownames(res_GFP_bact),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial RNA Fraction",
                subtitle = "GFP low humidity vs GFP control",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.7,
                pointSize = 3.0)

write.csv( as.data.frame(res_GFP_bact), file="GFP_DEGs_bacterial_fraction.csv")

# Subset metadata GFP with no  Bacterial fractions
coldata_GFP <- coldata[coldata$Treatment %in% c("GFP_UNT", "GFP_HUM"), ]
countdata_GFP <- countdata[, rownames(coldata_GFP)]
dds_GFP<- DESeqDataSetFromMatrix(
  countData = countdata_GFP,
  colData = coldata_GFP,
  design =  ~ Replicate+Treatment)
keep <- rowSums(counts(dds_GFP)) > 10
dds_GFP <- dds_GFP[keep,]
dds_GFP <- DESeq(dds_GFP)
res_GFP <- results(dds_GFP,contrast = c("Treatment", "GFP_HUM", "GFP_UNT"))
summary (res_GFP)
EnhancedVolcano(res_GFP,
                lab = rownames(res_GFP),
                x = "log2FoldChange",
                y = "padj",
                title = "No Bacterial RNA Fraction",
                subtitle = "GFP low humidity vs GFP control",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.7,
                pointSize = 3.0)

##############extracting CTL UNT and  CTL HUM 
# Subset metadata CTL with Bacterial fractions
coldata_CTL <- coldata[coldata$Treatment %in% c("CTL_UNT", "CTL_HUM"), ]
countdata_CTL <- countdata[, rownames(coldata_CTL)]
dds_CTL_bact <- DESeqDataSetFromMatrix(
  countData = countdata_CTL,
  colData = coldata_CTL,
  design =  ~ Bacterial_Fraction_Scaled + Replicate+Treatment
)
keep <- rowSums(counts(dds_CTL_bact)) > 10
dds_CTL_bact <- dds_CTL_bact[keep,]
dds_CTL_bact <- DESeq(dds_CTL_bact)
res_CTL_bact <- results(dds_CTL_bact,contrast = c("Treatment", "CTL_HUM", "CTL_UNT"))
summary (res_CTL_bact)
EnhancedVolcano(res_CTL_bact,
                lab = rownames(res_CTL_bact),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial RNA Fraction",
                subtitle = "Control low humidity vs Control CTL conditions",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.7,
                pointSize = 3.0)


write.csv( as.data.frame(res_CTL_bact), file="CTL_DEGs_bacterial_fraction.csv")

# Subset metadata CTL with no  Bacterial fractions
coldata_CTL <- coldata[coldata$Treatment %in% c("CTL_UNT", "CTL_HUM"), ]
countdata_CTL <- countdata[, rownames(coldata_CTL)]
dds_CTL<- DESeqDataSetFromMatrix(
  countData = countdata_CTL,
  colData = coldata_CTL,
  design =  ~ Replicate+Treatment)
keep <- rowSums(counts(dds_CTL)) > 10
dds_CTL <- dds_CTL[keep,]
dds_CTL <- DESeq(dds_CTL)
res_CTL <- results(dds_CTL,contrast = c("Treatment", "CTL_HUM", "CTL_UNT"))
summary (res_CTL)
EnhancedVolcano(res_CTL,
                lab = rownames(res_CTL),
                x = "log2FoldChange",
                y = "padj",
                title = "No Bacterial RNA Fraction",
                subtitle = "Control low humidity vs Control CTL conditions",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.7,
                pointSize = 3.0)

##############extracting YELL UNT and YELL HUM 
# Subset metadata YELLOW with titers
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

vsd_yell <- varianceStabilizingTransformation(dds_yell_titers, blind=TRUE,fitType="mean")

vsd_CTL <- varianceStabilizingTransformation(dds_CTL_titers, blind=TRUE,fitType="mean")
vsd_GFP <- varianceStabilizingTransformation(dds_GFP_titers, blind=TRUE,fitType="mean")
vsd_titers <- varianceStabilizingTransformation(dds_titers, blind=TRUE,fitType="mean")


pca_yell <- plotPCA(vsd_yell, intgroup = "Treatment", returnData = TRUE)
ggplot(pca_yell, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 4) +
  ggforce::geom_mark_hull(aes(label = Treatment, fill = Treatment), alpha = 0.1, show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "PCA - YELL samples",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )

pca_titers <- plotPCA(vsd_titers, intgroup = "Treatment", returnData = TRUE)
ggplot(pca_titers, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 4) +
  ggforce::geom_mark_hull(aes(label = Treatment, fill = Treatment), alpha = 0.1, show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "PCA - YELL samples",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )

pca_GFP <- plotPCA(vsd_GFP, intgroup = "Treatment", returnData = TRUE)
ggplot(pca_GFP, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 4) +
  ggforce::geom_mark_hull(aes(label = Treatment, fill = Treatment), alpha = 0.1, show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "PCA - YELL samples",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )

pca_CTL <- plotPCA(vsd_CTL, intgroup = "Treatment", returnData = TRUE)
ggplot(pca_CTL, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 4) +
  ggforce::geom_mark_hull(aes(label = Treatment, fill = Treatment), alpha = 0.1, show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "PCA - YELL samples",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )
library(vegan)		
expr_mat <- t(assay(vsd_titers))
meta <- as.data.frame(colData(dds_titers))
dist_matrix <- dist(expr_mat, method = "euclidean")

# Test for homogeneity of group dispersions
bd <- betadisper(dist_matrix, meta$Treatment)
anova(bd)

EnhancedVolcano(res_yell_titers,
                lab = rownames(res_yell_titers),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial Titers",
                subtitle = "Yellow low humidity vs Yellow control",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.5,
                pointSize = 3.0,
                xlim = c(-4, 4),
                ylim = c(0, 4)
)

write.csv( as.data.frame(res_yell_titers), file="YELL_DEGs_titers.csv")
write.csv( as.data.frame(res_GFP_titers), file="GFP_DEGs_titers.csv")
write.csv( as.data.frame(res_CTL_titers), file="CTL_DEGs_titers.csv")

##############extracting GFP UNT and GFP HUM 
# Subset metadata GFP with Titers
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
EnhancedVolcano(res_GFP_titers,
                lab = rownames(res_GFP_titers),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial Titers",
                subtitle = "GFP low humidity vs GFP control",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.5,
                pointSize = 3.0,
                xlim = c(-4, 4),
                ylim = c(0, 4)
)

##############extracting CTL UNT and  CTL HUM 
# Subset metadata CTL with Bacterial Titers
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
EnhancedVolcano(res_CTL_titers,
                lab = rownames(res_CTL_titers),
                x = "log2FoldChange",
                y = "padj",
                title = "Effect of Bacterial Titers",
                subtitle = "Control low humidity vs Control CTL conditions",
                pCutoff = 0.05,
                FCcutoff = 1,
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 0.5,
                pointSize = 3.0,
                xlim = c(-4, 4),
                ylim = c(0, 4)
)


##########heatmap YELL no bacterial fractions##################
ntd_yell <- normTransform(dds_yell)
data_ntd_yell<-assay(ntd_yell)

top_genes_yell <- head(order(rowMeans(data_ntd_yell), decreasing = TRUE), 230)
top_data_yell <- data_ntd_yell[top_genes_yell, ]
annotation_col_yell <- as.data.frame(colData(dds_yell)[, c("Treatment", "Replicate")])
pheatmap(top_data_yell,
         annotation_col_yell = annotation_col_yell,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardize each row (gene)
         fontsize_row = 6)

pheatmap(top_data_yell,
         annotation_col_yell = annotation_col_yell,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         fontsize_row = 6)

##########heatmap YELL bacterial fractions##################
ntd_yell_bact <- normTransform(dds_yell_bact)
data_ntd_yell_bact<-assay(ntd_yell_bact)

top_genes_yell_bact <- head(order(rowMeans(data_ntd_yell_bact), decreasing = TRUE), 88)
top_data_yell_bact <- data_ntd_yell_bact[top_genes_yell_bact, ]
annotation_col_yell_bact <- as.data.frame(colData(dds_yell_bact)[, c("Treatment", "Replicate")])
pheatmap(top_data_yell_bact,
         annotation_col_yell_bact = annotation_col_yell_bact,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",  # Standardise each row (gene)
         fontsize_row = 6)

pheatmap(top_data_yell_bact,
         annotation_col_yell_bact = annotation_col_yell_bact,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6)

##########heatmap YELL bacterial titers##################
ntd_yell_titers <- normTransform(dds_yell_titers)
data_ntd_yell_titers<-assay(ntd_yell_titers)

top_genes_yell_titers <- head(order(rowMeans(data_ntd_yell_titers), decreasing = TRUE), 230)
top_data_yell_titers <- data_ntd_yell_titers[top_genes_yell_titers, ]
annotation_col_yell_titers <- as.data.frame(colData(dds_yell_titers)[, c("Treatment", "Replicate")])
pheatmap(top_data_yell_titers,
         annotation_col_yell_titers = annotation_col_yell_titers,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardise each row (gene)
         fontsize_row = 6)


pheatmap(top_data_titers,
         annotation_col_titers = annotation_col_titers,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",  # Standardise each row (gene)
         fontsize_row = 6)
pheatmap(top_data_yell_titers,
         annotation_col_yell_titers = annotation_col_yell_titers,
         show_rownames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_row = 6)

#########BAR PLOTS FOR Genes with |log2FC|#################

library(ggplot2)

deg_counts <- data.frame(
  Condition = rep(c("CTL", "GFP", "YELL"), each = 2),
  Direction = rep(c("Downregulated", "Upregulated"), times = 3),
  Count = c(17, 36, 17, 20, 47, 31)
)

ggplot(deg_counts, aes(x = Condition, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Genes with |log2FC| > 1 after controlling for bacterial titers",
       y = "Number of genes",
       fill = "Direction")

###########STATS FOR CONTOL VS YELLOW################
# Downregulated gene counts (log2FC < -1)
yell_down <- 47
ctl_down <- 17

# Not downregulated
yell_not_down <- 78 - 47
ctl_not_down <- 53 - 17

# Create the contingency table
table_yell_ctl <- matrix(c(yell_down, yell_not_down,
                           ctl_down, ctl_not_down),
                         nrow = 2,
                         byrow = TRUE)

rownames(table_yell_ctl) <- c("YELL_HUM", "CTL_HUM")
colnames(table_yell_ctl) <- c("Downregulated", "Not_down")

# Run Fisher’s Exact Test
fisher.test(table_yell_ctl)

	Fisher's Exact Test for Count Data

data:  table_yell_ctl
p-value = 0.002338
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.453539 7.174989
sample estimates:
odds ratio 
  3.181033 
chisq.test(table_yell_ctl)

	Pearson's Chi-squared test with Yates' continuity correction

data:  table_yell_ctl
X-squared = 8.9338, df = 1, p-value = 0.002799


###########STATS FOR GFP VS YELLOW################
# Downregulated gene counts (log2FC < -1)
yell_down <- 47
gfp_down <- 17

# Not downregulated
yell_not_down <- 78 - 47
gfp_not_down <- 37 - 17

# Create the contingency table
table_yell_gfp <- matrix(c(yell_down, yell_not_down,
                           gfp_down, gfp_not_down),
                         nrow = 2,
                         byrow = TRUE)

rownames(table_yell_gfp) <- c("YELL_HUM", "GFP_HUM")
colnames(table_yell_gfp) <- c("Downregulated", "Not_down")

# Run Fisher’s Exact Test
fisher.test(table_yell_gfp)


 fisher.test(table_yell_gfp)

	Fisher's Exact Test for Count Data

data:  table_yell_gfp
p-value = 0.1646
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.751835 4.240915
sample estimates:
odds ratio 
  1.774581 
chisq.test(table_yell_gfp)

	Pearson's Chi-squared test with Yates' continuity correction

data:  table_yell_gfp
X-squared = 1.5429, df = 1, p-value = 0.2142


###########STATS FOR CONTOL VS YELLOW################
# Upregulated gene counts (log2FC < -1)
yell_up <- 31
ctl_up <- 36

# Not downregulated
yell_not_up <- 78 - 31
ctl_not_up <- 53 - 36

# Create the contingency table
table_yell_ctl_up <- matrix(c(yell_up, yell_not_up,
                           ctl_up, ctl_not_up),
                         nrow = 2,
                         byrow = TRUE)

rownames(table_yell_ctl_up) <- c("YELL_HUM", "CTL_HUM")
colnames(table_yell_ctl_up) <- c("Upregulated", "Not_up")

chisq.test(table_yell_ctl_up)

# Run Fisher’s Exact Test
fisher.test(table_yell_ctl_up)

	Fisher's Exact Test for Count Data

data:  table_yell_ctl_up
p-value = 0.002338
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.1393730 0.6879763
sample estimates:
odds ratio 
 0.3143633 
chisq.test(table_yell_ctl_up)

	Pearson's Chi-squared test with Yates' continuity correction

data:  table_yell_ctl_up
X-squared = 8.9338, df = 1, p-value = 0.002799


###########STATS FOR GFP VS YELLOW################
# Upregulated gene counts (log2FC < -1)
ctl_up <- 36
gfp_up <- 20

# Not unpregulated
ctl_not_up <- 53 - 36
gfp_not_up <- 37 - 20

# Create the contingency table
table_ctl_gfp_up <- matrix(c(ctl_up, ctl_not_up,
                           gfp_down, gfp_not_down),
                         nrow = 2,
                         byrow = TRUE)

rownames(table_ctl_gfp_up) <- c("YELL_HUM", "GFP_HUM")
colnames(table_ctl_gfp_up) <- c("Upregulated", "Not_down")

# Run Fisher’s Exact Test
fisher.test(table_ctl_gfp_up)

	Fisher's Exact Test for Count Data

data:  table_ctl_gfp_up
p-value = 0.05034
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.9608447 6.4839238
sample estimates:
odds ratio 
  2.465081 

chisq.test(table_ctl_gfp_up)

	Pearson's Chi-squared test with Yates' continuity correction

data:  table_ctl_gfp_up
X-squared = 3.4871, df = 1, p-value = 0.06185

##########BOTH UP AND DOWN
table_yell_ctl <- matrix(c(47, 31, 17, 36), nrow = 2,
                         dimnames = list(Direction = c("Down", "Up"),
                                         Group = c("YELL", "CTL")))

fisher.test(table_yell_ctl)

	Fisher's Exact Test for Count Data

data:  table_yell_ctl
p-value = 0.002338
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.453539 7.174989
sample estimates:
odds ratio 
  3.181033 

table_yell_gfp <- matrix(c(47, 31, 17, 20), nrow = 2,
                         dimnames = list(Direction = c("Down", "Up"),
                                         Group = c("YELL", "GFP")))

fisher.test(table_yell_gfp)

######COMBINED###############################
table_combined <- matrix(c(47, 31, 34, 56), nrow = 2,
                         dimnames = list(Direction = c("Down", "Up"),
                                         Group = c("YELL", "REF")))
fisher.test(table_combined)

	Fisher's Exact Test for Count Data

data:  table_combined
p-value = 0.005199
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.280892 4.881487
sample estimates:
odds ratio 
  2.483102 


################Box plot for rpoC###################

###################YELL##############

library(DESeq2)
library(ggplot2)

norm_counts_yell <- counts(dds_yell_titers, normalized = TRUE)
head(norm_counts_yell)
gene <- "rpoC"

if (!gene %in% rownames(norm_counts)) {
  stop("Gene not found — check ID format (case-sensitive?)")
}

rpoC_counts_yell <- data.frame(
  Expression = norm_counts_yell[gene, ],
  Sample = colnames(norm_counts_yell),
  Treatment = colData(dds_yell_titers)$Treatment
)

print(rpoC_counts_yell)


ggplot(rpoC_counts, aes(x = Treatment, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  theme_minimal() +
  labs(title = paste("Normalized expression of", gene),
       x = "Treatment",
       y = "DESeq2 normalized counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Wilcoxon rank-sum test (non-parametric)
wilcox_test <- wilcox.test(Expression ~ Treatment, data = rpoC_counts)

# Output result
print(wilcox_test)



###################GFP###################

norm_counts_GFP <- counts(dds_GFP_titers, normalized = TRUE)
head(norm_counts_GFP)
gene <- "rpoC"

if (!gene %in% rownames(norm_counts_GFP)) {
  stop("Gene not found — check ID format (case-sensitive?)")
}

rpoC_counts_GFP <- data.frame(
  Expression = norm_counts_GFP[gene, ],
  Sample = colnames(norm_counts_GFP),
  Treatment = colData(dds_GFP_titers)$Treatment
)

print(rpoC_counts_GFP)


ggplot(rpoC_counts, aes(x = Treatment, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  theme_minimal() +
  labs(title = paste("Normalized expression of", gene),
       x = "Treatment",
       y = "DESeq2 normalized counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
###################CTL###################

norm_counts_CTL <- counts(dds_CTL_titers, normalized = TRUE)
head(norm_counts_CTL)
gene <- "rpoC"

if (!gene %in% rownames(norm_counts_CTL)) {
  stop("Gene not found — check ID format (case-sensitive?)")
}

rpoC_counts_CTL <- data.frame(
  Expression = norm_counts[gene, ],
  Sample = colnames(norm_counts_CTL),
  Treatment = colData(dds_CTL_titers)$Treatment
)

print(rpoC_counts_CTL)


norm_counts <- counts(dds_titers, normalized = TRUE)
rpoC_counts <- data.frame(
  Expression = norm_counts[gene, ],
  Sample = colnames(norm_counts_CTL),
  Treatment = colData(dds_titers)$Treatment
)


ggplot(rpoC_counts_CTL, aes(x = Treatment, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  theme_minimal() +
  labs(title = paste("Normalized expression of", gene),
       x = "Treatment",
       y = "DESeq2 normalized counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


rpoC_all <- rbind(rpoC_counts_CTL, rpoC_counts_GFP, rpoC_counts_yell)
library(ggplot2)

ggplot(rpoC_all, aes(x = Treatment, y = Expression)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  theme_minimal() +
  labs(
    title = "Normalized expression of rpoC across conditions",
    y = "DESeq2 normalized counts",
    x = "Treatment Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Pastel1")

rpoC_stats_yell <- res_yell_titers["rpoC", ]
pval_yell <- rpoC_stats_yell$pvalue
padj_yell <- rpoC_stats_yell$padj
rpoC_stats_GFP <- res_GFP_titers["rpoC", ]
pval_GFP <- rpoC_stats_GFP$pvalue
padj_GFP <- rpoC_stats_GFP$padj
rpoC_stats_CTL <- res_CTL_titers["rpoC", ]
pval_CTL <- rpoC_stats_CTL$pvalue
padj_CTL <- rpoC_stats_CTL$padj


sig_label_yell <- pval_to_star(pval_yell)
sig_label_GFP <- pval_to_star(pval_GFP)
sig_label_CTL <- pval_to_star(pval_CTL)


library(ggplot2)

ggplot(rpoC_counts_yell, aes(x = Treatment, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  theme_minimal() +
  labs(
    title = "Expression of rpoC in YELL (DESeq2 normalized)",
    subtitle = paste("DESeq2 p-value:", signif(pval_yell, 3)),
    y = "Normalized counts"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") +
  # Add significance star
  annotate("text", x = 1.5, y = max(rpoC_counts_yell$Expression) * 1.05,
           label = sig_label, size = 6)

library(ggplot2)

ggplot(rpoC_counts_GFP, aes(x = Treatment, y = Expression, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3) +
  theme_minimal() +
  labs(
    title = "Expression of rpoC in GFP (DESeq2 normalized)",
    subtitle = paste("DESeq2 p-value:", signif(pval_GFP, 3)),
    y = "Normalized counts"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2") 


#########################HEATMAP YELLOW UNT YELLOW HUM
# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

coldata_yell <- coldata[coldata$Treatment %in% c("YELL_UNT", "YELL_HUM"), ]
countdata_yell <- countdata[, rownames(coldata_yell)]
coldata_yell$titers <- scale(as.numeric(coldata_yell$titers)) 

dds_yell_titers <- DESeqDataSetFromMatrix(countData = countdata_yell,
                                    colData = coldata_yell,
                                    design = ~ titers+Replicate+Treatment)
dds_yell_titers <- dds_yell_titers[rowSums(counts(dds_yell_titers)) > 10, ]
dds_yell_titers <- dds_yell_titers[keep,]
dds_yell_titers <- DESeq(dds_yell_titers)

ntd_yell_titers <- normTransform(dds_yell_titers)
data_ntd_yell_titers<-assay(ntd_yell_titers)

top_genes_yell_titers <- head(order(rowMeans(data_ntd_yell_titers), decreasing = TRUE), 300)
top_data_yell_titers <- data_ntd_yell_titers[top_genes_yell_titers, ]

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
pheatmap(top_data_yell_titers_sub_ordered,
         scale = "row",
         annotation_row = annotation_ordered,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         fontsize_row = 5,
         main = "Heatmap of Annotated Genes Ordered by Category",
         color = my_palette)
#########################HEATMAP ALL SAMPLES ########################
# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

coldata$titers <- scale(as.numeric(coldata$titers)) 

dds_titers <- DESeqDataSetFromMatrix(countData = countdata,
                                       colData = coldata,
                                       design = ~ titers + Replicate+Treatment)

keep <- rowSums(counts(dds_titers)) > 10
dds_titers <- dds_titers[keep,]
dds_titers <- DESeq(dds_titers)

ntd_titers <- normTransform(dds_titers)
data_ntd_titers<-assay(ntd_titers)

top_genes_titers <- head(order(rowMeans(data_ntd_titers), decreasing = TRUE), 100)
top_data_titers <- data_ntd_titers[top_genes_titers, ]

annotations <- read.delim("annotation_Stammera.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

annotations$gene <- trimws(annotations$gene)

common_genes <- intersect(rownames(top_data_titers), annotations$gene)

top_data_titers_sub <- top_data_titers[common_genes, ]
annot_sub <- annotations[match(common_genes, annotations$gene), ]

annotation_row <- data.frame(Category = annot_sub$Categories)
rownames(annotation_row) <- annot_sub$gene
ordering <- order(annotation_row$Category)
top_data_titers_sub_ordered <- top_data_titers_sub[ordering, ]
annotation_ordered <- annotation_row[ordering, , drop = FALSE]

annotation_ordered$Category <- factor(annotation_ordered$Category, 
                                      levels = sort(unique(annotation_ordered$Category)))
sorted_categories <- levels(annotation_ordered$Category)
cat_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(sorted_categories))
names(cat_colors) <- sorted_categories
ann_colors <- list(Category = cat_colors)

my_palette <- colorRampPalette(c("#a6dba0","#01665e","#8c510a"))(100)
pheatmap(top_data_titers_sub_ordered,
         scale = "row",
         annotation_row = annotation_ordered,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         fontsize_row = 5,
         main = "Heatmap of Annotated Genes Ordered by Category",
         color = my_palette)


#################YELL UNT YELL HUM REMOVING YELL HUM R3##################################
coldata_yell <- coldata[coldata$Treatment %in% c("YELL_UNT", "YELL_HUM"), ]
coldata_yell <- coldata_yell[!rownames(coldata_yell) %in% "YELL_HUM_R3", ]
countdata_yell <- countdata[, rownames(coldata_yell)]
coldata_yell$titers <- scale(as.numeric(coldata_yell$titers))
dds_yell_titers <- DESeqDataSetFromMatrix(
  countData = countdata_yell,
  colData = coldata_yell,
  design = ~ titers+Treatment
)

keep <- rowSums(counts(dds_yell_titers)) > 10
dds_yell_titers <- dds_yell_titers[keep, ]
dds_yell_titers <- DESeq(dds_yell_titers)
res_yell_titers <- results(dds_yell_titers, contrast = c("Treatment", "YELL_HUM", "YELL_UNT"))
summary(res_yell_titers)


##################AVERAGE##########
treatment_info <- coldata_yell$Treatment
names(treatment_info) <- rownames(coldata_yell)
# Collapse columns by treatment, taking rowMeans for each treatment group
avg_expr_by_treatment <- sapply(unique(treatment_info), function(trt) {
  rowMeans(data_ntd_yell_titers[, treatment_info == trt])
})
top_genes <- head(order(rowMeans(avg_expr_by_treatment), decreasing = TRUE), 80)
avg_expr_top <- avg_expr_by_treatment[top_genes, ]
avg_expr_top <- avg_expr_top[complete.cases(avg_expr_top), ]

pheatmap(avg_expr_top,
         scale = "row",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         main = "Average Expression Across Replicates")


dds_yell_noBact
dds_yell_titers

normalized_counts_yell<-counts(dds_yell, normalized=TRUE)
write.csv(as.data.frame(normalized_counts_yell),file="normalized_counts_DESeq2_yell.csv")
normalized_counts_yell_titers<-counts(dds_yell_titers, normalized=TRUE)
write.csv(as.data.frame(normalized_counts_yell_titers),file="normalized_counts_DESeq2_yell_titers.csv")
