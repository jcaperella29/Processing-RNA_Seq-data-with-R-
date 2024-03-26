library("airway")
library(dplyr)
library("ggplot2")
library("tximeta")
library("pheatmap")
library("RColorBrewer")
library("DESeq2")
library(randomForest)
library(varImp)
library(caret)

dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))

csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:2,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
file.exists(coldata$files)
coldata

se <- tximeta(coldata)

dim(se)

head(rownames(se))
gse <- summarizeToGene(se)
#dimension are reduced and row names are gene id(Emsembl ids)
dim(gse)


head(rownames(gse))
data(gse)

#isolating counts
assayNames(gse)
head(assay(gse), 3)


colSums(assay(gse))
#The rowRanges, when printed, shows the ranges for the first five and last five genes:


rowRanges(gse)
seqinfo(rowRanges(gse))

colData(gse)
#renaming labels0
gse$cell <- gse$donor
gse$dex <- gse$condition
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")

#differnetial expression analysis starting from counts and table of phenotype information



countdata <- round(assays(gse)[["counts"]])
head(countdata, 3)
coldata <- colData(gse)
dds <- DESeqDataSet(gse, design = ~ cell + dex)


ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)

#filtering data
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]
#vst transformatin
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)


dds <- estimateSizeFactors(dds)

df <- bind_rows(
  data.frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

sampleDists <- dist(t(assay(vsd)))
sampleDists

#examining simalarity between sampls using heatmap
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#PCA

plotPCA(vsd, intgroup = c("dex", "cell"))

pcaData <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))


ggplot(pcaData, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


dds <- DESeq(dds)

res <- results(dds)
res<-data.frame(res)
sorted_res<-res[order(res$padj),]
sorted_res$genes<-row.names(sorted_res)
filter_res<-filter(sorted_res, sorted_res$padj<=5e-8)


count_df<-data.frame(countdata)
count_df$genes<-row.names(count_df)
count_df<-count_df[count_df$genes %in% filter_res$genes,]

count_df<-subset(count_df,select=-genes)

trans_counts<-t(count_df)
count_pheno_df<-data.frame(trans_counts)
count_pheno_df$phenotype<-gse$condition
model <- randomForest( as.factor(count_pheno_df$phenotype)~.,data=count_pheno_df,importance=TRUE)
importance_df<-data.frame(importance(model))
importance_df$genenames<-row.names(importance_df)
importance_df<-importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE),]
Top_20_genes<-importance_df$genenames[1:20]

write.table(Top_20_genes,"path choosen by user)

#orginal code is from https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html.
#I added lines to narrow down hits  among the differentially expressed genes using feature selection via Random Forest
#As well as to output the 20 most important genes as chosen  by feature selection.
