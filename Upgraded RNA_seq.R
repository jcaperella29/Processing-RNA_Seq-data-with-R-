library(MLSeq)


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
nTest <- ceiling(ncol(gse) * 0.3)
ind <- sample(ncol(gse), nTest, FALSE)
data.train <- gse[ ,-ind]
data.test <- gse[ ,ind] 



#differnetial expression analysis starting from counts and table of phenotype information



train_countdata <- round(assays(data.train)[["counts"]])
test_countdata<-round(assays(data.test)[["counts"]])


train_coldata <- colData(data.train)
test_coldata<-colData(data.test)


ddsMat <- DESeqDataSetFromMatrix(countData = train_countdata,
                                 colData = train_coldata,
                                 design = ~ dex)

#filtering data
smallestGroupSize <- 4
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- ddsMat[keep,]


dds <- DESeq(ddsMat)



res <- results(dds)
res<-data.frame(res)
sorted_res<-res[order(res$padj),]
sorted_res$genes<-row.names(sorted_res)
filter_res<-filter(sorted_res, sorted_res$padj<=.05)


count_df<-as.data.frame(train_countdata)
count_df<-t(count_df)
count_df<-as.data.frame(count_df)
Phenotype<-train_coldata$dex

count_df$phenotype<-Phenotype


count_df<-count_df[,names(count_df) %in% c(filter_res$genes,"phenotype"),]





model <- randomForest( as.factor(count_df$phenotype) ~.,data=count_df,importance=TRUE)
importance_df<-data.frame(importance(model))
importance_df$genenames<-row.names(importance_df)
importance_df<-importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = TRUE),]
Top_20_genes<-importance_df$genenames[1:20]
TOP_df<-count_df[,names(count_df) %in% c(Top_20_genes,"phenotype"),]
model_RF <- randomForest( as.factor(TOP_df$phenotype) ~.,data=TOP_df)

predictions <- predict(model_RF, t(test_countdata), 
                      type="response")

confuse <- confusionMatrix(data=as.factor(predictions), reference = as.factor(test_coldata$dex))
con_matrix<-data.frame(confuse$byClass)

write.csv(con_matrix,"C:/Users/ccape/OneDrive/Documents/result20genes.csv")


write.table(Top_20_genes,"C:/Users/ccape/Downloads/20hits.txt")

#orginal code is from https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html.
#I add lines to  set up train and test sets, perform feature selection using random forest ,lastly I used the 20 best genes to build a classifier.