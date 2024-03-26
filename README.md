# Processing-RNA_Seq-data-with-R-
An R script that reads in RNA_Seq data  , preprocesses it , then performs differential expression analysis and feature selection.
#First needed libraries are imported.
#Sample data from the "airways" package is read in.
#the counts and phenotype data are isolated.
#data is then filtered based on Group size and undergo VST transformation.
#then the distance between samples is measured and a visualized via a heatmap.
#Next PCA analysis is performed .
#Then diffenrential expression  analysis is performed.
#results of differential expression are then sorted by adjusted p-value and p-values less 5e-8 are kept.
# After a bit of rearrangize the genes that produced the low p-values are isolated from the  counts matrix
# A new dataframe is prepared with columns of the counts and the phenotype.
# Using a Random Forest model , hits are narrowed down via feature selection
# Finally the top 20 genes are output as a text file.

