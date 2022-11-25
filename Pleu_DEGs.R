library(ggplot2)
library(grid)
library(gridExtra)
library(DESeq2)

# Set the working directory
directory <- "C:/TEMP/Bb_study"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "Peromyscus-perPleu_Inf-vs-Uninf__"

sampleFiles<- c("RNA-01-PL-C1_S1_Pleu_counts.txt",
                "RNA-02-PL-B1_S2_Pleu_counts.txt",
                "RNA-03-PL-A1_S3_Pleu_counts.txt",
				        "RNA-04-PL-C2_S4_Pleu_counts.txt",
                "RNA-05-PL-B2_S5_Pleu_counts.txt",
                "RNA-06-PL-A2_S6_Pleu_counts.txt"
                )

sampleNames <- c("Bb-PL_1","Bb-PL_2","Bb-PL_3","PL_1","PL_2","PL_3")
sampleCondition <- c("PL-infected","PL-infected","PL-infected","PL-UNinfected","PL-UNinfected", "PL-UNinfected")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments = c("PL-infected","PL-UNinfected")

##

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

dds <- DESeq(ddsHTSeq)
res <- results(dds)

res <- results(dds, contrast = c("condition","PL-infected","PL-UNinfected"))
resultsALL <- results(dds, contrast = c("condition","PL-infected","PL-UNinfected"))

res= subset(res, padj<0.1)
res <- res[order(res$padj),]  # order results by padj value (most significant to least)

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-padj0.1.csv"))

res= subset(res, padj<0.05)
res <- res[order(res$padj),]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-padj0.05.csv"))

## 

resultsALL <- results(dds)
resultsALLsorted <- resultsALL[order(resultsALL$padj),]
resultsALLsorted <- merge(as.data.frame(resultsALLsorted), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resultsALLsorted)[1] <- 'gene'
write.csv(resultsALLsorted, file = paste0(outputPrefix, "-results-ALL-with-normalizedCounts.csv"))

##
jpeg("C:/TEMP/Bb_study/Pleu_plotMA.jpeg", quality=100, width=960, height=960, bg="white")
plotMA(dds, xlab = "Mean of normalized counts", ylim=c(-4,4))
dev.off()

###### plot PCAs ##################
# Regularized log transformation 

rld <- rlogTransformation(dds)
jpeg("C:/TEMP/Bb_study/PCA-Pleu.jpeg", quality=100, width=3048, height=1648, bg="white", res=300)
plot_pca <- plot(plotPCA(rld, ntop = 500, returnData = FALSE)) +geom_point(size=15)  # increase dot size 
plot_pca <- plot_pca + geom_text(aes_string(label = "sampleNames"), color = "black") 
print(plot_pca)
dev.off()

#rm(list = ls())   ## clean out data

####################################
####################################

# sessionInfo()
# R version 4.0.5 (2021-03-31)
#
# other attached packages:
#  [1] DESeq2_1.30.1               SummarizedExperiment_1.20.0 Biobase_2.50.0              MatrixGenerics_1.2.1       
#  [5] matrixStats_0.58.0          GenomicRanges_1.42.0        GenomeInfoDb_1.26.4         IRanges_2.24.1             
#  [9] S4Vectors_0.28.1            BiocGenerics_0.36.0         gridExtra_2.3               ggplot2_3.3.3              




