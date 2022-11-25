library(ggplot2)
library(grid)
library(gridExtra)
library(DESeq2)

# Set the working directory
directory <- "C:/TEMP/Bb_study"
setwd(directory)

# Set the prefix for each output file name
outputPrefix <- "C3Hmapping_Inf-vs-Uninf__"

sampleFiles<- c("RNA-07-CH-C1_S7_counts-C3H.txt",
                "RNA-08-CH-B1_S8_counts-C3H.txt",
                "RNA-09-CH-A1_S9_counts-C3H.txt",
                "RNA-10-CH-C2_S10_counts-C3H.txt",
                "RNA-11-CH-B2_S11_counts-C3H.txt",
                "RNA-12-CH-C2_S12_counts-C3H.txt"
)

sampleNames <- c("Bb-C3H_1","Bb-C3H_2","Bb-C3H_3","C3H_1","C3H_2","C3H_3")
sampleCondition <- c("C3H-infected","C3H-infected","C3H-infected","C3H-UNinfected","C3H-UNinfected", "C3H-UNinfected")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

treatments = c("C3H-infected","C3H-UNinfected")

##

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

dds <- DESeq(ddsHTSeq)
res <- results(dds)

res <- results(dds, contrast = c("condition","C3H-infected","C3H-UNinfected"))
resultsALL <- results(dds, contrast = c("condition","C3H-infected","C3H-UNinfected"))

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
jpeg("C:/TEMP/Bb_study/C3H_plotMA.jpeg", quality=100, width=960, height=960, bg="white")
plotMA(dds, xlab = "Mean of normalized counts", ylim=c(-4,4))
dev.off()

###### plot PCAs ##################
# Regularized log transformation 

rld <- rlogTransformation(dds)
jpeg("C:/TEMP/Bb_study/PCA-C3H.jpeg", quality=100, width=3048, height=1648, bg="white", res=300)
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

