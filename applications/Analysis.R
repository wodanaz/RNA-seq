#################################################################
#################################################################
##### THIS IS AN ANALYSIS OF GENE EXPRESSION FOR      ###########
#####           TIME SERIES EXPERIMENT                ###########
#####                                                 ###########
#####   ( RESULTS FOR SNAIL/TWIST (invertebrate))     ###########
#####                                                 ###########
#################################################################


setwd("~/Desktop/snail")


library(arrayQualityMetrics)
library(DESeq)

library(VennDiagram)
library(genefilter)
library(gplots)
library(vegan)
library(plotrix)
library(rgl)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(reshape2)

source("uniHeatmap2.R") 
source("multiplot.R") 


annotations=read.delim("Lvar_master_top_hits_final.edited.txt",header=TRUE, row.names=1) 
dim(annotations)
head(annotations)

counts=read.delim("snail.twist_counts_final.txt",header=TRUE, row.names=1) 
counts$X.1 <- NULL
dim(counts)
head(counts)

colData <- colsplit(colnames(counts), "_", c("Sample", "Test", "Snumber","Lane"))
colData
time_series =c(rep("8hpf",9),rep("10hpf", 9), rep("12hpf",9), rep("14hpf",9),rep("16hpf",9))
length(time_series)

KnockOut=c(rep("Control",3),rep("Snail", 3), rep("Twist", 3),
           rep("Control",3),rep("Snail", 3), rep("Twist", 3),
           rep("Control",3),rep("Snail", 3), rep("Twist", 3),
           rep("Control",3),rep("Snail", 3), rep("Twist", 3),
           rep("Control",3),rep("Snail", 3), rep("Twist", 3))


colData$Condition <- KnockOut
colData$Time <- time_series




colData
colnames(counts) <- paste(colData$Condition, colData$Time, colData$Sample,  sep="_")

dim(counts)
head(counts)




#####use arrayQualityMetrics to find sample outliers
dev.off()
#-----------------
# creating counts data object, calculating size factors
rnaDesign <- colData
getwd()
adat= newCountDataSet(counts, rnaDesign)
str(rnaDesign)
str(counts)
adat=estimateSizeFactors(adat)
sizeFactors(adat)
# Plot size factors
plot(density(sizeFactors(adat)))
plot(sort(sizeFactors(adat)))

# Which are the outlayers???

sort(sizeFactors(adat))


# range from 0.5 to 2, 4-fold difference max - not too bad

#----------------
# quality metrics: see if we need to discard some samples as outliers
# For READs normalization!!!

adat=estimateDispersions(adat,method="blind")
vst=varianceStabilizingTransformation(adat)
# arrayQualityMetrics(expressionset = vsdBlind, outdir = "Report_for_vsdBlind",force =TRUE)
# prepare to wait about 5 min...
arrayQualityMetrics(vst, intgroup=c("Condition","Time"),outdir = "Report_for_vsd1",force=T)



# RESULTS OF THIS ANALYSIS INDICATE THERE WERE NO OUTLIERS FOUND


# if ribosomal genes are a problem, identify a list and remove it
ribo=read.delim("ribo_genes.tab", header=  F, stringsAsFactors = FALSE)  
head(ribo)

vector_ribo <- ribo$V1
cts.noRibo <- counts[!rownames(counts) %in% vector_ribo, ]
dim(counts)
dim(cts.noRibo)



##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
detach("package:DESeq", unload = TRUE)
detach("package:arrayQualityMetrics", unload = TRUE)
dev.off()

library(DESeq2)


means=apply(cts.noRibo,1,mean,na.rm=T)
table(means>1)
cts_clean=cts.noRibo[means>1,]
dim(cts_clean)
dim(colData)


write.table(cts_clean, file = "counts_clean_snail_twist.txt",  append = FALSE, row.names = T, quote=F,  sep="\t")
write.table(colData, file = "metadata.txt",  append = FALSE, row.names = T, quote=F,  sep="\t")



ddsFull <- DESeqDataSetFromMatrix(
  countData = cts_clean, 
  colData = colData, 
  design = ~ Condition+Time)
ddsFull



pheatmap(cor(assay(ddsFull)),border_color=NA, module_rows = TRUE,
         module_cols = TRUE)

ddsFull$Condition <- relevel(ddsFull$Condition, ref="Control")
ddsFull$Time <- relevel(ddsFull$Time, ref="8hpf")


levels(ddsFull$Condition)
levels(ddsFull$Time)

###########################################################################################################
###########################################################################################################
# TIME  is the variable of Interest
dds_time <- DESeq(ddsFull, test="LRT", reduced = ~   Condition  )                     # Time is the variable of Interest
vsd_time  <- vst(dds_time, nsub = 5000)

adonis(t(assay(vsd_time))~Condition*Time,data=colData,permutations = 1000,method="manhattan")  

###########################################################################################################
###########################################################################################################
# Condition  is the variable of Interest
dds_knockout <- DESeq(ddsFull, test="LRT", reduced = ~ Time  )                     # Knockout is the variable of Interest

vsd_knockout <- getVarianceStabilizedData(dds_knockout) 

adonis(t(assay(vsd_knockout))~Condition*Time,data=colData,permutations = 1000,method="manhattan")  

################################
### VISUALIZING DATA  

getwd()

dim(dds_time)
dim(cts.norib)


rv <- rowVars(assay(vsd_time))
head(rv)
length(rv)
ntop <- 1000
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                   length(rv)))]
select
length(select)
pca <- prcomp(t(assay(vsd_time)[select, ]))
pca
percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVar

dim(colData)

d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], colData = colData, name = colnames(vsd_time))
d


percentVar[1]
percentVar[2]
percentVar[3]
percentVar[4]
percentVar[5]
percentVar[6]

library(ggplot2)
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "colData.Condition")) + 
  #geom_point(size = 2) +
  geom_text(label=d$name, size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  coord_fixed() + theme_bw()

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "colData.Time")) + 
  geom_point(size = 2) +
  #geom_text(label=d$name, size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  coord_fixed() + theme_bw()


ggplot(data = d, aes_string(x = "PC1", y = "PC2",  shape = "colData.Time",color = "colData.Condition")  ) + 
  geom_point(size = 2) + scale_shape_manual(values=1:nlevels(d$colData.time)) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  coord_fixed() + theme_bw()



ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "colData.Time", shape = "colData.Condition")  ) + 
  geom_point(size = 2) +
  #geom_text(label=d$name, size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  coord_fixed() + theme_bw()



ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "colData.Condition")) + 
  #geom_point(size = 0.5) +
  geom_text(label=d$name, size = 3) +
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
  coord_fixed() + theme_bw()


ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "colData.Condition")) + 
  #geom_point(size = 0.5) +
  geom_text(label=d$name, size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
  coord_fixed() + theme_bw()


