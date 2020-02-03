## HEADER
rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")

# PACKAGES
library(Glimma)
library(edgeR)

countdata <- read.csv2('data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('data/Mouse_desc_20191024.csv')
View(descript)

# create group 
descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))
View(head(descript))

rownames(descript) <- paste0("X", descript$Sample.Name)

# ter certeza que sao as mesmas amostras

common.samples <- intersect(colnames(countdata), rownames(descript))
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))


# create DEGlist
y <- DGEList(counts = countdata, genes = row.names(countdata), group= descript$Group)
#class(y)
# filrter
keep <- rowSums(cpm(y)>2)>=6
y.1 <- y[keep,]
dim(y.1) # 13153    60

# normalize
y.1 <- calcNormFactors(y.1) # normaliza as reads ver vst

#BiocManager::install('vegan')
library(ggplot2)
library(vegan)

cp <- log2(cpm(y.1)+1) # normalizou de vdd
sample.var <- apply(cp, 2, var)
sample.var <- sort(sample.var)
#barplot(sample.var, las=2)    # NEEDS SCALATION! PLEASE ATTENTION! - Do quality control if it is possible


boxplot(x = y.1)


boxplot(x = cp, use.cols = TRUE)

#View(cp)
pca<-prcomp(t(cp))
#pca<-prcomp(t(cp), center = FALSE, scale. = TRUE)  # transpor e escalar
View(t(cp))

pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

percentVar<- round(100* percentVar)

group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

ggplot(pcapanel, aes(PC2, PC3, color=descript$Group, shape = as.character(descript$Time.point)))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ labs(shape="Time Point", group = "group")+
  scale_color_manual(values = group.colors)+ geom_text(label=descript$Sample.Name, nudge_x = 0.9, nudge_y = 1)+
  xlab(paste0("PC2: ",percentVar[2],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance"))
 
# FIND PERTURBED GENES IN CONTROL SAMPLES

pert_genes <- data.frame(sort(abs(pca$rotation[,"PC1"]), decreasing=TRUE)[1:5000])
list(pert_genes)
genes_to_remove <- rownames(pert_genes)

View(genes_to_remove)
setdiff(genes_to_remove, rownames(countdata)) # 0


# DO THE ANALYSIS AGAIN, WITHOUT THESE GENES

nrow(countdata) # 23930
countdata_after <- countdata[setdiff(rownames(countdata), genes_to_remove),]
nrow(countdata_after) 
dim(countdata_after) # 18930    60

y <- DGEList(counts = countdata_after, genes = row.names(countdata_after), group= descript$Group)

keep_after <- rowSums(cpm(y)>2)>=6
y.1 <- y[keep_after,]
dim(y.1) # 13103    60

# normalize
y.1 <- calcNormFactors(y.1) # normaliza as reads ver vst

cp <- log2(cpm(y.1)+1) # normalizou de vdd
pca<-prcomp(t(cp))
pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

percentVar<- round(100* percentVar)

group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

ggplot(pcapanel, aes(PC1, PC2, color=descript$Group, shape = as.character(descript$Time.point)))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ labs(shape="Time Point", group = "group")+
  scale_color_manual(values = group.colors)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



##### PCA without controls

View(descript)
descript_OUT <- descript[descript$Group != "NI_NS",]

common.samples <- intersect(colnames(countdata), rownames(descript_OUT))
countdata_out <- countdata[, common.samples]
View(countdata_out)

y <- DGEList(counts = countdata_out, genes = row.names(countdata_out), group= descript_OUT$Group)
#class(y)
# filrter
keep <- rowSums(cpm(y)>2)>=6
y.1 <- y[keep,]
dim(y.1) # 13153    60

# normalize
y.1 <- calcNormFactors(y.1) # normaliza as reads ver vst

#BiocManager::install('vegan')
library(ggplot2)
library(vegan)

cp <- log2(cpm(y.1)+1) # normalizou de vdd
#sample.var <- apply(cp, 2, var)
#sample.var <- sort(sample.var)
#barplot(sample.var, las=2)    # NEEDS SCALATION! PLEASE ATTENTION! - Do quality control if it is possible

boxplot(x = y.1)

boxplot(x = cp, use.cols = TRUE)

#View(cp)
pca<-prcomp(t(cp))
#pca<-prcomp(t(cp), center = FALSE, scale. = TRUE)  # transpor e escalar
View(t(cp))

pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

percentVar<- round(100* percentVar)

group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

ggplot(pcapanel, aes(PC1, PC2, color=descript_OUT$Group, shape = as.character(descript_OUT$Time.point)))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ labs(shape="Time Point", group = "group")+
  scale_color_manual(values = group.colors)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



############################## PCA SCRIPT ENDS HERE!!!!!! ################################
  
  
# Signal to noise ratio - QUALITY CONTROL - Not part of the PCA script
  
gene_var <- apply(cp, 1, var)
gene_mean <- rowMeans(cp)
gene_mean[1:5]
gene_var[1:5]
identical(names(gene_var), names(gene_mean))

plot(gene_mean, gene_var)

# DESeq 
library(DESeq2)
BiocManager::install('vsn')
library(vsn)
rownames(descript)[1:5]
colnames(countdata)[1:5]
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = descript, design = ~Group )

keep <- rowSums(counts(dds)>2)>=4
dds <- dds[keep,]
dim(dds)

dds <- DESeq(dds)

norm_data <- vst(dds, blind = F)

norm_data <- assay(dds)
norm_data[1:5, 1:5]

class(norm_data)
meanSdPlot(cp)


###################
countdata <- read.csv2('data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('data/Mouse_desc_20191024.csv')

# create group 
descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))
rownames(descript) <- paste0("X", descript$Sample.Name)

# ter certeza que sao as mesmas amostras

common.samples <- intersect(colnames(countdata), rownames(descript))
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))

off <- c("X1_4_NS", "X1_6_NS", "X10_2_NS")

desc_off <- descript[setdiff(rownames(descript), off),]
common.samples <- intersect(colnames(countdata), rownames(desc_off))
countdata_off <- countdata[, common.samples]

# y <- DGEList(counts = countdata_off, genes = row.names(countdata_off), group= desc_off$Group)
# #class(y)
# # filrter
# keep <- rowSums(cpm(y)>2)>=4
# y.1 <- y[keep,]
# dim(y.1) # 14299    57

library(DESeq2)

countdata_off <- as.matrix(countdata_off)

count_vst <- varianceStabilizingTransformation(countdata_off)
boxplot(count_vst)




############# without controls
rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")

library(Glimma)
library(edgeR)
library(ggplot2)

group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

countdata <- read.csv2('data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('data/Mouse_desc_20191024.csv')

descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))
View(descript)

rownames(descript) <- paste0("X", descript$Sample.Name)
common.samples <- intersect(colnames(countdata), rownames(descript))
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))

desc_out <- descript[which(descript$Group != "NI_NS"),]
common.samples <- intersect(colnames(countdata), rownames(desc_out))

setdiff(colnames(countdata), rownames(desc_out))
setdiff(rownames(desc_out), colnames(countdata))

count_out <- countdata[, common.samples]
View(count_out)

y <- DGEList(counts = count_out, genes = row.names(count_out), group= desc_out$Group)
class(y)
#filrter
keep <- rowSums(cpm(y)>2)>=6
y.1 <- y[keep,]
dim(y.1) # 14299    57

y.1 <- calcNormFactors(y.1) # normaliza as reads ver vst

#BiocManager::install('vegan')
library(ggplot2)
library(vegan)

cp <- log2(cpm(y.1)+1) # normalizou de vdd
pca<-prcomp(t(cp))
View(t(cp))

pcapanel<- as.data.frame(pca$x)
percentVar <- pca$sdev^2/sum(pca$sdev^2)

percentVar<- round(100* percentVar)

group.colors <- c(NI_NS = "#66C2A5", NI_T = "#8DA0CB", I_NS ="#FC8D62", I_T = "#E78AC3")

Group <- desc_out$Group
Samples <- desc_out$Sample.Name

ggplot(pcapanel, aes(PC1, PC2, color= Group, shape = as.character(desc_out$Time.point)))+
  geom_point(size = 5)+ theme_minimal()+ theme(legend.position = "top")+ labs(shape="Time Point", group = "group")+
  scale_color_manual(values = group.colors)+ geom_text(aes(label= Samples),hjust=0, vjust=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

###### test
cp <- t(cp)

trans.pca <- pca(cp, ncomp = 10, center = T, scale = T)

plotIndiv(trans.pca, group= desc_out$Group, legend = TRUE, title = 'PCA')

biplot(trans.pca, cex = 0.7,
               xlabs = paste(desc_out$Group, 1:nrow(cp)))






