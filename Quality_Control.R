#### HEADER ####
rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")

#BiocManager::install("DEGreport")

library(DEGreport)
library(edgeR)
library(DESeq2)

# Normalized count table

countdata <- read.csv2('data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('data/Mouse_desc_20191024.csv')

descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))

rownames(descript) <- paste0("X", descript$Sample.Name)
common.samples <- intersect(colnames(countdata), rownames(descript))
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))
identical(colnames(countdata), rownames(descript)) # TRUE

y <- DGEList(counts = countdata, genes = row.names(countdata), group= descript$Group)
keep <- rowSums(cpm(y)>2)>=6
y.1 <- y[keep,]
dim(y.1) # 13153    60
y.1 <- calcNormFactors(y.1) 
cp <- log2(cpm(y.1)+1) 
dim(cp) # 13153    60

# cpp <- read.csv2('Intermediate/Normalized_Count_cp')
# View(cpp)
# dim(cpp) # 13153    61
# 
# rownames(cpp) <- cpp$X
# cpp$X <- NULL
# View(cpp)



#### DESeq2 #### 
browseVignettes("DESeq2")
View(descript)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = descript,
                              design = ~ Group)

dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds)

counts <- counts(dds, normalized = TRUE)
design <- as.data.frame(colData(dds))

#--- check if it is scaled/normalized
library(reshape2)
library(ggplot2)
exprs_data <- melt(data.frame(cp))
ggplot(exprs_data, aes(x = value, group = variable)) + geom_density(show.legend = F) + theme_minimal(base_size = 15) +
  labs(x = "Expression Value", y = 'Density')


#### BARPLOT ####
# library size
names <- as.character(colnames(countdata))
barplot(colSums(countdata) * 1e-6, ylab = "Library Size (millions)")
# mean and variance
sample.var <- apply(cp, 2, var)
sample.var <- sort(sample.var)
barplot(sample.var, las=2, main = "Mean - Variance")
        
##### HEATMAP #### 
# BEFORE REMOVING OUTLIERS
#BiocManager::install('ComplexHeatmap')
library(ComplexHeatmap)
heat <- dist(t(cp), method = "euclidean") # Euclidean distance matrix.
Heatmap(as.matrix(heat))

# AFTER REMOVING OUTLIERS
remove <- c("X1_4_NS", "X1_6_NS", "X10_2_NS")
desc_three <- descript[setdiff(rownames(descript), remove),]
common.samples <- intersect(colnames(countdata), rownames(desc_three))
table(colnames(countdata) == rownames(desc_three))
setdiff(colnames(countdata), rownames(desc_three))
setdiff(rownames(desc_three), colnames(countdata))
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))
setdiff(colnames(countdata), rownames(desc_both)) #  "X1_4_NS"  "X1_6_NS"  "X10_2_NS"
setdiff(rownames(desc_three), colnames(countdata)) # 0
common.samples <- intersect(colnames(countdata), rownames(desc_three))
count_out_three <- countdata[, common.samples]
identical(colnames(count_out_three), rownames(desc_three))
y <- DGEList(counts = count_out_three, genes = row.names(count_out_three), group= desc_three$Group)
keep <- rowSums(cpm(y)>2)>=4
y.1 <- y[keep,]
dim(y.1) # 14299    57
y.1 <- calcNormFactors(y.1) 
cp_norm <- log2(cpm(y.1)+1) 
dim(cp_norm) # 14299    57

heat <- dist(t(cp_norm), method = "euclidean") # Euclidean distance matrix.
Heatmap(as.matrix(heat))

##### BOXPLOT ####
ggplot(cp, aes(x=, y=len)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)


#### SIZE FACTOR ####
# not good for too much sample, did for the controls
View(countdata)
degCheckFactors(countdata[, c(1:4, 17, 43)])
# did not save


#### MEAN-VARIANCE PLOT #####
# P-value distribution gives an idea on how well you model is capturing the input data.
# Whether it could be some problem for some set of genes.
# We expect to have a flat distribution with peaks at 0 and 1. 
# Add the mean count information to check if any set of genes are enriched in any specific p-value range.
# Variation (dispersion) and average expression relationship shouldn't be a factor among the differentially expressed genes. 
# When plotting average mean and standard deviation, significant genes should be randomly distributed.
# In this case, it would be good to look at the ones that are totally outside the expected correlation.

degQC(counts, design[["group"]], pvalue = res[["pvalue"]])


#### MDP ####
#BiocManager::install('mdp')
library(mdp)
#View(descript)
descript_mdp <- descript[, c("Sample.Name", "Group")]
colnames(descript_mdp) <- c("Sample", "Class")
descript_mdp$Sample <- paste0("X", descript_mdp$Sample)
View(descript_mdp)
cp_mdp <- as.data.frame(cp)
# View(cp_mdp)
# expression data has gene names in the rows
# pheno data needs a Sample and Class column
identical(descript_mdp$Sample, colnames(cp_mdp))
head(descript_mdp)
cp_mdp[1:5,1:5]
class(descript_mdp)
mdp.results <- mdp(data=cp_mdp, pdata=descript_mdp, control_lab = "NI_NS")



#### tSNE ####
#install.packages("Rtsne")
library(Rtsne)
#BiocManager::install("M3C")



#### RLE - Relative Log Expression ####
#plotRLE(countdata, output_filepath = 'Images/RLE', plot_type_name = "RLE Count", plot_width = 12, plot_height = 7, plot_flag = T)






