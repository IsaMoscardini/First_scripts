# BiocManager::install("rJava")
# Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8.0_251/')

rm(list = ls())
options(stringsAsFactors = F)
setwd("PROJECTS/Spleen/")
library(S4Vectors)
library(rJava)
library(survival)
library(DaMiRseq)
library(lattice)

countdata <- read.csv2('Data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('Data/Mouse_desc_20191024.csv')
descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))

rownames(descript) <- paste0("X", descript$Sample.Name)
common.samples <- intersect(colnames(countdata), rownames(descript))
table(colnames(countdata) == rownames(descript))
setdiff(colnames(countdata), rownames(descript)) #"X6_2_.NS"
setdiff(rownames(descript), colnames(countdata)) #"X6_2_NS"
colnames(countdata) <- gsub("X6_2_.NS", "X6_2_NS", colnames(countdata))
identical(colnames(countdata), rownames(descript))

# Remove outliers from descript table
remove <- c("X1_4_NS", "X1_6_NS", "X10_2_NS")
descript_clean <- descript[setdiff(rownames(descript), remove),]
nrow(descript_clean)
descript_clean$Sample.Name <- rownames(descript_clean)
descript_clean$Barcode.ID <- NULL
descript_clean$Mapped.Reads <- NULL
descript_clean$Targets.Detected <- NULL
descript_clean$Valid.Reads <- NULL
colnames(descript_clean) <- c("Sample", "Chip", "Vitro", "Vivo", "TP", "class")
View(descript_clean)

# Remove outliers from countdata table
common.samples <- intersect(colnames(countdata), rownames(descript_clean))
countdata_clean <- countdata[, common.samples]
identical(colnames(countdata_clean), rownames(descript_clean))
View(countdata_clean)


######## DaMiR com controles ########

colnames(descript) <- c("Barcode", "Sample", "Mapped_reads", "Valid_reads", "targets_detected", "Chip", "Vitro", "Vivo", "TP", "class")
View(descript)
descript[] <- lapply(descript[], as.factor)
descript$Valid_reads <- gsub("%", "", descript$Valid_reads)
descript$targets_detected <- gsub("%", "", descript$targets_detected)
descript <- as.data.frame(descript)
descript$class <- as.factor(descript$class)
descript$Valid_reads <- as.numeric(descript$Valid_reads)
descript$targets_detected <- as.numeric(descript$targets_detected)
descript$Barcode <- NULL
descript$Mapped_reads <- NULL
descript$targets_detected <- NULL
descript$Valid_reads <- NULL
str(descript)
View(descript)

# Criar Summarized Experiment object
SE_all <- DaMiR.makeSE(x = countdata, y = descript) # 23930 Features
SE= SE_all[,which(SE_all$Vitro == "NS")]
table(SE$Vitro)
SE$Vitro <- NULL

# Filtering and applying VST normalization
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "no")             # 10776 Features remained

# Quality check on the correlation of samples
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.9)

sv <- DaMiR.SV(data_filt)
DaMiR.corrplot(sv, colData(data_filt), sig.level = 0.01)
dev.off()
#pdf("Figures/Quality_Check_DaMiRseq.pdf")
DaMiR.Allplot(data_filt, colData(data_filt))
#dev.off()

# x <- countdata[1:10,2]
# targets_detected <- function(x) sum(x > 0)/length(x) * 100
# targets_detected()
# tg <- apply(countdata, 2, targets_detected)
# tg[1:10]
# plot(as.numeric(descript$Mapped_reads), descript$targets_detected, col = as.numeric(descript$class))
# class(descript$targets_detected)
# identical(names(tg), rownames(descript))
# View(descript)

data_adjust <- DaMiR.SVadjust(data_filt, sv, n.sv=7)
assay(data_adjust[c(1:5), c(1:5, 21:25)])

pdf("Figures/Quality_Check_DaMiRseq_pos_QC.pdf")
DaMiR.Allplot(data_adjust, colData(data_adjust))
dev.off()

set.seed(12345)
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)

data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.4)
data_reduced <- DaMiR.FReduct(data_reduced$data)

DaMiR.MDSplot(data_reduced, df)

df.importance <- DaMiR.FSort(data_reduced, df)


######## DaMiR sem outliers ########

#colnames(descript_clean) <- c("Barcode", "Sample", "Mapped_reads", "Valid_reads", "targets_detected", "Chip", "Vitro", "Vivo", "TP", "class")
View(descript_clean)
descript_clean[] <- lapply(descript_clean[], as.factor)
descript_clean$Valid_reads <- gsub("%", "", descript_clean$Valid_reads)
descript_clean$targets_detected <- gsub("%", "", descript_clean$targets_detected)
descript_clean <- as.data.frame(descript_clean)
descript_clean$class <- as.factor(descript_clean$class)
descript_clean$Valid_reads <- as.numeric(descript_clean$Valid_reads)
descript_clean$targets_detected <- as.numeric(descript_clean$targets_detected)
descript_clean$Barcode <- NULL
str(descript_clean)
View(descript_clean)


# Criar Summarized Experiment object
SE <- DaMiR.makeSE(x = countdata_clean, y = descript_clean) # 23930 Features

# Filtering and applying VST normalization
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "no")             # 10774 Features remained

# Quality check on the correlation of samples
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.5) # DEVERIA retirar  X1_1_NS X1_3_NS 

sv <- DaMiR.SV(data_filt)
DaMiR.corrplot(sv, colData(data_filt), sig.level = 0.01)

#pdf("Figures/Quality_Check_DaMiRseq_sem_outliers.pdf")
DaMiR.Allplot(data_filt, colData(data_filt))
#dev.off()

data_adjust <- DaMiR.SVadjust(data_filt, sv, n.sv=7)
assay(data_adjust[c(1:5), c(1:5, 21:25)])

#pdf("Figures/Quality_Check_DaMiRseq_sem_outlier_pos_QC.pdf")
DaMiR.Allplot(data_adjust, colData(data_adjust))
#dev.off()

set.seed(12345)
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)

data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.4)
data_reduced <- DaMiR.FReduct(data_reduced$data)

DaMiR.MDSplot(data_reduced, df)

df.importance <- DaMiR.FSort(data_reduced, df)
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,
                                 n.pred = 10)
DaMiR.Clustplot(selected_features$data, df)
