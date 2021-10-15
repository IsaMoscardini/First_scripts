#### HEADER ####

rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")

#install.packages("tmod")
library(dplyr)
library(Glimma)
library(edgeR)
library(data.table)
library(tmod)
library(tidyr)
library(DESeq2)

countdata <- read.csv2('data/Mouse_counts_20191024.csv', row.names = 1)
descript <- read.csv2('data/Mouse_desc_20191024.csv')

# Create group
descript$Group <- factor(paste(descript$Vivo, descript$Vitro, sep="_"))

# Check samples and row/column names
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

# Check
setdiff(colnames(countdata), rownames(descript_clean)) #  "X1_4_NS"  "X1_6_NS"  "X10_2_NS"
setdiff(rownames(descript_clean), colnames(countdata)) # 0

# Remove outliers from countdata table
common.samples <- intersect(colnames(countdata), rownames(descript_clean))
countdata_clean <- countdata[, common.samples]
identical(colnames(countdata_clean), rownames(descript_clean))



#############------------ DAY 1 ------------############

# Day 1

day_one <- descript_clean[descript_clean$Time.point %in% c(0,1),]
#View(day_one)

common.samples1 <- intersect(colnames(countdata), rownames(day_one))
length(common.samples1)  # 22
count_one <- countdata[, common.samples1] 
#View(count_one)

identical(rownames(day_one), colnames(count_one))

# create DEGlist
y <- DGEList(counts = count_one, genes = row.names(count_one), group= day_one$Group)

# filrter
keep <- rowSums(cpm(y)>2)>=4
y_1 <- y[keep,]
dim(y_1) # 14068    22

# normalize
y_1 <- calcNormFactors(y_1) # normaliza as reads MAAAAAS... ver vst

# Create design table to perform GLM (general linear model)
design <- model.matrix(~0 + day_one$Group)
row.names(design) <- day_one$Sample.Name
colnames(design) <- levels(day_one$Group)

# Dispersion and Fit estimation
fit <- estimateDisp(y_1, design)
qlfit <- glmQLFit(fit, design)

# Contrasts
my.contrasts <- makeContrasts(A= I_NS - NI_NS,
                              B= NI_T - NI_NS,
                              C= I_T - I_NS,
                              D= I_T - NI_T,
                              E= I_T - NI_NS,
                              levels=design)

#dir.create('Results')
#dir.create('Results/DEGs')

# Differential analysis
ql<- glmQLFTest(qlfit, contrast=my.contrasts[,"A"]) 
A<- topTags(ql, n=Inf)
write.csv2(A, "Results/DEGs/TP_1_Non stimulated A.csv")
summary (dup<- decideTestsDGE(ql))
#        1*I_NS -1*NI_NS
# Down              5195
# NotSig            6598
# Up                2275
# TOTAL 14068

ql_B <- glmQLFTest(qlfit, contrast=my.contrasts[,"B"]) 
B <- topTags(ql_B, n=Inf)
write.csv2(B, "Results/DEGs/TP_1_Non infected B.csv")
summary (dup<- decideTestsDGE(ql_B))
#        -1*NI_NS 1*NI_T
# Down              5338
# NotSig            6412
# Up                2318
# TOTAL 14068

ql_C <- glmQLFTest(qlfit, contrast=my.contrasts[,"C"]) 
C <- topTags(ql_C, n=Inf)
write.csv2(C, "Results/DEGs/TP_1_Infected C.csv")
summary (dup<- decideTestsDGE(ql_C))
#        -1*I_NS 1*I_T
# Down             840
# NotSig         12515
# Up               713
# TOTAL 14068

ql_D <- glmQLFTest(qlfit, contrast=my.contrasts[,"D"]) 
D <- topTags(ql_D, n=Inf)
write.csv2(D, "Results/DEGs/TP_1_Stimulated D.csv")
summary (dup<- decideTestsDGE(ql_D))
#        1*I_T -1*NI_T
# Down             680
# NotSig         12768
# Up               620
# TOTAL 14068

ql_E<- glmQLFTest(qlfit, contrast=my.contrasts[,"E"]) 
E<- topTags(ql_E, n=Inf)
write.csv2(E, "Results/DEGs/TP_1_Infected_Non infected E.csv")
summary (dup<- decideTestsDGE(ql_E))
#        1*I_T -1*NI_NS
# Down             5682
# NotSig           5541
# Up               2845
# TOTAL 14068

#grafici 
# glMDPlot(ql2,  status=dup, main="Infetti",
#          counts=y.1$counts,  groups=tab_desc$Group, transform = T,
#          search.by="genes", launch=TRUE)
# 
# plotSmear(ql3, de.tags = rownames(D)[which(D$table$FDR <0.05)], 
#           cex = 0.5, main="Stimulated (D)", cex.main=1.8)
# 
# abline(h = c(-0.5, 0.5), col = "blue")


#############------------ DAY 2 ------------############

# Day 2

day_two <- descript_clean[descript_clean$Time.point %in% c(0, 2) ,]
View(day_two)

common.samples2 <- intersect(colnames(countdata), rownames(day_two))
length(common.samples2) # 22
count_two <- countdata[, common.samples2] 
View(count_two)

identical(rownames(day_two), colnames(count_two))

# create DEGlist
y <- DGEList(counts = count_two, genes = row.names(count_two), group= day_two$Group)
#class(y)
# filrter
keep <- rowSums(cpm(y)>2)>=4
y_2 <- y[keep,]
dim(y_2) #   14042    22

# normalize
y_2 <- calcNormFactors(y_2) 

# Create design table to perform GLM (general linear model)
design <- model.matrix(~0 + day_two$Group)
row.names(design) <- day_two$Sample.Name
colnames(design) <- levels(day_two$Group)

# Dispersion and Fit estimation
fit <- estimateDisp(y_2, design)
qlfit <- glmQLFit(fit, design)

# Contrasts
my.contrasts <- makeContrasts(A= I_NS - NI_NS,
                              B= NI_T - NI_NS,
                              C= I_T - I_NS,
                              D= I_T - NI_T,
                              E= I_T - NI_NS,
                              levels=design)


# Differential analysis
ql_2_A <- glmQLFTest(qlfit, contrast=my.contrasts[,"A"]) 
A_2 <- topTags(ql_2_A, n=Inf)
write.csv2(A_2, "Results/EdgeR/DEGs/DEGs_each_timepoint/TP_2_Non stimulated A.csv")
summary (dup<- decideTestsDGE(ql_2_A))
#        1*I_NS -1*NI_NS
# Down              5347
# NotSig            5946
# Up                2749
# TOTAL 14042,

ql_2_B <- glmQLFTest(qlfit, contrast=my.contrasts[,"B"]) 
B_2 <- topTags(ql_2_B, n=Inf)
write.csv2(B_2, "Results/DEGs/TP_2_Non infected B.csv")
summary (dup<- decideTestsDGE(ql_2_B))
#        -1*NI_NS 1*NI_T
# Down              5339
# NotSig            6198
# Up                2505
# TOTAL 14042,

ql_2_C <- glmQLFTest(qlfit, contrast=my.contrasts[,"C"]) 
C_2 <- topTags(ql_2_C, n=Inf)
write.csv2(C_2, "Results/DEGs/TP_2_Infected C.csv")
summary (dup<- decideTestsDGE(ql_2_C))
#        -1*I_NS 1*I_T
# Down            1416
# NotSig         11409
# Up              1217
# TOTAL 14042,

ql_2_D <- glmQLFTest(qlfit, contrast=my.contrasts[,"D"]) 
D_2 <- topTags(ql_2_D, n=Inf)
write.csv2(D_2, "Results/DEGs/TP_2_Stimulated D.csv")
summary (dup<- decideTestsDGE(ql_2_D))
#        1*I_T -1*NI_T
# Down            1409
# NotSig         10948
# Up              1685
# TOTAL 14042,

ql_2_E<- glmQLFTest(qlfit, contrast=my.contrasts[,"E"]) 
E_2 <- topTags(ql_2_E, n=Inf)
write.csv2(E_2, "Results/DEGs/TP_2_Infected_Non infected E.csv")
summary (dup<- decideTestsDGE(ql_2_E))
#        1*I_T -1*NI_NS
# Down             5776
# NotSig           4866
# Up               3400
# TOTAL 14042,

#############------------ DAY 4 ------------############

# Day 4

day_four <- descript_clean[descript_clean$Time.point %in% c(0, 4) ,]
View(day_four)

common.samples4 <- intersect(colnames(countdata), rownames(day_four))
length(common.samples4) # 22
count_four <- countdata[, common.samples4] 
View(count_four)

identical(rownames(day_four), colnames(count_four))

# create DEGlist
y <- DGEList(counts = count_four, genes = row.names(count_four), group= day_four$Group)

# filrter
keep <- rowSums(cpm(y)>2)>=4
y_4 <- y[keep,]
dim(y_4) # 14082    22

# normalize
y_4 <- calcNormFactors(y_4)

# Create design table to perform GLM (general linear model)
design <- model.matrix(~0 + day_four$Group)
row.names(design) <- day_four$Sample.Name
colnames(design) <- levels(day_four$Group)

# Dispersion and Fit estimation
fit <- estimateDisp(y_4, design)
qlfit <- glmQLFit(fit, design)

# Contrasts
my.contrasts <- makeContrasts(A= I_NS - NI_NS,
                              B= NI_T - NI_NS,
                              C= I_T - I_NS,
                              D= I_T - NI_T,
                              E= I_T - NI_NS,
                              levels=design)


# Differential analysis
ql_4_A <- glmQLFTest(qlfit, contrast=my.contrasts[,"A"]) 
A_4 <- topTags(ql_4_A, n=Inf)
write.csv2(A_4, "Results/EdgeR/DEGs/DEGs_each_timepoint/TP_4_Non stimulated A.csv")
summary (dup<- decideTestsDGE(ql_4_A))
#        1*I_NS -1*NI_NS
# Down              4878
# NotSig            6588
# Up                2616
# TOTAL 14082

ql_4_B <- glmQLFTest(qlfit, contrast=my.contrasts[,"B"]) 
B_4 <- topTags(ql_4_B, n=Inf)
write.csv2(B_4, "Results/DEGs/TP_4_Non infected B.csv")
summary (dup<- decideTestsDGE(ql_4_B))
#       -1*NI_NS 1*NI_T
# Down             4965
# NotSig           7146
# Up               1971
# TOTAL 14082

ql_4_C <- glmQLFTest(qlfit, contrast=my.contrasts[,"C"]) 
C_4 <- topTags(ql_4_C, n=Inf)
write.csv2(C_4, "Results/DEGs/TP_4_Infected C.csv")
summary (dup<- decideTestsDGE(ql_4_C))
#        -1*I_NS 1*I_T
# Down             982
# NotSig         12280
# Up               820
# TOTAL 14082

ql_4_D <- glmQLFTest(qlfit, contrast=my.contrasts[,"D"]) 
D_4 <- topTags(ql_4_D, n=Inf)
write.csv2(D_4, "Results/DEGs/TP_4_Stimulated D.csv")
summary (dup<- decideTestsDGE(ql_4_D))
#        1*I_T -1*NI_T
# Down            1507
# NotSig         11007
# Up              1568
# TOTAL 14082

ql_4_E<- glmQLFTest(qlfit, contrast=my.contrasts[,"E"]) 
E_4 <- topTags(ql_4_E, n=Inf)
write.csv2(E_4, "Results/DEGs/TP_4_Infected_Non infected E.csv")
summary (dup<- decideTestsDGE(ql_4_E))
#        1*I_T -1*NI_NS
# Down             5586
# NotSig           5393
# Up               3103
# TOTAL 14082

#############------------ DAY 7 ------------############

# Day 7

day_seven <- descript_clean[descript_clean$Time.point %in% c(0, 7) ,]
View(day_seven)

common.samples7 <- intersect(colnames(countdata), rownames(day_seven))
length(common.samples7) # 21
count_seven <- countdata[, common.samples7] 
View(count_seven)

identical(rownames(day_seven), colnames(count_seven))

# create DEGlist
y <- DGEList(counts = count_seven, genes = row.names(count_seven), group= day_seven$Group)

# filrter
keep <- rowSums(cpm(y)>2)>=4
y_7 <- y[keep,]
dim(y_7) # 14066    21

# normalize
y_7 <- calcNormFactors(y_7)

# Create design table to perform GLM (general linear model)
design <- model.matrix(~0 + day_seven$Group)
row.names(design) <- day_seven$Sample.Name
colnames(design) <- levels(day_seven$Group)

# Dispersion and Fit estimation
fit <- estimateDisp(y_7, design)
qlfit <- glmQLFit(fit, design)

# Contrasts
my.contrasts <- makeContrasts(A= I_NS - NI_NS,
                              B= NI_T - NI_NS,
                              C= I_T - I_NS,
                              D= I_T - NI_T,
                              E= I_T - NI_NS,
                              levels=design)


# Differential analysis
ql_7_A <- glmQLFTest(qlfit, contrast=my.contrasts[,"A"]) 
A_7 <- topTags(ql_7_A, n=Inf)
write.csv2(A_7, "Results/DEGs/TP_7_Non stimulated A.csv")
summary (dup<- decideTestsDGE(ql_7_A))
#        1*I_NS -1*NI_NS
# Down              4970
# NotSig            7300
# Up                1796
# TOTAL 14066

ql_7_B <- glmQLFTest(qlfit, contrast=my.contrasts[,"B"]) 
B_7 <- topTags(ql_7_B, n=Inf)
write.csv2(B_7, "Results/DEGs/TP_7_Non infected B.csv")
summary (dup<- decideTestsDGE(ql_7_B))
#       -1*NI_NS 1*NI_T
# Down              5042
# NotSig            7786
# Up                1238
# TOTAL 14066

ql_7_C <- glmQLFTest(qlfit, contrast=my.contrasts[,"C"]) 
C_7 <- topTags(ql_7_C, n=Inf)
write.csv2(C_7, "Results/DEGs/TP_7_Infected C.csv")
summary (dup<- decideTestsDGE(ql_7_C))
#        -1*I_NS 1*I_T
# Down             535
# NotSig         13133
# Up               398
# TOTAL 14066

ql_7_D <- glmQLFTest(qlfit, contrast=my.contrasts[,"D"]) 
D_7 <- topTags(ql_7_D, n=Inf)
write.csv2(D_7, "Results/DEGs/TP_7_Stimulated D.csv")
summary (dup<- decideTestsDGE(ql_7_D))
#        1*I_T -1*NI_T
# Down             788
# NotSig         12313
# Up               965
# TOTAL 14066

ql_7_E<- glmQLFTest(qlfit, contrast=my.contrasts[,"E"]) 
E_7 <- topTags(ql_7_E, n=Inf)
write.csv2(E_7, "Results/DEGs/TP_7_Infected_Non infected E.csv")
summary (dup<- decideTestsDGE(ql_7_E))
#        1*I_T -1*NI_NS
# Down             5418
# NotSig           6233
# Up               2415
# TOTAL 14066






