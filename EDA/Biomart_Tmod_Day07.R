#### HEADER #### 

rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")
library(biomaRt)
library(dplyr)
library(tmod)

############## BIOMART ##############

convert.mmGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  return(genesV2)
  muosex <- unique(genesV2[, 1])
}

##### DAY 7 ######

#### Comparison A: I_NS X NI_NS

A7 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_7_Non stimulated A.csv")
A7$table$genes <- as.character(A7$table$genes)
mmGenes_A7 <- as.data.frame(A7$table$genes)
hsGenesA7 <- convert.mmGeneList(mmGenes_A7)

A7_H <- inner_join(hsGenesA7, as.data.frame(A7), by=c("MGI.symbol"="genes")) # put mice genes on table
A7_H <- subset(A7_H, !duplicated(A7_H$HGNC.symbol))
index_A7 <- order(A7_H$FDR)
A7_H <- A7_H[index_A7,] # order by FDR
A7_H <- A7_H[, c(2, 4, 8)]
View(A7_H)
nrow(A7_H) # 11239
write.csv2(A7_H, "Results/EdgeR/Biomart/A7_H.csv", row.names = FALSE)

# tmod

martA7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_A7_H.csv")
index_A7 <- order(martA7$FDR)
martA7 <- martA7[index_A7,]
A7_c <- tmodCERNOtest(martA7$HGNC.symbol)
A7_LI <- list("Unstimulated"=A7_c)
write.csv2(A7_LI, file = "Results/EdgeR/Tmod/A7_LI.csv", row.names = FALSE)
TabelaA7 <- read.csv2("Results/EdgeR/Tmod/A7_LI.csv")
View(TabelaA7)

#### Comparison B: NI_T X NI_NS

B7 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_7_Non infected B.csv")
B7$table$genes <- as.character(B7$table$genes)
mmGenes_B7 <- as.data.frame(B7$table$genes)
hsGenesB7 <- convert.mmGeneList(mmGenes_B7)

B7_H <- inner_join(hsGenesB7, as.data.frame(B7), by=c("MGI.symbol"="genes")) # put mice genes on table
B7_H <- subset(B7_H, !duplicated(B7_H$HGNC.symbol))
index_B7 <- order(B7_H$FDR)
B7_H <- B7_H[index_B7,] # order by FDR
View(B7_H)
B7_H <- B7_H[, c(2, 4, 8)]
View(B7_H)
nrow(B7_H) # 11239
write.csv2(B7_H, "Results/EdgeR/Biomart/B7_H.csv", row.names = FALSE)

# tmod

martB7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_B7_H.csv")
index_B7 <- order(martB7$FDR)
martB7 <- martB7[index_B7,]
B7_c <- tmodCERNOtest(martB7$HGNC.symbol)
B7_LI<- list("Non_Infected"=B7_c)
write.csv2(B7_LI, file = "Results/EdgeR/Tmod/B7_LI.csv", row.names = F)
TabelaB7 <- read.csv2('Results/EdgeR/Tmod/B7_LI.csv')
View(TabelaB7)

#### Comparison C: I_T X I_NS

C7 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_7_Infected C.csv")
C7$table$genes <- as.character(C7$table$genes)
mmGenes_C7 <- as.data.frame(C7$table$genes)
hsGenesC7 <- convert.mmGeneList(mmGenes_C7)

C7_H <- inner_join(hsGenesC7, as.data.frame(C7), by=c("MGI.symbol"="genes")) # put mice genes on table
C7_H <- subset(C7_H, !duplicated(C7_H$HGNC.symbol))
index_C7 <- order(C7_H$FDR)
C7_H <- C7_H[index_C7,] # order by FDR
C7_H <- C7_H[, c(2, 4, 8)]
View(C7_H)
nrow(C7_H) # 11239
write.csv2(C7_H, "Results/EdgeR/Biomart/C7_H.csv", row.names = FALSE)

# tmod

martC7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_C7_H.csv")
index_C7 <- order(martC7$FDR)
martC7 <- martC7[index_C7,]
C7_c <- tmodCERNOtest(martC7$HGNC.symbol)
C7_LI<- list("Infected"=C7_c)
write.csv2(C7_LI, file = "Results/EdgeR/Tmod/C7_LI.csv", row.names = F)
TabelaC7 <- read.csv2('Results/EdgeR/Tmod/C7_LI.csv')
View(TabelaC7)

#### Comparison D: I_T X NI_T

D7 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_7_Stimulated D.csv")
D7$table$genes <- as.character(D7$table$genes)
mmGenes_D7 <- as.data.frame(D7$table$genes)
hsGenesD7 <- convert.mmGeneList(mmGenes_D7)

D7_H <- inner_join(hsGenesD7, as.data.frame(D7), by=c("MGI.symbol"="genes")) # put mice genes on table
D7_H <- subset(D7_H, !duplicated(D7_H$HGNC.symbol))
index_D7 <- order(D7_H$FDR)
D7_H <- D7_H[index_D7,] # order by FDR
D7_H <- D7_H[, c(2, 4, 8)]
View(D7_H)
nrow(D7_H) # 11239
write.csv2(D7_H, "Results/EdgeR/Biomart/D7_H.csv", row.names = FALSE)

# tmod

martD7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_D7_H.csv")
index_D7 <- order(martD7$FDR)
martD7 <- martD7[index_D7,]
D7_c <- tmodCERNOtest(martD7$HGNC.symbol)
D7_LI<- list("Stimulated"=D7_c)
write.csv2(D7_LI, file = "Results/EdgeR/Tmod/D7_LI.csv", row.names = F)
TabelaD7 <- read.csv2('Results/EdgeR/Tmod/D7_LI.csv')
View(TabelaD7)

#### Comparison E: I_T X NI_NS

E7 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_7_Infected_Non infected E.csv")
E7$table$genes <- as.character(E7$table$genes)
mmGenes_E7 <- as.data.frame(E7$table$genes)
hsGenesE7 <- convert.mmGeneList(mmGenes_E7)

E7_H <- inner_join(hsGenesE7, as.data.frame(E7), by=c("MGI.symbol"="genes")) # put mice genes on table
E7_H <- subset(E7_H, !duplicated(E7_H$HGNC.symbol))
index_E7 <- order(E7_H$FDR)
E7_H <- E7_H[index_E7,] # order by FDR
E7_H <- E7_H[, c(2, 4, 8)]
View(E7_H)
nrow(E7_H) # 11215
write.csv2(E7_H, "Results/EdgeR/Biomart/E7_H.csv", row.names = FALSE)

# tmod

martE7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_E7_H.csv")
index_E7 <- order(martE7$FDR)
martE7 <- martE7[index_E7,]
E7_c <- tmodCERNOtest(martE7$HGNC.symbol)
E7_LI<- list("Infected/Non_Infected"=E7_c)
write.csv2(E7_LI, file = "Results/EdgeR/Tmod/E7_LI.csv", row.names = FALSE)
TabelaE7 <- read.csv2('Results/EdgeR/Tmod/E7_LI.csv')
View(TabelaE7)

########## tmod pies ##########
# Create Pies
pieA<- tmodDecideTests(A7_H$HGNC.symbol, lfc=A7_H$logFC, pval = A7_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieB<- tmodDecideTests(B7_H$HGNC.symbol, lfc=B7_H$logFC, pval = B7_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieC<- tmodDecideTests(C7_H$HGNC.symbol, lfc=C7_H$logFC, pval = C7_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieD<- tmodDecideTests(D7_H$HGNC.symbol, lfc=D7_H$logFC, pval = D7_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieE<- tmodDecideTests(E7_H$HGNC.symbol, lfc=E7_H$logFC, pval = E7_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

# transform in data frame
pieA<- as.data.frame(pieA)
pieB<- as.data.frame(pieB)
pieC<- as.data.frame(pieC)
pieD<- as.data.frame(pieD)
pieE<- as.data.frame(pieE)

# Remove the X if there is any
colnames(pieA) <- gsub("X.*\\.", "", colnames(pieA))
colnames(pieB) <- gsub("X.*\\.", "", colnames(pieB))
colnames(pieC) <- gsub("X.*\\.", "", colnames(pieC))
colnames(pieD) <- gsub("X.*\\.", "", colnames(pieD))
colnames(pieE) <- gsub("X.*\\.", "", colnames(pieE))

pieC<- list("Infected"=pieC)
pieD<- list("Stimulated"= pieD)
pieA<- list("Unstimulated"= pieA)
pieB<- list("Not Infected"=pieB)
pieE<- list("Infected/Non infected"=pieE)

A7_LI<- list("Unstimulated"=A7_c)
B7_LI<- list("Not Infected"=B7_c)
C7_LI<- list("Infected"=C7_c)
D7_LI<- list("Stimulated"=D7_c)
E7_LI<- list("Infected/Non infected"=E7_c)

# Create the panels

# I_NS NI_NS vs I_T NI_NS
panel_list_AE <- c(A7_LI, E7_LI)
is.list(panel_list_AE)
pie_list_AE <- c(pieA, pieE)

tmodPanelPlot(panel_list_AE, pval.thr = 10^-2, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_AE, pie.style = "pie")


# NI_T NI_NS vs I_T I_NS

panel_list_BC <- c(B7_LI, C7_LI)
is.list(panel_list_BC)
pie_list_BC <- c(pieB, pieC)


tmodPanelPlot(panel_list_BC, pval.thr = 10^-4, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_BC, pie.style = "pie")


# NI_T NI_NS vs I_T I_NS

panel_list_DE <- c(D7_LI, E7_LI)
is.list(panel_list_DE)
pie_list_DE <- c(pieD, pieE)

tmodPanelPlot(panel_list_DE, pval.thr = 10^-2, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_DE, pie.style = "pie")


