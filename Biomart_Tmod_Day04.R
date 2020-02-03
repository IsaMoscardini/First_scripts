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

##### DAY 4 ######

#### Comparison A: I_NS X NI_NS

A4 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_4_Non stimulated A.csv")
A4$table$genes <- as.character(A4$table$genes)
mmGenes_A4 <- as.data.frame(A4$table$genes)
hsGenesA4 <- convert.mmGeneList(mmGenes_A4)

A4_H <- inner_join(hsGenesA4, as.data.frame(A4), by=c("MGI.symbol"="genes")) # put mice genes on table
A4_H <- subset(A4_H, !duplicated(A4_H$HGNC.symbol))
index_A4 <- order(A4_H$FDR)
A4_H <- A4_H[index_A4,] # order by FDR
A4_H <- A4_H[, c(2, 4, 8)]
View(A4_H)
nrow(A4_H) # 11260
write.csv2(A4_H, "Results/EdgeR/Biomart/A4_H.csv", row.names = FALSE)

# tmod

martA4 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_A4_H.csv")
index_A4 <- order(martA4$FDR)
martA4 <- martA4[index_A4,]
A4_c <- tmodCERNOtest(martA4$HGNC.symbol)
A4_LI <- list("Unstimulated"=A4_c)
write.csv2(A4_LI, file = "Results/EdgeR/Tmod/A4_LI.csv", row.names = FALSE)
TabelaA4 <- read.csv2("Results/EdgeR/Tmod/A4_LI.csv")
View(TabelaA4)

#### Comparison B: NI_T X NI_NS

B4 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_4_Non infected B.csv")
B4$table$genes <- as.character(B4$table$genes)
mmGenes_B4 <- as.data.frame(B4$table$genes)
hsGenesB4 <- convert.mmGeneList(mmGenes_B4)

B4_H <- inner_join(hsGenesB4, as.data.frame(B4), by=c("MGI.symbol"="genes")) # put mice genes on table
B4_H <- subset(B4_H, !duplicated(B4_H$HGNC.symbol))
index_B4 <- order(B4_H$FDR)
B4_H <- B4_H[index_B4,] # order by FDR
B4_H <- B4_H[, c(2, 4, 8)]
View(B4_H)
nrow(B4_H) # 11260
write.csv2(B4_H, "Results/EdgeR/Biomart/B4_H.csv", row.names = FALSE)

# tmod

martB4 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_B4_H.csv")
index_B4 <- order(martB4$FDR)
martB4 <- martB4[index_B4,]
B4_c <- tmodCERNOtest(martB4$HGNC.symbol)
B4_LI<- list("Non_Infected"=B4_c)
write.csv2(B4_LI, file = "Results/EdgeR/Tmod/B4_LI.csv", row.names = F)
TabelaB4 <- read.csv2('Results/EdgeR/Tmod/B4_LI.csv')
View(TabelaB4)

#### Comparison C: I_T X I_NS

C4 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_4_Infected C.csv")
C4$table$genes <- as.character(C4$table$genes)
mmGenes_C4 <- as.data.frame(C4$table$genes)
hsGenesC4 <- convert.mmGeneList(mmGenes_C4)

C4_H <- inner_join(hsGenesC4, as.data.frame(C4), by=c("MGI.symbol"="genes")) # put mice genes on table
C4_H <- subset(C4_H, !duplicated(C4_H$HGNC.symbol))
index_C4 <- order(C4_H$FDR)
C4_H <- C4_H[index_C4,] # order by FDR
C4_H <- C4_H[, c(2, 4, 8)]
View(C4_H)
nrow(C4_H) # 11260
write.csv2(C4_H, "Results/EdgeR/Biomart/C4_H.csv", row.names = FALSE)

# tmod

martC4 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_C4_H.csv")
index_C4 <- order(martC4$FDR)
martC4 <- martC4[index_C4,]
C4_c <- tmodCERNOtest(martC4$HGNC.symbol)
C4_LI<- list("Infected"=C4_c)
write.csv2(C4_LI, file = "Results/EdgeR/Tmod/C4_LI.csv", row.names = F)
TabelaC4 <- read.csv2('Results/EdgeR/Tmod/C4_LI.csv')
View(TabelaC4)

#### Comparison D: I_T X NI_T

D4 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_4_Stimulated D.csv")
D4$table$genes <- as.character(D4$table$genes)
mmGenes_D4 <- as.data.frame(D4$table$genes)
hsGenesD4 <- convert.mmGeneList(mmGenes_D4)

D4_H <- inner_join(hsGenesD4, as.data.frame(D4), by=c("MGI.symbol"="genes")) # put mice genes on table
D4_H <- subset(D4_H, !duplicated(D4_H$HGNC.symbol))
index_D4 <- order(D4_H$FDR)
D4_H <- D4_H[index_4,] # order by FDR
D4_H <- D4_H[, c(2, 4, 8)]
View(D4_H)
nrow(D4_H) # 11215
write.csv2(D4_H, "Results/EdgeR/Biomart/D4_H.csv", row.names = FALSE)

# tmod

martD4 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_D4_H.csv")
index_D4 <- order(martD4$FDR)
martD4 <- martD4[index_D4,]
D4_c <- tmodCERNOtest(martD4$HGNC.symbol)
D4_LI<- list("Stimulated"=D4_c)
write.csv2(D4_LI, file = "Results/EdgeR/Tmod/D4_LI.csv", row.names = F)
TabelaD4 <- read.csv2('Results/EdgeR/Tmod/D4_LI.csv')
View(TabelaD4)

#### Comparison E: I_T X NI_NS

E4 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_4_Infected_Non infected E.csv")
E4$table$genes <- as.character(E4$table$genes)
mmGenes_E4 <- as.data.frame(E4$table$genes)
hsGenesE4 <- convert.mmGeneList(mmGenes_E4)

E4_H <- inner_join(hsGenesE4, as.data.frame(E4), by=c("MGI.symbol"="genes")) # put mice genes on table
E4_H <- subset(E4_H, !duplicated(E4_H$HGNC.symbol))
index_E4 <- order(E4_H$FDR)
E4_H <- E4_H[index_E4,] # order by FDR
E4_H <- E4_H[, c(2, 4, 8)]
View(E4_H)
nrow(E4_H) # 11260
write.csv2(E4_H, "Results/EdgeR/Biomart/E4_H.csv", row.names = FALSE)

# tmod

martE4 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_E4_H.csv")
index_E4 <- order(martE4$FDR)
martE4 <- martE4[index_E4,]
E4_c <- tmodCERNOtest(martE4$HGNC.symbol)
E4_LI<- list("Infected/Non_Infected"=E4_c)
write.csv2(E4_LI, file = "Results/EdgeR/Tmod/E4_LI.csv", row.names = FALSE)
TabelaE4 <- read.csv2('Results/EdgeR/Tmod/E4_LI.csv')
View(TabelaE4)

########## tmod pies ##########
# Create Pies
pieA<- tmodDecideTests(A4_H$HGNC.symbol, lfc=A4_H$logFC, pval = A4_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieB<- tmodDecideTests(B4_H$HGNC.symbol, lfc=B4_H$logFC, pval = B4_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieC<- tmodDecideTests(C4_H$HGNC.symbol, lfc=C4_H$logFC, pval = C4_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieD<- tmodDecideTests(D4_H$HGNC.symbol, lfc=D4_H$logFC, pval = D4_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieE<- tmodDecideTests(E4_H$HGNC.symbol, lfc=E4_H$logFC, pval = E4_H$FDR, 
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

A4_LI<- list("Unstimulated"=A4_c)
B4_LI<- list("Not Infected"=B4_c)
C4_LI<- list("Infected"=C4_c)
D4_LI<- list("Stimulated"=D4_c)
E4_LI<- list("Infected/Non infected"=E4_c)

# Create the panels

# I_NS NI_NS vs I_T NI_NS
panel_list_AE <- c(A4_LI, E4_LI)
is.list(panel_list_AE)
pie_list_AE <- c(pieA, pieE)

tmodPanelPlot(panel_list_AE, pval.thr = 10^-2, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_AE, pie.style = "pie")


# NI_T NI_NS vs I_T I_NS

panel_list_BC <- c(B4_LI, C4_LI)
is.list(panel_list_BC)
pie_list_BC <- c(pieB, pieC)


tmodPanelPlot(panel_list_BC, pval.thr = 10^-4, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_BC, pie.style = "pie")


# NI_T NI_NS vs I_T I_NS

panel_list_DE <- c(D1_LI, E1_LI)
is.list(panel_list_DE)
pie_list_DE <- c(pieD, pieE)

tmodPanelPlot(panel_list_DE, pval.thr = 10^-4, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_DE)

# NOTHING TO PLOT
