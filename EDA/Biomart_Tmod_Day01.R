#### HEADER #### 

rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")
library(biomaRt)
library(tmod)
library(dplyr)

############## BIOMART ##############

convert.mmGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  return(genesV2)
  muosex <- unique(genesV2[, 1])
}

  ##### DAY 1 ######

    #### Comparison A: I_NS X NI_NS

A1 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_1_Non stimulated A.csv")
A1$table$genes <- as.character(A1$table$genes)
mmGenes_A1 <- as.data.frame(A1)
hsGenesA1 <- convert.mmGeneList(mmGenes_A1$genes)

A1_H <- inner_join(hsGenesA1, as.data.frame(A1), by=c("MGI.symbol"="genes")) # put mice genes on table
#A1_H <- subset(A1_H, !duplicated(A1_H$HGNC.symbol))
index_A1 <- order(A1_H$FDR)
A1_H <- A1_H[index_A1,] # order by FDR
A1_H <- A1_H[, c(2, 4, 8)]
View(A1_H)
nrow(A1_H) # before remova duplicated 12187, after 11215
write.csv2(A1_H, "Results/EdgeR/Biomart/A1_H.csv", row.names = FALSE)

## LOOP
# tmod

martA1 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_A1_H.csv")
View(martA1)
index_A1 <- order(martA1$FDR)
martA1 <- martA1[index_A1,]
A_c <- tmodCERNOtest(martA1$HGNC.symbol)
A_LI <- list("Unstimulated"=A_c)
write.csv2(A_LI, file = "Results/EdgeR/Tmod/A1_LI.csv", row.names = FALSE)
TabelaA1 <- read.csv2("Results/EdgeR/Tmod/A1_LI.csv")
View(TabelaA1)

    #### Comparison B: NI_T X NI_NS

B1 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_1_Non infected B.csv")
B1$table$genes <- as.character(B1$table$genes)
mmGenes_B1 <- as.data.frame(B1$table$genes)
hsGenesB1 <- convert.mmGeneList(mmGenes_B1)

B1_H <- inner_join(hsGenesB1, as.data.frame(B1), by=c("MGI.symbol"="genes")) # put mice genes on table
B1_H <- subset(B1_H, !duplicated(B1_H$HGNC.symbol))
index_B1 <- order(B1_H$FDR)
B1_H <- B1_H[index_B1,] # order by FDR
#View(B1_H)
B1_H <- B1_H[, c(2, 4, 8)]
View(B1_H)
nrow(B1_H) # 11215
write.csv2(B1_H, "Results/EdgeR/Biomart/B1_H.csv", row.names = FALSE)

## LOOP
# tmod

martB1 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_B1_H.csv")
index_B1 <- order(martB1$FDR)
martB1 <- martB1[index_B1,]
B_c <- tmodCERNOtest(martB1$HGNC.symbol)
B_LI <- list("Unstimulated"=B_c)
write.csv2(B_LI, file = "Results/EdgeR/Tmod/B1_LI.csv", row.names = FALSE)
TabelaB1 <- read.csv2("Results/EdgeR/Tmod/B1_LI.csv")
View(TabelaB1)

    #### Comparison C: I_T X I_NS

C1 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_1_Infected C.csv")
C1$table$genes <- as.character(C1$table$genes)
mmGenes_C1 <- as.data.frame(C1$table$genes)
hsGenesC1 <- convert.mmGeneList(mmGenes_C1)

C1_H <- inner_join(hsGenesC1, as.data.frame(C1), by=c("MGI.symbol"="genes")) # put mice genes on table
C1_H <- subset(C1_H, !duplicated(C1_H$HGNC.symbol))
index_C1 <- order(C1_H$FDR)
C1_H <- C1_H[index_C1,] # order by FDR
C1_H <- C1_H[, c(2, 4, 8)]
View(C1_H)
nrow(C1_H) # 11215
write.csv2(C1_H, "Results/EdgeR/Biomart/C1_H.csv", row.names = FALSE)

## LOOP
# tmod

martC1 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_C1_H.csv")
index_C1 <- order(martC1$FDR)
martC1 <- martC1[index_C1,]
C_c <- tmodCERNOtest(martC1$HGNC.symbol)
C1_LI<- list("Infected"=C_c)
write.csv2(C1_LI, file = "Results/EdgeR/Tmod/C1_LI.csv", row.names = F)
TabelaC1 <- read.csv2('Results/EdgeR/Tmod/C1_LI.csv')
View(TabelaC1)

    #### Comparison D: I_T X NI_T

D1 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_1_Stimulated D.csv")
D1$table$genes <- as.character(D1$table$genes)
mmGenes_D1 <- as.data.frame(D1$table$genes)
hsGenesD1 <- convert.mmGeneList(mmGenes_D1)

D1_H <- inner_join(hsGenesD1, as.data.frame(D1), by=c("MGI.symbol"="genes")) # put mice genes on table
D1_H <- subset(D1_H, !duplicated(D1_H$HGNC.symbol))
index_D1 <- order(D1_H$FDR)
D1_H <- D1_H[index_D1,] # order by FDR
D1_H <- D1_H[, c(2, 4, 8)]
View(D1_H)
nrow(D1_H) # 11215
write.csv2(D1_H, "Results/EdgeR/Biomart/D1_H.csv", row.names = FALSE)

## LOOP
# tmod

martD1 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_D1_H.csv")
index_D1 <- order(martD1$FDR)
martD1 <- martD1[index_D1,]
D_c <- tmodCERNOtest(martD1$HGNC.symbol)
D1_LI<- list("Stimulated"=D_c)
write.csv2(D1_LI, file = "Results/EdgeR/Tmod/D1_LI.csv", row.names = F)
TabelaD1 <- read.csv2('Results/EdgeR/Tmod/D1_LI.csv')
View(TabelaD1)

    #### Comparison E: I_T X NI_NS

E1 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_1_Infected_Non infected E.csv")
E1$table$genes <- as.character(E1$table$genes)
mmGenes_E1 <- as.data.frame(E1$table$genes)
hsGenesE1 <- convert.mmGeneList(mmGenes_E1)

E1_H <- inner_join(hsGenesE1, as.data.frame(E1), by=c("MGI.symbol"="genes")) # put mice genes on table
E1_H <- subset(E1_H, !duplicated(E1_H$HGNC.symbol))
index_E1 <- order(E1_H$FDR)
E1_H <- E1_H[index_E1,] # order by FDR
E1_H <- E1_H[, c(2, 4, 8)]
View(E1_H)
nrow(E1_H) # 11215
write.csv2(E1_H, "Results/EdgeR/Biomart/E1_H.csv", row.names = FALSE)

# tmod

martE1 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_E1_H.csv")
index_E1 <- order(martE1$FDR)
martE1 <- martE1[index_E1,]
E_c <- tmodCERNOtest(martE1$HGNC.symbol)
E1_LI<- list("Infected/Non_Infected"=E_c)
write.csv2(E1_LI, file = "Results/EdgeR/Tmod/E1_LI.csv", row.names = FALSE)
TabelaE1 <- read.csv2('Results/EdgeR/Tmod/E1_LI.csv')
View(TabelaE1)

########## tmod pies ##########
# Create Pies
pieA1<- tmodDecideTests(martA1$HGNC.symbol, lfc=martA1$logFC, pval = martA1$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieB1<- tmodDecideTests(martB1$HGNC.symbol, lfc=martB1$logFC, pval = martB1$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieC1<- tmodDecideTests(martC1$HGNC.symbol, lfc=martC1$logFC, pval = martC1$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieD1<- tmodDecideTests(martD1$HGNC.symbol, lfc=martD1$logFC, pval = martD1$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieE1<- tmodDecideTests(martE1$HGNC.symbol, lfc=martE1$logFC, pval = martE1$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

# transform in data frame
pieA<- as.data.frame(pieA1)
pieB<- as.data.frame(pieB1)
pieC<- as.data.frame(pieC1)
pieD<- as.data.frame(pieD1)
pieE<- as.data.frame(pieE1)
View(pieE)

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

A1_LI<- list("Unstimulated"=A_c)
B1_LI<- list("Not Infected"=B_c)
C1_LI<- list("Infected"=C_c)
D1_LI<- list("Stimulated"=D_c)
E1_LI<- list("Infected/Non infected"=E_c)

# Create the panels

# I_NS NI_NS vs I_T NI_NS
panel_list_A_E <- c(A1_LI, E1_LI)
is.list(panel_list_A_E)
pie_list_A_E <- c(pieA, pieE)

tmodPanelPlot(panel_list_A_E, pval.thr = 10^-2, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.60, clust = "qval", pie =pie_list_A_E, pie.style = "pie")


# NI_T NI_NS vs I_T I_NS

panel_listBC <- c(B1_LI, C1_LI)
is.list(panel_listBC)
pie_listBC <- c(pieB, pieC)


tmodPanelPlot(panel_listBC, pval.thr = 10^-4, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_listBC, pie.style = "pie")


# NI_T NI_NS vs I_T I_NS

panel_list_DE <- c(D1_LI, E1_LI)
is.list(panel_list_DE)
pie_list_DE <- c(pieD, pieE)

tmodPanelPlot(panel_list_DE, pval.thr = 10^-2, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_DE, pie.style = "pie")

