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

##### DAY 2 ######

#### Comparison A: I_NS X NI_NS

A2 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_2_Non stimulated A.csv")
A2$table$genes <- as.character(A2$table$genes)
mmGenes_A2 <- as.data.frame(A2$table$genes)
hsGenesA2 <- convert.mmGeneList(mmGenes_A2)

A2_H <- inner_join(hsGenesA2, as.data.frame(A2), by=c("MGI.symbol"="genes")) # put mice genes on table
A2_H <- subset(A2_H, !duplicated(A2_H$HGNC.symbol))
index_A2 <- order(A2_H$FDR)
A2_H <- A2_H[index_A2,] # order by FDR
A2_H <- A2_H[, c(2, 4, 8)]
View(A2_H)
nrow(A2_H) # 11222
write.csv2(A2_H, "Results/EdgeR/Biomart/A2_H.csv", row.names = FALSE)

# tmod

martA2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_A2_H.csv")
index_A2 <- order(martA2$FDR)
martA2 <- martA2[index_A2,]
A2_c <- tmodCERNOtest(martA2$HGNC.symbol)
A2_LI <- list("Unstimulated"=A2_c)
write.csv2(A2_LI, file = "Results/EdgeR/Tmod/A2_LI.csv", row.names = FALSE)
TabelaA2 <- read.csv2("Results/EdgeR/Tmod/A2_LI.csv")
View(TabelaA2)

#### Comparison B: NI_T X NI_NS

B2 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_2_Non infected B.csv")
B2$table$genes <- as.character(B2$table$genes)
mmGenes_B2 <- as.data.frame(B2$table$genes)
hsGenesB2 <- convert.mmGeneList(mmGenes_B2)

B2_H <- inner_join(hsGenesB2, as.data.frame(B2), by=c("MGI.symbol"="genes")) # put mice genes on table
B2_H <- subset(B2_H, !duplicated(B2_H$HGNC.symbol))
index_B2 <- order(B2_H$FDR)
B2_H <- B2_H[index_B2,] # order by FDR
#View(B2_H)
B2_H <- B2_H[, c(2, 4, 8)]
View(B2_H)
nrow(B2_H) # 11222
write.csv2(B2_H, "Results/EdgeR/Biomart/B2_H.csv", row.names = FALSE)

# tmod

martB2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_B2_H.csv")
index_B2 <- order(martB2$FDR)
martB2 <- martB2[index_B2,]
B2_c <- tmodCERNOtest(martB2$HGNC.symbol)
B2_LI<- list("Non_Infected"=B2_c)
write.csv2(B2_LI, file = "Results/EdgeR/Tmod/B2_LI.csv", row.names = F)
TabelaB2 <- read.csv2('Results/EdgeR/Tmod/B2_LI.csv')
View(TabelaB2)

#### Comparison C: I_T X I_NS

C2 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_2_Infected C.csv")
C2$table$genes <- as.character(C2$table$genes)
mmGenes_C2 <- as.data.frame(C2$table$genes)
hsGenesC2 <- convert.mmGeneList(mmGenes_C2)

C2_H <- inner_join(hsGenesC2, as.data.frame(C2), by=c("MGI.symbol"="genes")) # put mice genes on table
C2_H <- subset(C2_H, !duplicated(C2_H$HGNC.symbol))
index_C2 <- order(C2_H$FDR)
C2_H <- C2_H[index_C2,] # order by FDR
C2_H <- C2_H[, c(2, 4, 8)]
View(C2_H)
nrow(C2_H) # 11222
write.csv2(C2_H, "Results/EdgeR/Biomart/C2_H.csv", row.names = FALSE)

# tmod

martC2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_C2_H.csv")
View(martC2)
index_C2 <- order(martC2$FDR)
martC2 <- martC2[index_C2,]
C2_c <- tmodCERNOtest(martC2$HGNC.symbol)
C2_LI<- list("Infected"=C2_c)
write.csv2(C2_LI, file = "Results/EdgeR/Tmod/C2_LI.csv", row.names = F)
TabelaC2 <- read.csv2('Results/EdgeR/Tmod/C2_LI.csv')
View(TabelaC2)

#### Comparison D: I_T X NI_T

D2 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_2_Stimulated D.csv")
D2$table$genes <- as.character(D2$table$genes)
mmGenes_D2 <- as.data.frame(D2$table$genes)
hsGenesD2 <- convert.mmGeneList(mmGenes_D2)

D2_H <- inner_join(hsGenesD2, as.data.frame(D2), by=c("MGI.symbol"="genes")) # put mice genes on table
D2_H <- subset(D2_H, !duplicated(D2_H$HGNC.symbol))
index_D2 <- order(D2_H$FDR)
D2_H <- D2_H[index_D2,] # order by FDR
D2_H <- D2_H[, c(2, 4, 8)]
View(D2_H)
nrow(D2_H) # 11222
write.csv2(D2_H, "Results/EdgeR/Biomart/D2_H.csv", row.names = FALSE)

# tmod

martD2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_D2_H.csv")
View(martD2)
index_D2 <- order(martD2$FDR)
martD2 <- martD2[index_D2,]
D2_c <- tmodCERNOtest(martD2$HGNC.symbol)
D2_LI<- list("Stimulated"=D2_c)
write.csv2(D2_LI, file = "Results/EdgeR/Tmod/D2_LI.csv", row.names = F)
TabelaD2 <- read.csv2('Results/EdgeR/Tmod/D2_LI.csv')
View(TabelaD2)

#### Comparison E: I_T X NI_NS

E2 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_2_Infected_Non infected E.csv")
E2$table$genes <- as.character(E2$table$genes)
mmGenes_E2 <- as.data.frame(E2$table$genes)
hsGenesE2 <- convert.mmGeneList(mmGenes_E2)

E2_H <- inner_join(hsGenesE2, as.data.frame(E2), by=c("MGI.symbol"="genes")) # put mice genes on table
E2_H <- subset(E2_H, !duplicated(E2_H$HGNC.symbol))
index_E2 <- order(E2_H$FDR)
E2_H <- E2_H[index_E2,] # order by FDR
E2_H <- E2_H[, c(2, 4, 8)]
View(E2_H)
nrow(E2_H) # 11222
write.csv2(E2_H, "Results/EdgeR/Biomart/E2_H.csv", row.names = FALSE)

# tmod

martE2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_E2_H.csv")
index_E2 <- order(martE2$FDR)
martE2 <- martE2[index_E2,]
E2_c <- tmodCERNOtest(martE2$HGNC.symbol)
E2_LI<- list("Infected/Non_Infected"=E2_c)
write.csv2(E2_LI, file = "Results/EdgeR/Tmod/E2_LI.csv", row.names = FALSE)
TabelaE2 <- read.csv2('Results/EdgeR/Tmod/E2_LI.csv')
View(TabelaE2)

########## tmod pies ##########
# Create Pies
pieA<- tmodDecideTests(A2_H$HGNC.symbol, lfc=A2_H$logFC, pval = A2_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieB<- tmodDecideTests(B2_H$HGNC.symbol, lfc=B2_H$logFC, pval = B2_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieC<- tmodDecideTests(C2_H$HGNC.symbol, lfc=C2_H$logFC, pval = C2_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieD<- tmodDecideTests(D2_H$HGNC.symbol, lfc=D2_H$logFC, pval = D2_H$FDR, 
                       mset = "LI", pval.thr = 0.05, lfc.thr = 0)

pieE<- tmodDecideTests(E2_H$HGNC.symbol, lfc=E2_H$logFC, pval = E2_H$FDR, 
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

A2_LI<- list("Unstimulated"=A2_c)
B2_LI<- list("Not Infected"=B2_c)
C2_LI<- list("Infected"=C2_c)
D2_LI<- list("Stimulated"=D2_c)
E2_LI<- list("Infected/Non infected"=E2_c)

# Create the panels

# I_NS NI_NS vs I_T NI_NS
panel_listBC2 <- c(B2_LI, C2_LI)
is.list(panel_listBC2)
pie_listBC2 <- c(pieB, pieC)

tmodPanelPlot(panel_listBC2, pval.thr = 10^-4, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_listBC2, pie.style="pie")


# NI_T NI_NS vs I_T I_NS

panel_list_B_C <- c(B1_LI, C1_LI)
is.list(panel_list_B_C)
pie_list_BC <- c(pieB, pieC)


tmodPanelPlot(panel_list_B_C, pval.thr = 10^-4, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_BC, pie.style="pie")


# NI_T NI_NS vs I_T I_NS

panel_list_DE <- c(D2_LI, E2_LI)
is.list(panel_list_DE)
pie_list_DE <- c(pieD, pieE)

tmodPanelPlot(panel_list_DE, pval.thr = 10^-2, 
              pval.cutoff= 10^-30, filter.unknown = T, 
              text.cex =0.55, clust = "qval", pie =pie_list_DE, pie.style="pie")

# NOTHING TO PLOT
