rm(list = ls())
options(stringsAsFactors = F)
setwd("~/PROJECTS/Streptococcus pneumoniae")
library(biomaRt)
library(tmod)
library(dplyr)

msig=tmodImportMSigDB(file = "Reactome_2016", format = "gmt")
C2 <- msig$MODULES$ID

#### Comparison A: I_NS X NI_NS
# tmod

martA2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_A2_H.csv")
#View(martA1)
index_A2 <- order(martA2$FDR)
martA2 <- martA2[index_A2,]
A2_c <- tmodCERNOtest(martA2$HGNC.symbol, mset= msig[C2])
View(A2_c)
write.csv2(A2_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_A2_c.csv", row.names = FALSE)
TabelaA2 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_A2_c.csv")
View(TabelaA2)

#### Comparison B: NI_T X NI_NS
# tmod

martB2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_B2_H.csv")
index_B2 <- order(martB2$FDR)
martB2 <- martB2[index_B2,]
B2_c <- tmodCERNOtest(martB2$HGNC.symbol, mset= msig[C2])
View(B2_c)
write.csv2(B2_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_B2_c.csv", row.names = FALSE)
TabelaB2 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_B2_c.csv")
View(TabelaB2)

#### Comparison C: I_T X I_NS
# tmod

martC2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_C2_H.csv")
index_C2 <- order(martC2$FDR)
martC2 <- martC2[index_C2,]
C2_c <- tmodCERNOtest(martC2$HGNC.symbol, mset= msig[C2])
View(C2_c)
write.csv2(C2_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_C2_c.csv", row.names = FALSE)
TabelaC2 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_C2_c.csv")
View(TabelaC2)

#### Comparison D: I_T X NI_T
# tmod

martD2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_D2_H.csv")
index_D2 <- order(martD2$FDR)
martD2 <- martD2[index_D2,]
D2_c <- tmodCERNOtest(martD2$HGNC.symbol, mset= msig[C2])
View(D2_c)
write.csv2(D2_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_D2_c.csv", row.names = FALSE)
TabelaD2 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_D2_c.csv")
View(TabelaD2)

#### Comparison E: I_T X NI_NS
# tmod

martE2 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_E2_H.csv")
index_E2 <- order(martE2$FDR)
martE2 <- martE2[index_E2,]
E2_c <- tmodCERNOtest(martE2$HGNC.symbol, mset= msig[C2])
View(E2_c)
write.csv2(E2_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_E2_c.csv", row.names = FALSE)
TabelaE2 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_E2_c.csv")
View(TabelaE2)

########## tmod pies ##########
# Create Pies
pieA1<- tmodDecideTests(martA1$HGNC.symbol, lfc=martA1$logFC, pval = martA1$FDR, 
                        mset = "all", pval.thr = 0.05, lfc.thr = 0)
View(pieA1)

pieB1<- tmodDecideTests(martB1$HGNC.symbol, lfc=martB1$logFC, pval = martB1$FDR, 
                        mset = "Reactome_2016", pval.thr = 0.05, lfc.thr = 0)

pieC1<- tmodDecideTests(martC1$HGNC.symbol, lfc=martC1$logFC, pval = martC1$FDR, 
                        mset = "all", pval.thr = 0.05, lfc.thr = 0)

pieD1<- tmodDecideTests(martD1$HGNC.symbol, lfc=martD1$logFC, pval = martD1$FDR, 
                        mset = "all", pval.thr = 0.05, lfc.thr = 0)

pieE1<- tmodDecideTests(martE1$HGNC.symbol, lfc=martE1$logFC, pval = martE1$FDR, 
                        mset = "all", pval.thr = 0.05, lfc.thr = 0)

# transform in data frame
pieA1<- as.data.frame(pieA1)
pieB1<- as.data.frame(pieB1)
pieC1<- as.data.frame(pieC1)
pieD1<- as.data.frame(pieD1)
pieE1<- as.data.frame(pieE1)
#View(pieE)

# Remove the X if there is any
colnames(pieA1) <- gsub("X.*\\.", "", colnames(pieA1))
colnames(pieB1) <- gsub("X.*\\.", "", colnames(pieB1))
colnames(pieC1) <- gsub("X.*\\.", "", colnames(pieC1))
colnames(pieD1) <- gsub("X.*\\.", "", colnames(pieD1))
colnames(pieE1) <- gsub("X.*\\.", "", colnames(pieE1))

pieC<- list("Infected"=pieC1)
pieD<- list("Stimulated"= pieD1)
pieA<- list("Unstimulated"= pieA1)
pieB<- list("Not Infected"=pieB1)
pieE<- list("Infected/Non infected"=pieE1)

A1_LI<- list("Unstimulated"=A1_c)
B1_LI<- list("Not Infected"=B1_c)
C1_LI<- list("Infected"=C1_c)
D1_LI<- list("Stimulated"=D1_c)
E1_LI<- list("Infected/Non infected"=E1_c)

# Create the panels

# I_NS NI_NS vs I_T NI_NS
panel_list_A_E <- c(A1_LI, E1_LI)
is.list(panel_list_A_E)
pie_list_A_E <- c(pieA, pieE)

tmodPanelPlot(panel_list_A_E, pval.thr = 10^-4, 
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

