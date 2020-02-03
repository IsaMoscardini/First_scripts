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

martA7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_A7_H.csv")
#View(martA1)
index_A7 <- order(martA7$FDR)
martA7 <- martA7[index_A7,]
A7_c <- tmodCERNOtest(martA7$HGNC.symbol, mset= msig[C2])
View(A7_c)
write.csv2(A7_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_A7_c.csv", row.names = FALSE)
TabelaA7 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_A7_c.csv")
View(TabelaA7)

#### Comparison B: NI_T X NI_NS
# tmod

martB7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_B7_H.csv")
index_B7 <- order(martB7$FDR)
martB7 <- martB7[index_B7,]
View(martB7)
B7_c <- tmodCERNOtest(martB7$HGNC.symbol, mset= msig[C2])
View(B7_c)
write.csv2(B7_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_B7_c.csv", row.names = FALSE)
TabelaB7 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_B7_c.csv")
View(TabelaB7)

#### Comparison C: I_T X I_NS
# tmod

martC7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_C7_H.csv")
index_C7 <- order(martC7$FDR)
martC7 <- martC7[index_C7,]
C7_c <- tmodCERNOtest(martC7$HGNC.symbol, mset= msig[C2])
View(C7_c)
write.csv2(C7_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_C7_c.csv", row.names = FALSE)
TabelaC7 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_C7_c.csv")
View(TabelaC7)

#### Comparison D: I_T X NI_T
# tmod

martD7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_D7_H.csv")
index_D7 <- order(martD7$FDR)
martD7 <- martD7[index_D7,]
D7_c <- tmodCERNOtest(martD7$HGNC.symbol, mset= msig[C2])
View(D7_c)
write.csv2(D7_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_D7_c.csv", row.names = FALSE)
TabelaD7 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_D7_c.csv")
View(TabelaD7)

#### Comparison E: I_T X NI_NS
# tmod

martE7 <- read.csv2("Results/EdgeR/Tables_remove_biomart_duplicates/Filter_duplic_by_FDR/by_FDR_SCORE_E7_H.csv")
index_E7 <- order(martE7$FDR)
martE7 <- martE7[index_E7,]
E7_c <- tmodCERNOtest(martE7$HGNC.symbol, mset= msig[C2])
View(E7_c)
write.csv2(E7_c, file = "Results/EdgeR/Tmod/Reactome/Reactome_E7_c.csv", row.names = FALSE)
TabelaE7 <- read.csv2("Results/EdgeR/Tmod/Reactome/Reactome_E7_c.csv")
View(TabelaE7)

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

