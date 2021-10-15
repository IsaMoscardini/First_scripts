#############------------ PRE PROCESS ------------############

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
#install.packages("VennDiagram")
library(VennDiagram)

intersezioni<-function(confronti,l,conc1,conc2,conc3=3,conc4=4,conc5=5){
  nomi<-c(conc1,conc2,conc3,conc4,conc5)
  nomi<-as.data.frame(nomi)
  n<-1
  n1<-1
  n2<-1
  n3<-1
  while(n<l+1){
    while(n1<l+1){
      if (n!=n1){
        comuni<-confronti[,n] %in% confronti[,n1]
        com<- confronti[,n][comuni]
        if (n<n1) write.table(com ,paste (nomi[n,1], "vs", nomi[n1,1], ".txt", sep = ""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
        while (n2<l+1){
          if (n2!=n & n2!=n1){
            co<-confronti[,n2]%in% com
            c<-confronti[,n2][co]
            salvataggio<-paste0(nomi[n,1], "vs", nomi[n1,1], "vs", nomi[n2,1], sep="")
            if (n<n1 & n1<n2) write.table(c ,paste ( salvataggio, ".txt", sep = ""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
            if (l>3){
              while (n3<l+1){
                if (n3!=n & n3!=n1 & n3!=n2){
                  ccc<-confronti[,n3] %in% c
                  cc<-confronti[,n3][ccc]
                  salvataggio<-paste0 (nomi[n,1], "vs", nomi[n1,1], "vs", nomi[n2,1], "vs", nomi[n3,1], sep="")
                  if(n<n1 & n1<n2 & n2<n3)write.table(cc ,paste ( salvataggio, ".txt", sep = ""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
                }
                n3<-n3+1
              }
              n3<-1
            }
          }
          n2<-n2+1
        }
      }
      n2<-1
      n1<-n1+1
    }
    n1<-1
    n<-n+1
    if (l==5){
      comunitutti<-confronti[,1] %in% cc
      comun<-confronti[,1][comunitutti]
      salvataggio<-paste0(conc1, "vs", conc2, "vs", conc3, "vs", conc4, "vs", conc5, sep = "")
      write.table(comun ,paste (salvataggio, ".txt", sep = ""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
    }
  }
}

#############------------ DEG Analysis ------------############

# Differential expression analysis - cutoff
ff <- list.files(path="Results/EdgeR/DEGs/DEGs_each_timepoint/", full.names=TRUE)
tabelas <- lapply(ff, read.csv2) %>% setNames(basename(ff))

# create a loop to run only once
ff <- list.files(path="Results/EdgeR/DEGs/DEGs_each_timepoint/", full.names=TRUE)
tabelas <- lapply(ff, read.csv2) %>% setNames(basename(ff))
for (i in names(tabelas)){
  tabela = tabelas[[i]]
  degs <-  tabela[which(abs(tabela$logFC) >= 0.322 & tabela$FDR < .01),]
  write.csv2(degs, paste0('DEGs_',i))
}

#tabelaA <- read.csv2("DEGs_TP_1_Infected C.csv")
#View(tabelaA) # apparently it is fine


#############------------ Venn Diagram ------------#############

ff <- list.files(path="Results/EdgeR/DEGs/DEGs_each_timepoint/", full.names=TRUE)
tabelas <- lapply(ff, read.csv2) %>% setNames(basename(ff))

for (i in names(tabelas)){
  tabela = tabelas[[i]]
  degs <-  tabela[which(tabela$FDR < .05),]
  write.csv2(degs, paste0('Venn',i))
}

##### Day 01 ####

con<-read.csv("Results/Table_Venn_Timepoint1.csv",  check.names = FALSE, stringsAsFactors = FALSE)
View(con)
#comando per richiamare la funzione
intersezioni(con,5,1,2,3,4,5)

n12 <- read.csv("Results/Intersections/TP 1/1vs2.txt")
View(n12)
n13 <- read.csv("Results/Intersections/TP 1/1vs3.txt")
View(n13)
n14 <- read.csv("Results/Intersections/TP 1/1vs4.txt")
View(n14)
n15 <- read.csv("Results/Intersections/TP 1/1vs5.txt")
View(n15)
#####
n23 <- read.csv("Results/Intersections/TP 1/2vs3.txt")
View(n23)
n24 <- read.csv("Results/Intersections/TP 1/2vs4.txt")
View(n24)
n25 <- read.csv("Results/Intersections/TP 1/2vs5.txt")
View(n25)
n34 <- read.csv("Results/Intersections/TP 1/3vs4.txt")
View(n34)
#####
n35 <- read.csv("Results/Intersections/TP 1/3vs5.txt")
View(n35)
n45 <- read.csv("Results/Intersections/TP 1/4vs5.txt")
View(n45)
n123 <- read.csv("Results/Intersections/TP 1/1vs2vs3.txt")
View(n123)
n124 <- read.csv("Results/Intersections/TP 1/1vs2vs4.txt")
View(n124)
#####
n125 <- read.csv("Results/Intersections/TP 1/1vs2vs5.txt")
View(n125)
n134 <- read.csv("Results/Intersections/TP 1/1vs3vs4.txt")
View(n134)
n135 <- read.csv("Results/Intersections/TP 1/1vs3vs5.txt")
View(n135)
n145 <- read.csv("Results/Intersections/TP 1/1vs4vs5.txt")
View(n145)
n234 <- read.csv("Results/Intersections/TP 1/2vs3vs4.txt")
View(n234)
n235 <- read.csv("Results/Intersections/TP 1/2vs3vs5.txt")
View(n235)
n245 <- read.csv("Results/Intersections/TP 1/2vs4vs5.txt")
View(n245)
#####
n345 <- read.csv("Results/Intersections/TP 1/3vs4vs5.txt")
View(n345)
n1234 <- read.csv("Results/Intersections/TP 1/1vs2vs3vs4.txt")
View(n1234)
n1235 <- read.csv("Results/Intersections/TP 1/1vs2vs3vs5.txt")
View(n1235)
n1245 <- read.csv("Results/Intersections/TP 1/1vs2vs4vs5.txt")
View(n1245)
n1345 <- read.csv("Results/Intersections/TP 1/1vs3vs4vs5.txt")
View(n1345)
#####
n2345 <- read.csv("Results/Intersections/TP 1/2vs3vs4vs5.txt")
View(n2345)
n12345 <- read.csv("Results/Intersections/TP 1/1vs2vs3vs4vs5.txt")
View(n12345)
#####

area1 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Non stimulated A.csv")
View(area1$genes)
area2 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Non infected B.csv")
View(area2$genes)
area3 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Infected C.csv")
View(area3$genes)
area4<- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Stimulated D.csv")
View(area4$genes)
area5 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Infected_Non infected E.csv")
View(area5$genes)

venn.plot <- draw.quintuple.venn(
  area1 = 7470,
  area2 = 7656,
  area3 = 1553,
  area4 = 1300,
  area5 = 8527,
  n12 = 5905,
  n13 = 611,
  n14 = 582,
  n15 = 6586,
  n23 = 996,
  n24 = 434,
  n25 = 6616,
  n34 = 208,
  n35 = 1142,
  n45 = 809,
  n123 = 409,
  n124 = 204,
  n125 = 5640,
  n134 = 90,
  n135 = 368,
  n145 = 511,
  n234 = 104,
  n235 = 899,
  n245 = 203,
  n345 = 168,
  n1234 = 75,
  n1235 = 344,
  n1245 = 170,
  n1345 = 69,
  n2345 = 85,
  n12345 = 59,
  category = c("A", "B", "C", "D", "E"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  margin = 0.05,
  ind = TRUE
);


##### Day 02

n12 <- read.csv("Results/Intersections/TP 1/1vs2.txt")
View(n12)
n13 <- read.csv("Results/Intersections/TP 1/1vs3.txt")
View(n13)
n14 <- read.csv("Results/Intersections/TP 1/1vs4.txt")
View(n14)
n15 <- read.csv("Results/Intersections/TP 1/1vs5.txt")
View(n15)
##### Day 02 ####
con<-read.csv("Results/Table_Venn_Timepoint2.csv",  check.names = FALSE, stringsAsFactors = FALSE)
View(con)
#comando per richiamare la funzione
intersezioni(con,5,1,2,3,4,5)

#####
n12 <- read.csv("Results/Intersections/TP 2/1vs2.txt")
View(n12)
n13 <- read.csv("Results/Intersections/TP 2/1vs3.txt")
View(n13)
n14 <- read.csv("Results/Intersections/TP 2/1vs4.txt")
View(n14)
n15 <- read.csv("Results/Intersections/TP 2/1vs5.txt")
View(n15)
n23 <- read.csv("Results/Intersections/TP 2/2vs3.txt")
View(n23)
n24 <- read.csv("Results/Intersections/TP 2/2vs4.txt")
View(n24)
n25 <- read.csv("Results/Intersections/TP 2/2vs5.txt")
View(n25)
n34 <- read.csv("Results/Intersections/TP 2/3vs4.txt")
View(n34)
#####
n35 <- read.csv("Results/Intersections/TP 2/3vs5.txt")
View(n35)
n45 <- read.csv("Results/Intersections/TP 2/4vs5.txt")
View(n45)
n123 <- read.csv("Results/Intersections/TP 2/1vs2vs3.txt")
View(n123)
n124 <- read.csv("Results/Intersections/TP 2/1vs2vs4.txt")
View(n124)
#####
n125 <- read.csv("Results/Intersections/TP 2/1vs2vs5.txt")
View(n125)
n134 <- read.csv("Results/Intersections/TP 2/1vs3vs4.txt")
View(n134)
n135 <- read.csv("Results/Intersections/TP 2/1vs3vs5.txt")
View(n135)
n145 <- read.csv("Results/Intersections/TP 2/1vs4vs5.txt")
View(n145)
n234 <- read.csv("Results/Intersections/TP 2/2vs3vs4.txt")
View(n234)
n235 <- read.csv("Results/Intersections/TP 2/2vs3vs5.txt")
View(n235)
n245 <- read.csv("Results/Intersections/TP 2/2vs4vs5.txt")
View(n245)
#####
n345 <- read.csv("Results/Intersections/TP 2/3vs4vs5.txt")
View(n345)
n1234 <- read.csv("Results/Intersections/TP 2/1vs2vs3vs4.txt")
View(n1234)
n1235 <- read.csv("Results/Intersections/TP 2/1vs2vs3vs5.txt")
View(n1235)
n1245 <- read.csv("Results/Intersections/TP 2/1vs2vs4vs5.txt")
View(n1245)
n1345 <- read.csv("Results/Intersections/TP 2/1vs3vs4vs5.txt")
View(n1345)
#####
n2345 <- read.csv("Results/Intersections/TP 2/2vs3vs4vs5.txt")
View(n2345)
n12345 <- read.csv("Results/Intersections/TP 2/1vs2vs3vs4vs5.txt")
View(n12345)
#####
area12 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Non stimulated A.csv")
View(area1$genes)
area2 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Non infected B.csv")
View(area2$genes)
area3 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Infected C.csv")
View(area3$genes)
area4<- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Stimulated D.csv")
View(area4$genes)
area5 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Infected_Non infected E.csv")
View(area5$genes)

venn.plot <- draw.quintuple.venn(
  area1 = 8096,
  area2 = 7844,
  area3 = 2633,
  area4 = 3094,
  area5 = 9176,
  n12 = 6117,
  n13 = 1366,
  n14 = 1699,
  n15 = 7120,
  n23 = 1636,
  n24 = 1240,
  n25 = 6723,
  n34 = 839,
  n35 = 1967,
  n45 = 2185,
  n123 = 987,
  n124 = 743,
  n125 = 5809,
  n134 = 414,
  n135 = 965,
  n145 = 1485,
  n234 = 422,
  n235 = 1430,
  n245 = 749,
  n345 = 678,
  n1234 = 297,
  n1235 = 848,
  n1245 = 616,
  n1345 = 308,
  n2345 = 319,
  n12345 = 216,
  category = c("A", "B", "C", "D", "E"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  margin = 0.05,
  ind = TRUE
);



##### Day 04 ####
con<-read.csv("Results/Table_Venn_Timepoint4.csv",  check.names = FALSE, stringsAsFactors = FALSE)
View(con)
#comando per richiamare la funzione
intersezioni(con,5,1,2,3,4,5)

#####
n12 <- read.csv("Results/Intersections/TP 4/1vs2.txt")
View(n12) # 
n13 <- read.csv("Results/Intersections/TP 4/1vs3.txt")
View(n13) # 
n14 <- read.csv("Results/Intersections/TP 4/1vs4.txt")
View(n14) # 
n15 <- read.csv("Results/Intersections/TP 4/1vs5.txt")
View(n15) # 
n23 <- read.csv("Results/Intersections/TP 4/2vs3.txt")
View(n23) # 
n24 <- read.csv("Results/Intersections/TP 4/2vs4.txt")
View(n24) # 
n25 <- read.csv("Results/Intersections/TP 4/2vs5.txt")
View(n25) # 
#####
n34 <- read.csv("Results/Intersections/TP 4/3vs4.txt")
View(n34) # 
n35 <- read.csv("Results/Intersections/TP 4/3vs5.txt")
View(n35) # 
n35 <- intersect(con$C, con$E) 
View(con$C)


n45 <- read.csv("Results/Intersections/TP 4/4vs5.txt")
View(n45) # 
n123 <- read.csv("Results/Intersections/TP 4/1vs2vs3.txt")
View(n123) # 
n124 <- read.csv("Results/Intersections/TP 4/1vs2vs4.txt")
View(n124) # 
#####
n125 <- read.csv("Results/Intersections/TP 4/1vs2vs5.txt")
View(n125) # 
n134 <- read.csv("Results/Intersections/TP 4/1vs3vs4.txt")
View(n134) # 
n135 <- read.csv("Results/Intersections/TP 4/1vs3vs5.txt")
View(n135) # 
n145 <- read.csv("Results/Intersections/TP 4/1vs4vs5.txt")
View(n145) # 
n234 <- read.csv("Results/Intersections/TP 4/2vs3vs4.txt")
View(n234) # 
n235 <- read.csv("Results/Intersections/TP 4/2vs3vs5.txt")
View(n235) # 
n245 <- read.csv("Results/Intersections/TP 4/2vs4vs5.txt")
View(n245) # 
#####
n345 <- read.csv("Results/Intersections/TP 4/3vs4vs5.txt")
View(n345) # 
n1234 <- read.csv("Results/Intersections/TP 4/1vs2vs3vs4.txt")
View(n1234) # 
n1235 <- read.csv("Results/Intersections/TP 4/1vs2vs3vs5.txt")
View(n1235) # 
n1245 <- read.csv("Results/Intersections/TP 4/1vs2vs4vs5.txt")
View(n1245) # 
n1345 <- read.csv("Results/Intersections/TP 4/1vs3vs4vs5.txt")
View(n1345) # 
#####
n2345 <- read.csv("Results/Intersections/TP 4/2vs3vs4vs5.txt")
View(n2345) # 
n12345 <- read.csv("Results/Intersections/TP 4/1vs2vs3vs4vs5.txt")
View(n12345) # 
#####
area1 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Non stimulated A.csv")
View(area1$genes) # 
area2 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Non infected B.csv")
View(area2$genes) # 
area3 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Infected C.csv")
View(area3$genes) # 
area4<- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Stimulated D.csv")
View(area4$genes) # 
area5 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Infected_Non infected E.csv")
View(area5$genes) # 
#####

venn.plot <- draw.quintuple.venn(
  area1 = 7494,
  area2 = 6936,
  area3 = 1802,
  area4 = 3075,
  area5 = 8689,
  n12 = 5441,
  n13 = 753,
  n14 = 1578,
  n15 = 6600,
  n23 = 973,
  n24 = 1100,
  n25 = 6041,
  n34 = 714,
  n35 = 1356,
  n45 = 2195,
  n123 = 487,
  n124 = 668,
  n125 = 5269,
  n134 = 321,
  n135 = 499,
  n145 = 1454,
  n234 = 348,
  n235 = 891,
  n245 = 727,
  n345 = 614,
  n1234 = 244,
  n1235 = 428,
  n1245 = 605,
  n1345 = 270,
  n2345 = 309,
  n12345 = 212,
  category = c("A", "B", "C", "D", "E"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  margin = 0.05,
  ind = TRUE
  );



##### Day 07 ####
con<-read.csv("Results/Table_Venn_Timepoint7.csv",  check.names = FALSE, stringsAsFactors = FALSE)
View(con)
#comando per richiamare la funzione
intersezioni(con,5,1,2,3,4,5)

#####
n12 <- read.csv("Results/Intersections/TP 7/1vs2.txt")
dim(n12) # 
n13 <- read.csv("Results/Intersections/TP 7/1vs3.txt")
dim(n13) # 
n14 <- read.csv("Results/Intersections/TP 7/1vs4.txt")
dim(n14) # 
n15 <- read.csv("Results/Intersections/TP 7/1vs5.txt")
dim(n15) # 
n23 <- read.csv("Results/Intersections/TP 7/2vs3.txt")
dim(n23) # 
n24 <- read.csv("Results/Intersections/TP 7/2vs4.txt")
dim(n24) # 
n25 <- read.csv("Results/Intersections/TP 7/2vs5.txt")
dim(n25) # 
#####
n34 <- read.csv("Results/Intersections/TP 7/3vs4.txt")
dim(n34) # 
n35 <- read.csv("Results/Intersections/TP 7/3vs5.txt")
dim(n35) # 
n45 <- read.csv("Results/Intersections/TP 7/4vs5.txt")
dim(n45) # 
n123 <- read.csv("Results/Intersections/TP 7/1vs2vs3.txt")
dim(n123) # 
n124 <- read.csv("Results/Intersections/TP 7/1vs2vs4.txt")
dim(n124) # 
#####
n125 <- read.csv("Results/Intersections/TP 7/1vs2vs5.txt")
dim(n125) # 
n134 <- read.csv("Results/Intersections/TP 7/1vs3vs4.txt")
dim(n134) # 
n135 <- read.csv("Results/Intersections/TP 7/1vs3vs5.txt")
dim(n135) # 
n145 <- read.csv("Results/Intersections/TP 7/1vs4vs5.txt")
dim(n145) # 
n234 <- read.csv("Results/Intersections/TP 7/2vs3vs4.txt")
dim(n234) # 
n235 <- read.csv("Results/Intersections/TP 7/2vs3vs5.txt")
dim(n235) # 
n245 <- read.csv("Results/Intersections/TP 7/2vs4vs5.txt")
dim(n245) # 
#####
n345 <- read.csv("Results/Intersections/TP 7/3vs4vs5.txt")
dim(n345) # 
n1234 <- read.csv("Results/Intersections/TP 7/1vs2vs3vs4.txt")
dim(n1234) # 
n1235 <- read.csv("Results/Intersections/TP 7/1vs2vs3vs5.txt")
dim(n1235) # 
n1245 <- read.csv("Results/Intersections/TP 7/1vs2vs4vs5.txt")
dim(n1245) # 
n1345 <- read.csv("Results/Intersections/TP 7/1vs3vs4vs5.txt")
dim(n1345) # 
#####
n2345 <- read.csv("Results/Intersections/TP 7/2vs3vs4vs5.txt")
dim(n2345) # 
n12345 <- read.csv("Results/Intersections/TP 7/1vs2vs3vs4vs5.txt")
dim(n12345) # 
#####
area1 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Non stimulated A.csv")
View(area1$genes) # 
area2 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Non infected B.csv")
View(area2$genes) # 
area3 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Infected C.csv")
View(area3$genes) # 
area4<- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Stimulated D.csv")
View(area4$genes) # 
area5 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Infected_Non infected E.csv")
View(area5$genes) # 
#####

venn.plot <- draw.quintuple.venn(
  area1 = 6766,
  area2 = 6280,
  area3 = 933,
  area4 = 1753,
  area5 = 7833,
  n12 = 5000,
  n13 = 298,
  n14 = 863,
  n15 = 6051,
  n23 = 597,
  n24 = 480,
  n25 = 5539,
  n34 = 218,
  n35 = 708,
  n45 = 1170,
  n123 = 194,
  n124 = 245,
  n125 = 4880,
  n134 = 84,
  n135 = 164,
  n145 = 769,
  n234 = 120,
  n235 = 531,
  n245 = 257,
  n345 = 179,
  n1234 = 70,
  n1235 = 153,
  n1245 = 210,
  n1345 = 57,
  n2345 = 96,
  n12345 = 50,
  category = c("A", "B", "C", "D", "E"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  margin = 0.05,
  ind = TRUE
);


##############################################################################################################
#############------------ Venn Diagram BC 2x2 ------------#############
##### Day 01 #####
venn.plot <- draw.pairwise.venn(
  area1 = 7656,
  area2 = 1553,
  cross.area = 996,
  category = c("B", "C"),
  fill = c("orchid3", "darkorange1"),
  cat.col = c("orchid3", "darkorange1"),
  margin = 0.05,
  ind = TRUE
);

##### Day 02 #####
venn.plot <- draw.pairwise.venn(
  area1 = 7844,
  area2 = 2633,
  cross.area = 1636,
  category = c("B", "C"),
  fill = c("orchid3", "darkorange1"),
  cat.col = c("orchid3", "darkorange1"),
  margin = 0.05,
  ind = TRUE
);
##### Day 04 #####
venn.plot <- draw.pairwise.venn(
  area1 = 6936,
  area2 = 1802,
  cross.area = 973,
  category = c("B", "C"),
  fill = c("orchid3", "darkorange1"),
  cat.col = c("orchid3", "darkorange1"),
  margin = 0.05,
  ind = TRUE
);

##### Day 07 #####
venn.plot <- draw.pairwise.venn(
  area1 = 6280,
  area2 = 933,
  cross.area = 597,
  category = c("B", "C"),
  fill = c("orchid3", "darkorange1"),
  cat.col = c("orchid3", "darkorange1"),
  margin = 0.05,
  ind = TRUE
);

###############################################################################################################
###############################################################################################################
#############------------ Find Interesting genes - with DEGs in venn ------------############
# Which genes are down regulated in comparisons B and up regulated in C? - different pattern after stimulation
##### Day1 #####
B1 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Non infected B.csv") # stimulation without infection
C1 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_1_Infected C.csv") # stimulation with infection

Down_B1 <- B1[B1$logFC < 0,]
Up_B1 <- B1[B1$logFC > 0,]
Down_C1 <- C1[C1$logFC < 0,]
Up_C1 <- C1[C1$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec1 <- intersect(Down_B1$genes, Up_C1$genes)
genes_up_infec1 <-  as.data.frame(genes_up_infec1)
View(genes_up_infec1)
dim(genes_up_infec1) # 49 genes
# good_genes_venn1 <- c("Cyp27b1", "Kcnd1", "Gm19510", "Qrich2", "Zfp787", "", "E130012A19Rik", "Tsen15", "Ints5", "Lsmem2", "Cxcr1", "Rrs1", "1810026J23Rik",
#                        "Fem1a", "Myo19", "2310039H08Rik", "Eid3", "Cyp1a1", "Zfp783", "Utf1", "Slc29a2", "Qtrt1", "Dhx37", "Acy1", "Nsun5", "Mrps10", "H2afx",
#                        "Amigo2", "Fhl3", "Slc25a22", "Upp1", "Slc25a39", "Noc4l", "Znhit6", "Gar1", "Ssbp4", "G0s2", "Alg8", "Larp7", "Elovl1", "Cd2bp2", "Mthfd1",
#                        "Tomm40", "Naa 38", "Nob1", "Csf3", "Sh2d5", "Trib3", "Il23a", "Wdr55")
good_genes_venn1 <- genes_up_infec1
View(good_genes_venn1)

# genes that are UP without infection and DOWN with infection 
genes_down_infec1 <- intersect(Up_B1$genes, Down_C1$genes)
genes_down_infec1 <- as.data.frame(genes_down_infec1)
View(genes_down_infec1)
dim(genes_down_infec1) # 15 genes

##### Day 2 #####

B2 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Non infected B.csv") # stimulation without infection
C2 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_2_Infected C.csv") # stimulation with infection

Down_B2 <- B2[B2$logFC < 0,]
Up_B2 <- B2[B2$logFC > 0,]
Down_C2 <- C2[C2$logFC < 0,]
Up_C2 <- C2[C2$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec2 <- intersect(Down_B2$genes, Up_C2$genes)
genes_up_infec2 <-  as.data.frame(genes_up_infec2)
View(genes_up_infec2)
dim(genes_up_infec2) # 70 genes
good_genes_venn2 <- genes_up_infec2  

# genes that are UP without infection and DOWN with infection 
genes_down_infec2 <- intersect(Up_B2$genes, Down_C2$genes)
genes_down_infec2 <- as.data.frame(genes_down_infec2)
View(genes_down_infec2)
dim(genes_down_infec2) # 71 genes

#####
# genes with same direction
#genes_same_direction2 <- intersect(Up_B2$genes, Up_C2$genes)
#genes_same_direction2 <-  as.data.frame(genes_same_direction2)
#View(genes_same_direction2)
#dim(genes_same_direction2) # 657 genes

#genes_same_direction_down2 <- intersect(Down_B2$genes, Down_C2$genes)
#genes_same_direction_down2  <-  as.data.frame(genes_same_direction_down2)
#View(genes_same_direction_down2)
#dim(genes_same_direction_down2)  # 839 genes


##### Day 4 ######
B4 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Non infected B.csv") # stimulation without infection
C4 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_4_Infected C.csv") # stimulation with infection

Down_B4 <- B4[B4$logFC < 0,]
Up_B4 <- B4[B4$logFC > 0,]
Down_C4 <- C4[C4$logFC < 0,]
Up_C4 <- C4[C4$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec4 <- intersect(Down_B4$genes, Up_C4$genes)
genes_up_infec4 <-  as.data.frame(genes_up_infec4)
View(genes_up_infec4)
dim(genes_up_infec4) # 73 genes
good_genes_venn4 <- genes_up_infec4

# genes that are UP without infection and DOWN with infection 
genes_down_infec4 <- intersect(Up_B4$genes, Down_C4$genes)
genes_down_infec4 <- as.data.frame(genes_down_infec4)
View(genes_down_infec4)
dim(genes_down_infec4) # 20 genes
#####
# genes with same direction
#genes_same_direction4 <- intersect(Up_B4$genes, Up_C4$genes)
#genes_same_direction4 <-  as.data.frame(genes_same_direction4)
#View(genes_same_direction4)
#dim(genes_same_direction4) # 470 genes

#genes_same_direction_down4 <- intersect(Down_B4$genes, Down_C4$genes)
#genes_same_direction_down4  <-  as.data.frame(genes_same_direction_down4)
#View(genes_same_direction_down4)
#dim(genes_same_direction_down4)  # 411 genes



##### Day 7 #####

B7 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Non infected B.csv") # stimulation without infection
C7 <- read.csv2("Results/EdgeR/DEGs/DEGs_Venn/VennTP_7_Infected C.csv") # stimulation with infection

Down_B7 <- B7[B7$logFC < 0,]
Up_B7 <- B7[B7$logFC > 0,]
Down_C7 <- C7[C7$logFC < 0,]
Up_C7 <- C7[C7$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec7 <- intersect(Down_B7$genes, Up_C7$genes)
genes_up_infec7 <-  as.data.frame(genes_up_infec7)
View(genes_up_infec7)
dim(genes_up_infec7) # 24 genes: 
good_genes_venn7 <- genes_up_infec7

# genes that are UP without infection and DOWN with infection 
genes_down_infec7 <- intersect(Up_B7$genes, Down_C7$genes)
genes_down_infec7 <- as.data.frame(genes_down_infec7)
View(genes_down_infec7)
dim(genes_down_infec7) # 3 genes

#####
# genes with same direction
# genes_same_direction7 <- intersect(Up_B7$genes, Up_C7$genes)
# genes_same_direction7  <-  as.data.frame(genes_same_direction7)
# View(genes_same_direction7)
# dim(genes_same_direction7) # 237 genes
# 
# genes_same_direction_down7 <- intersect(Down_B7$genes, Down_C7$genes)
# genes_same_direction_down7  <-  as.data.frame(genes_same_direction_down7)
# View(genes_same_direction_down77)
# dim(genes_same_direction_down7)  # 334 genes

#####

inter_venn12 <- intersect(good_genes_venn1$genes_up_infec1, good_genes_venn2$genes_up_infec2) # 
intervenn124 <- intersect(inter_venn12, good_genes_venn4$genes_up_infec4) # "Gm19510" E130012A19Rik" "Cxcr1" "Rrs1" "Myo19" "2310039H08Rik" 
                                                                          # "Cyp1a1" "Amigo2" "Upp1" "Slc25a39" "Noc4l" "Znhit6" "Gar1" 

inter1247 <- intersect(intervenn124, good_genes_venn7$genes_up_infec7) # "Cxcr1"  "Rrs1"   "Myo19"  "Amigo2" "Upp1"   "Noc4l" "Znhit6"

inter24 <- intersect(good_genes_venn2$genes_up_infec2, good_genes_venn4$genes_up_infec4)

##############################################################################################################
##############################################################################################################
#############------------ Find Interesting genes ------------############

# Which genes are down regulated in comparisons B and up regulated in C? - different pattern after stimulation
##### Day 1 #####

B1 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_1_Non infected B.csv") # stimulation without infection
C1 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_1_Infected C.csv") # stimulation with infection

Down_B1 <- B1[B1$logFC < 0,]
Up_B1 <- B1[B1$logFC > 0,]
Down_C1 <- C1[C1$logFC < 0,]
Up_C1 <- C1[C1$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec1 <- intersect(Down_B1$genes, Up_C1$genes)
genes_up_infec1 <-  as.data.frame(genes_up_infec1)
View(genes_up_infec1)
dim(genes_up_infec1) # 20 genes: Gm19510, Tsen15, Lsmem2, Cxcr1, Rrs1, Znhit6, Gar1, Ssbp4, G0s2, Myo19, Eid3, Cyp1a1, Utf1, Qtrt1, Mrps10, Amigo2, Slc25a22, Upp1, Slc25a39, Noc4l
good_genes1 <- c("Gm19510", "Tsen15", "Lsmem2", "Cxcr1", "Rrs1", "Znhit6", "Gar1", "Ssbp4", "G0s2", "Myo19", "Eid3", "Cyp1a1", "Utf1", "Qtrt1", "Mrps10", "Amigo2", "Slc25a22", "Upp1", "Slc25a39", "Noc4l")

# genes that are UP without infection and DOWN with infection 
genes_down_infec1 <- intersect(Up_B1$genes, Down_C1$genes)
genes_down_infec1 <- as.data.frame(genes_down_infec1)
View(genes_down_infec1)
dim(genes_down_infec1) # 0 genes

# genes with same direction
genes_same_direction1 <- intersect(Up_B1$genes, Up_C1$genes)
genes_same_direction1 <-  as.data.frame(genes_same_direction1)
View(genes_same_direction1)
dim(genes_same_direction1) # 253 genes

genes_same_direction_down1 <- intersect(Down_B1$genes, Down_C1$genes)
genes_same_direction_down1  <-  as.data.frame(genes_same_direction_down1)
View(genes_same_direction_down1)
dim(genes_same_direction_down1)  # 309 genes

##### Day 2 #####

B2 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_2_Non infected B.csv") # stimulation without infection
C2 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_1_Infected C.csv") # stimulation with infection

Down_B2 <- B2[B2$logFC < 0,]
Up_B2 <- B2[B2$logFC > 0,]
Down_C2 <- C2[C2$logFC < 0,]
Up_C2 <- C2[C2$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec2 <- intersect(Down_B2$genes, Up_C2$genes)
genes_up_infec2 <-  as.data.frame(genes_up_infec2)
View(genes_up_infec2)
dim(genes_up_infec2) # 20 genes: Cyp1a1, LOC105245439, 1110004E09Rik, 2310039H08Rik, Eid3, Rrs1, Usp16, Emc6, Amigo2, Kti12, Dimt1, Znhit6, Ndufc1, Cdk5r1, Ankrd28, Polr3d, Utp14a, Upp1, Snhg5, G0s2
good_genes2 <- c("Cyp1a1", "LOC105245439", "1110004E09Rik", "2310039H08Rik", "Eid3", "Rrs1", "Usp16", "Emc6", "Amigo2", "Kti12", "Dimt1", "Znhit6", "Ndufc1", "Cdk5r1", "Ankrd28", "Polr3d", "Utp14a", "Upp1", "Snhg5", "G0s2")  

# genes that are UP without infection and DOWN with infection 
genes_down_infec2 <- intersect(Up_B2$genes, Down_C2$genes)
genes_down_infec2 <- as.data.frame(genes_down_infec2)
View(genes_down_infec2)
dim(genes_down_infec2) # 1 gene

# genes with same direction
genes_same_direction2 <- intersect(Up_B2$genes, Up_C2$genes)
genes_same_direction2 <-  as.data.frame(genes_same_direction2)
View(genes_same_direction2)
dim(genes_same_direction2) # 254 genes

genes_same_direction_down2 <- intersect(Down_B2$genes, Down_C2$genes)
genes_same_direction_down2  <-  as.data.frame(genes_same_direction_down2)
View(genes_same_direction_down2)
dim(genes_same_direction_down2)  # 299 genes

##### Day 4 #####

B4 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_4_Non infected B.csv") # stimulation without infection
C4 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_4_Infected C.csv") # stimulation with infection

Down_B4 <- B4[B4$logFC < 0,]
Up_B4 <- B4[B4$logFC > 0,]
Down_C4 <- C4[C4$logFC < 0,]
Up_C4 <- C4[C4$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec4 <- intersect(Down_B4$genes, Up_C4$genes)
genes_up_infec4 <-  as.data.frame(genes_up_infec4)
View(genes_up_infec4)
dim(genes_up_infec4) # 22 genes
good_genes4 <- genes_up_infec4

# genes that are UP without infection and DOWN with infection 
genes_down_infec4 <- intersect(Up_B4$genes, Down_C4$genes)
genes_down_infec4 <- as.data.frame(genes_down_infec4)
View(genes_down_infec4)
dim(genes_down_infec4) # 3 genes

# genes with same direction
genes_same_direction4 <- intersect(Up_B4$genes, Up_C4$genes)
genes_same_direction4 <-  as.data.frame(genes_same_direction4)
View(genes_same_direction4)
dim(genes_same_direction4) # 302 genes

genes_same_direction_down4 <- intersect(Down_B4$genes, Down_C4$genes)
genes_same_direction_down4  <-  as.data.frame(genes_same_direction_down4)
View(genes_same_direction_down4)
dim(genes_same_direction_down4)  # 324 genes


##### Day 7 #####

B7 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_7_Non infected B.csv") # stimulation without infection
C7 <- read.csv2("Results/EdgeR/DEGs/DEGs_selected/DEGs_TP_7_Infected C.csv") # stimulation with infection

Down_B7 <- B7[B7$logFC < 0,]
Up_B7 <- B7[B7$logFC > 0,]
Down_C7 <- C7[C7$logFC < 0,]
Up_C7 <- C7[C7$logFC > 0,]

# genes that are DOWN without infection and UP with infection 
genes_up_infec7 <- intersect(Down_B7$genes, Up_C7$genes)
genes_up_infec7 <-  as.data.frame(genes_up_infec7)
View(genes_up_infec7)
dim(genes_up_infec7) # 5 genes: Ifng, Rrs1, Cxcr1, Amigo2, Upp1
good_genes7 <- genes_up_infec7

# genes that are UP without infection and DOWN with infection 
genes_down_infec7 <- intersect(Up_B7$genes, Down_C7$genes)
genes_down_infec7 <- as.data.frame(genes_down_infec7)
View(genes_down_infec7)
dim(genes_down_infec7) # 1 genes

# genes with same direction
genes_same_direction7 <- intersect(Up_B7$genes, Up_C7$genes)
genes_same_direction7  <-  as.data.frame(genes_same_direction7)
View(genes_same_direction7)
dim(genes_same_direction7) # 139 genes

genes_same_direction_down7 <- intersect(Down_B7$genes, Down_C7$genes)
genes_same_direction_down7  <-  as.data.frame(genes_same_direction_down7)
View(genes_same_direction_down7)
dim(genes_same_direction_down7)  # 212 genes

##### GOOD GENES ANALYSIS ##### 

inter_12 <- intersect(good_genes1, good_genes2) # "Rrs1" "Znhit6" "G0s2" "Eid3" "Cyp1a1" "Amigo2" "Upp1"
inter124 <- intersect(inter_12, good_genes4) # "Rrs1" "Cyp1a1" "Amigo2" "Upp1"
inter1247 <- intersect(inter124, good_genes7) # "Rrs1" "Amigo2" "Upp1"  
