rm(list = ls())

pkgs <- c('tidyverse', 'data.table','FactoMineR')
sapply(pkgs, require, character.only = T)

# Data
counts <- fread('pneumococcus/data/RNA_filtered.csv',data.table = F)
cytok <- fread('pneumococcus/data/Cytok_filtered.csv',data.table = F)

# Preprocessing
counts = counts %>% column_to_rownames('V1') %>% as.matrix
cytok = cytok %>% column_to_rownames('V1')
pheno <- cytok[,c('X.2','TP')]
cytok <- cytok %>% select(-c(X.2,TP)) %>% as.matrix
colnames(cytok) <- gsub('\\.','_',gsub('\\.\\..*','',colnames(cytok)))

identical(rownames(counts), rownames(cytok))
identical(rownames(counts), rownames(pheno))

boxplot(cytok, las = 2)
min(cytok)
cytok <- log10(cytok + .5)
boxplot(cytok, las = 2)

# Analysis
