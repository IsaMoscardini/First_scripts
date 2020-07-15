rm(list = ls())

pkgs <- c('tidyverse', 'data.table','FactoMineR','ggplot2')
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

# Test if tables are congruent
identical(rownames(counts), rownames(cytok))
identical(rownames(counts), rownames(pheno))

# Preview Cytokines
boxplot(cytok, las = 2)
min(cytok)
cytok <- log10(cytok + .5)
boxplot(cytok, las = 2)

# Preview pheno
pheno = pheno %>% rename(Group = X.2) 
head(pheno)

# Preview counts
counts.melt = counts %>% as.data.frame %>% rownames_to_column() %>% 
  gather(gene, value, -rowname) %>% 
  merge(pheno, by.y = 'row.names', by.x = 'rowname') 

counts.melt %>% ggplot(aes(value, group = rowname, fill = Group)) + 
  geom_density(alpha = .6) + theme_linedraw() + theme(legend.position = 'top') +
  facet_grid(TP ~., scales = 'free') + labs(x= 'Expression Value')

#--- Prepare data for FactoMine/Shiny
data <- Reduce(cbind, list(pheno, cytok, counts))
(ncol(pheno) + ncol(cytok) + ncol(counts)) == ncol(data)

# Analysis
Factoshiny::Factoshiny(res = data)

cytokines_idx <- which(colnames(data) %in% colnames(cytok))
pheno_idx <- 1:2

res.PCA <- PCA(data,
               quali.sup=pheno_idx, quanti.sup=cytokines_idx, graph=T)
summary(res.PCA)
res.PCA$quali.sup
plot.PCA(res.PCA, choix='var', select='cos2  1.2',
         unselect=0, cex=0.5, cex.main=0.5, cex.axis=0.5,
         col.quanti.sup='#0000FF')
plotellipses(res.PCA, keepvar=1,invisible=c('ind.sup'),select='cos2  0',label =c('quali'))
dimdesc(res.PCA)
