#####
library(dplyr)
library(tidyr)
library(ggplot2)
#####

A1 <- read.csv2("Results/EdgeR/DEGs/DEGs_each_timepoint/TP_1_Non stimulated A.csv")
qplot(A1$PValue, geom="histogram",fill=I("blue"), col=I("darkblue"), alpha=I(.2),
      main = paste0("Histogram of P-value Frequency"), ylab = "Frequency", xlab = "P-value")
ggsave(plot, file=paste0("Histogram",".png"))

ff <- list.files(path="Results/EdgeR/DEGs/DEGs_each_timepoint/", full.names=TRUE)
tabelas <- lapply(ff, read.csv2) %>% setNames(basename(ff))
#View(tabelas[[1]])

for (i in names(tabelas)){
  tabela = tabelas[[i]]
  plot = qplot(tabela$PValue, geom="histogram", fill=I("blue"), col=I("darkblue"), alpha=I(.2),
        main = paste0("Histogram of P-value Frequency ",i), ylab = "Frequency", xlab = "P-value")
  ggsave(plot, file=paste0("Histogram_",i,".png"))
  }

