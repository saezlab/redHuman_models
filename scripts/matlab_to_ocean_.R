library(R.matlab)
library(stringr)
library(readr)

redHuman <- readMat("redHuman_recon3smin4aurelien.mat")

recon2_redhuman <- list()

pathway <- as.data.frame(matrix(NA,1732,1))
for(i in 1:1732)
{
  if(length(redHuman[[1]][[29]][[i]][[1]]) > 0)
  {
    pathway[i,1] <- redHuman[[1]][[29]][[i]]
  } else
  {
    pathway[i,1] <- NA
  }
  
}

names(pathway) <- "X1"

recon2_redhuman$pathway <- pathway

minmax <- as.data.frame(cbind(redHuman[[1]][[13]],redHuman[[1]][[14]]))
names(minmax) <- c("X1","X2")
minmax$X2 <- ifelse(minmax$X2 >= 0, minmax$X2, 0)
minmax$X1 <- ifelse(minmax$X1 <= 0, minmax$X1, 0)

minmax$direction <- ifelse(minmax$X1 < 0, #negative flux
                           ifelse(minmax$X2 > 0, # also positive flux
                                  ifelse(minmax$X2 + minmax$X1 > 0, #positive flux is higher 
                                         ifelse(abs(log10(minmax$X2) - log10(abs(minmax$X1))) > 1.3, # positive flux is 20 times higher
                                                1, #thus positive
                                                0 #thus bidirectional
                                                ), 
                                         ifelse(minmax$X2 + minmax$X1 == 0, #positive and negative flux are equal strenght
                                                0, #thus bidirecitonal
                                                ifelse(abs(log10(minmax$X2) - log10(abs(minmax$X1))) > 1.3, #negative flux is 20 time higher 
                                                       -1, #thus negative
                                                       0) #thus bidirectional
                                                )
                                         )
                                  , -1 #thus negative
                                  ),
                           ifelse(minmax$X2 > 0, 1, 0) #THIS SHOULD NOT HAPPEN
                           )

recon2_redhuman$minmax <- minmax

rules <- as.data.frame(matrix(NA,1732,1))
for(i in 1:1732)
{
  if(length(redHuman[[1]][[27]][[i]][[1]]) > 0)
  {
    rules[i,1] <- redHuman[[1]][[27]][[i]]
  } else
  {
    rules[i,1] <- NA
  }
  
}

names(rules) <- "X1"

recon2_redhuman$rules <- rules

rxns <- data.frame(unlist(redHuman[[1]][[2]]))
names(rxns) <- "X1"

recon2_redhuman$rxns <- rxns

stochio <- as.data.frame(redHuman[[1]][[1]])
names(stochio) <- paste("X",1:1732,sep = "")

recon2_redhuman$stochio <- stochio

metab <- data.frame(unlist(redHuman[[1]][[3]]))
names(metab) <- "X1"

recon2_redhuman$metab <- metab

library(readr)

gene_mapping <- as.data.frame(read_csv("gene_mapping.csv")) #entrez to symbol
metab_mapping <- as.data.frame(read_csv("metab_mapping.csv")) #bigg name to kegg coumpound

recon2_redhuman$gene_mapping <- gene_mapping
recon2_redhuman$metab_mapping <- metab_mapping

save(recon2_redhuman, file = "recon2_redhuman.RData")
