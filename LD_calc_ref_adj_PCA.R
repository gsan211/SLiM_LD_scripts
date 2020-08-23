#script for calculating mean LD between deleterious and neutral mutations from SLiM script, includes steps for filtering out migrants via PCA for admixture simulations

#!/usr/bin/env Rscript

library(matrixStats)
library(data.table)
library(magrittr)
library(windowscanr)
library(dplyr)
library(gdsfmt)
library(SNPRelate)


#read file with paths of all simulation vcf outputs
ll = fread("/plas1/george.sandler/syn_epi/equal_pops/admixture100k/list_ad200k.list", header = F)

analyze  <- function(line_filename){
  line = nth(ll$V1, line_filename)
  
  
  vcf.fn <- line
  # Reformat
  snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")
  #snpgdsSummary("test.gds")
  genofile <- snpgdsOpen("test.gds", readonly = F)
  pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F)
  pc.percent <- pca$varprop*100
  head(round(pc.percent, 2))
  tab <- data.frame(sample.id = pca$sample.id,
                    EV1 = pca$eigenvect[,1], # the first eigenvector
                    EV2 = pca$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE)
                    
  #take output of SNPrelate PCA analysis and filter out individuals more then 1 SD away from mean PC1 or PC2, only neccessary to filter recent migrants in admixure sims              
  inliers = tab[tab$EV1 < sd(tab$EV1) & tab$EV1 >-sd(tab$EV1) & tab$EV2 < sd(tab$EV2) & tab$EV2 > -sd(tab$EV2),]
  kep = c(inliers$sample.id)
  closefn.gds(genofile)  

  multi2 = fread(line)
  df <- multi2 
  df <- lapply(df, gsub, pattern = "1|0", replacement = "1", fixed = TRUE)
  df <- lapply(df, gsub, pattern = "0|1", replacement = "1", fixed = TRUE)
  df <- lapply(df, gsub, pattern = "1|1", replacement = "2", fixed = TRUE)
  df <- lapply(df, gsub, pattern = "0|0", replacement = "0", fixed = TRUE)
  df2 <- data.frame(df)
  
  #separate selected/neutral mutations based on MT tag in VCF, MT=2 is deleterious
  dfs = df2[grep("MT=2", df2$INFO),]  
  dfs = dfs[,c(-1:-9)]
  dfs = subset(dfs, select=kep) 
  cols = c(1:(ncol(dfs)))
  dfs[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
  dfs = dfs[(rowSums(dfs)) < 200,] #remove fixed mutation from LD calculation
  #dfs = dfs[(rowSums(dfs) < 6),]
  
  #to replicate how LD is calculated in real world pops, derived mutations with a frequency of over 50% are set as the reference allele in the sims in this block of code
  dfs$sum = rowSums(dfs) 
  norms = dfs[dfs$sum < 100,]
  chans = dfs[dfs$sum > 99,]
  chans <- lapply(chans, gsub, pattern = "0", replacement = "zero", fixed = TRUE)
  chans <- lapply(chans, gsub, pattern = "2", replacement = "two", fixed = TRUE)
  chans <- lapply(chans, gsub, pattern = "zero", replacement = "2", fixed = TRUE)
  chans <- lapply(chans, gsub, pattern = "two", replacement = "0", fixed = TRUE)
  chans <- data.frame(chans)
  cols2s = c(1:(ncol(chans)))
  chans[,cols2s] %<>% lapply(function(x) as.numeric(as.character(x)))
  df3s = rbind(norms,chans)
  df3s = df3s[df3s$sum < 6,]
  df3s$sum = NULL
  
  #convert mutation counts to matrix then caluclate average LD by calculating covariane in mutation burden divided by the number of pairs of mutations (n choose 2)
  mat_s <- data.matrix(df3s) 
  netld_s=(var(colSums(mat_s, na.rm = T)) - sum(rowVars(mat_s, na.rm = T)))/choose(nrow(mat_s), 2)
  burden_s=mean(colSums(mat_s), na.rm = T)
  
  #same as above but for neutral mutations MT=1
  dfn = df2[grep("MT=1", df2$INFO),]
  dfn = dfn[,c(-1:-9)]
  dfn = subset(dfn, select=kep) 
  cols = c(1:(ncol(dfn)))
  dfn[,cols] %<>% lapply(function(x) as.numeric(as.character(x)))
  dfn = dfn[(rowSums(dfn)) < 200,]
  #dfn = dfn[(rowSums(dfn) < 6),]
  
  dfn$sum = rowSums(dfn)
  normn = dfn[dfn$sum < 100,]
  chann = dfn[dfn$sum > 99,]
  chann <- lapply(chann, gsub, pattern = "0", replacement = "zero", fixed = TRUE)
  chann <- lapply(chann, gsub, pattern = "2", replacement = "two", fixed = TRUE)
  chann <- lapply(chann, gsub, pattern = "zero", replacement = "2", fixed = TRUE)
  chann <- lapply(chann, gsub, pattern = "two", replacement = "0", fixed = TRUE)
  chann <- data.frame(chann)
  cols2n = c(1:(ncol(chann)))
  chann[,cols2n] %<>% lapply(function(x) as.numeric(as.character(x)))
  df3n = rbind(normn,chann)
  df3n = df3n[df3n$sum < 6,]
  df3n$sum = NULL
  
  mat_n <- data.matrix(df3n)  
  netld_n=(var(colSums(mat_n, na.rm = T)) - sum(rowVars(mat_n, na.rm = T)))/choose(nrow(mat_n), 2)
  burden_n=mean(colSums(mat_n), na.rm = T)


  #summarize results in a string including average LD between selected/neutral mutations, average burden of each mutation type and tag for type of simualtion
  line2=(noquote(paste(netld_s, netld_n, burden_s, burden_n, "200k_ad")))
                             
  #output file for writing out results
  write(line2,file="/plas1/george.sandler/syn_epi/equal_pops/results_100k_netLD_admix_ref_adj_PCA_maflow.txt",append=TRUE)
}

#loop above code over each simulation replicate defined in file on top, can use multiple cores to speed up
library(foreach)
library(doParallel)
registerDoParallel(cores=1)
foreach(i=1:nrow(ll), .errorhandling="pass") %dopar% analyze(i)
