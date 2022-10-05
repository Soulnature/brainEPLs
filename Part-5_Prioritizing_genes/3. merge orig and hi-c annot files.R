#==============================
# Merge orig+re
#==============================
rm(list = ls())
setwd("./")


#function of merging snp_to_gene1 and snp_to_gene2
merge_snp_to_gene = function(annot1, gene_bp1, snp_to_gene2, j){
  gene = c()
  bp = c()
  snp = c()
  for(i in 1:length(snp_to_gene2)){
    temp = unlist(strsplit(snp_to_gene2[i], "\t"))
    gene[i] = temp[1]
    bp[i] = temp[2]
    snp[i] = paste(temp[-2:-1],collapse="\t")
  }
  annot2 = data.frame(gene=gene,snp=snp,stringsAsFactors = FALSE)
  gene_bp2 = data.frame(gene=gene,bp=bp,stringsAsFactors = FALSE)
  
  
  ##Merge
  annot_comb = rbind(annot1, annot2)
  #If you have retained coding genes when you generating annotation files
  gene_bp = rbind(gene_bp1, gene_bp2)
  gene_bp = gene_bp[!duplicated(gene_bp$gene), ]
  
  library(dplyr)
  annot = annot_comb %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = "\t"))
  annot = merge(gene_bp, annot, by="gene",all=F)
  
  intermediate_file = paste(j, ".Intermediate.genes.annot", sep = '')
  write.table(annot, intermediate_file, row.names = F, col.names=F, sep="\t", quote=F)  
  
  #Remove duplicate SNPs on the same gene
  data = readLines(intermediate_file)
  data = as.list(data)
  temp = list()
  for(i in 1:length(data)){
    temp[[i]] = unlist(strsplit(data[[i]], "\t")) #\t
    temp[[i]] = as.vector(temp[[i]])
    temp[[i]] = unique(temp[[i]])
  }
  file.remove(intermediate_file)
  return(temp)
}


#-------------------
#The first file
#-------------------
snp_to_gene1 = readLines("annot/snp_to_gene.2kb.genes.annot") 
gene = c()
bp = c()
snp = c()
for(i in 1:length(snp_to_gene1)){
  temp = unlist(strsplit(snp_to_gene1[i], "\t"))
  gene[i] = temp[1]
  bp[i] = temp[2]
  snp[i] = paste(temp[-2:-1], collapse="\t") 
}
annot1 = data.frame(gene=gene, snp=snp, stringsAsFactors = FALSE)
gene_bp1 = data.frame(gene=gene, bp=bp, stringsAsFactors = FALSE)


#-------------------
#The second file
#-------------------
path = "annot/hi-c annot/"
file = list.files(path) 
brain_tissue = c()
for (i in 31:85){
  tissue_txt = paste('State_0', i, '.hi-c.annot', sep = '')
  brain_tissue = c(brain_tissue, tissue_txt)
}
file = subset(file, file %in% brain_tissue) 
dir = paste(path, file, sep = "")
n = length(dir)

write_temp = function(i){
  snp_to_gene2 = readLines(dir[i])
  temp = merge_snp_to_gene(annot1, gene_bp1, snp_to_gene2, i)
  
  file_name = gsub('hi-c', 'orig+hi-c.genes', file[i])
  file_name = paste('annot/orig_hi-c_annot/', file_name, sep = '')
  
  #写入orig+hi-c.genes.annot
  for(j in 1:length(temp)){
    cat(temp[[j]], file = file_name, sep = "\t", fill = FALSE, labels = NULL, append = TRUE)
    cat("\n", file = file_name, sep = "\t", fill = FALSE, labels = NULL, append = TRUE)
  }
  
  cat(i, '\n')
  return(0)
}


#----------
#并行计算
#----------
library(doParallel)
library(foreach)
num_cores = detectCores() #detectCores: 查看自己多少核
cl = makeCluster(4) 
registerDoParallel(cl)
clusterExport(cl, c("annot1", "gene_bp1")) # 导入全局环境中的数据(这步不能省)
t0 = Sys.time()
result = foreach(i = 1:n, .combine = 'rbind') %dopar% write_temp(i)  # %dopar并行，%do串行
Sys.time() - t0

stopCluster(cl)










