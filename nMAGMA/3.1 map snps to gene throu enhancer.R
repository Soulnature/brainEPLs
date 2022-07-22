#==============================
# Create Hi-C annotate files
#==============================
rm(list = ls())
setwd("F:/类脑/1_实操文件/yyc_postgwas")
library(data.table)
library(dplyr) 
library(tidyr)

# assign the SNPs to genes via SNP-to-RE annotations and RE-gene regulatory pairs
snp_to_gene = function(re_to_gene, gene_bp){
  re = c()
  snp = c()
  for(i in 1:length(snp_to_re)){
    temp = unlist(strsplit(snp_to_re[i], "\t"))
    re[i] = temp[1]
    snp[i] = paste(temp[-2:-1],collapse=" ")
  }
  snp_to_re2 = data.frame(re=re, snp=snp, stringsAsFactors = FALSE)
  
  b = intersect(re_to_gene[,1], snp_to_re2[,1]) #re取交集
  b = as.vector(b)
  re_to_gene = subset(re_to_gene, re_to_gene$V1%in%b)
  snp_to_re2 = subset(snp_to_re2, snp_to_re2[,1]%in%b)
  snp_to_gene = merge(re_to_gene, snp_to_re2, by.x="V1", by.y="re", all.x=T)
  snp_to_gene = snp_to_gene[,-1] #去掉re
  
  #合并同一个基因上的snp
  names(snp_to_gene) = c("gene", "snp")
  snp_to_gene = snp_to_gene %>% group_by(gene) %>% summarise(snp = paste(snp, collapse = " "))
  
  # add gene coordinates
  names(gene_bp) = c("gene", "chr", "begin", "end")
  gene_bp = tidyr::unite(gene_bp, "bp", chr, begin, sep = ":")
  gene_bp = tidyr::unite(gene_bp, "bp", bp, end, sep = ":")
  gene_bp = gene_bp[, c(1,2)]
  
  final_snp_to_gene = merge(snp_to_gene, gene_bp, by.x="gene", by.y="gene")
  final_snp_to_gene = final_snp_to_gene[, c(1,3,2)]
  write.table(final_snp_to_gene,"Intermediate.genes.annot", row.names = F, col.names=F, sep=" ", quote=F)
  
  #Remove duplicate SNPs on the same gene
  data = readLines("Intermediate.genes.annot")
  data = as.list(data)
  temp = list()
  for(i in 1:length(data)){
    temp[[i]] = unlist(strsplit(data[[i]], " "))
    temp[[i]] = as.vector(temp[[i]])
    temp[[i]] = unique(temp[[i]])
  }
  file.remove('Intermediate.genes.annot')
  return(temp)
}


snp_to_re = readLines("annot/snp_to_enhancer.genes.annot")  #要去掉annot文件开头的两行注释
gene_bp = read.table("loc_file/protein_coding_genes.gencode_v38lift37.chr1-22.loc", header = F,  sep = "\t",stringsAsFactors = F)
gene_bp = gene_bp[,1:4]   #去掉最后一列加了基因名字的

path = "RE_gene/"
file = list.files(path) 
#脑组织文件名
brain_tissue = c()
for (i in 31:85){
  tissue_txt = paste('RE_gene.State_0', i, '.txt', sep = '')
  brain_tissue = c(brain_tissue, tissue_txt)
}
file = subset(file, file %in% brain_tissue) #暂时只做脑组织
dir = paste(path, file, sep = "") #用paste命令构建路径变量dir  
n = length(dir) #文件夹下的文件个数

for (i in 1:n){
  re_to_gene = read.table(file = dir[i], header = F, sep = "\t",stringsAsFactors = F)   #同时有enhancer和promoter的
  temp = snp_to_gene(re_to_gene, gene_bp)  #返回Intermediate.genes.annot
  
  file_name = gsub('.txt', '', file[i])
  file_name = gsub('RE_gene.', '', file_name)
  file_name = paste('annot/hi-c annot/', file_name, '.hi-c.annot', sep = '')
  
  #写入snp_to_gene
  for(j in 1:length(temp)){
    cat(temp[[j]], file = file_name, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
    cat("\n", file = file_name, sep = "\t", fill = FALSE, labels = NULL,append = TRUE)
  }
  cat(i, '\n')
}

















