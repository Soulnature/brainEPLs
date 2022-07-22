#--------------------------
#提取基因坐标
#--------------------------
rm(list = ls())
setwd("F:/类脑/1_实操文件/yyc_postgwas")
library(data.table)
library(tidyverse)

# Gencode v38官网自己转成hg19的文件
f1 = fread("loc_file/gencode.v38lift37.annotation.gtf", header = F, data.table = F, check.names = F, sep = "\t")
f2 = subset(f1, f1$V3 == "gene")
f2 = f2[, -6:-8]
#!!不用X和Y染色体
chr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
        "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chr19","chr22","chr21")
f2 = subset(f2, f2$V1%in%chr)
f2$V1 = gsub('chr', '', f2$V1)

f22 = separate(data = f2, col = V9, into = c("gene_id", "gene_type","gene_name","gene_version","X1", "X2"),sep = ";")
f3 = f22[,c(-2,-3,-9:-11)]
f3$gene_id = gsub("gene_id ","",f3$gene_id); f3$gene_id = gsub("\"","",f3$gene_id)
f3$gene_type = gsub(" gene_type ","",f3$gene_type); f3$gene_type = gsub("\"","",f3$gene_type)
f3$gene_name = gsub(" gene_name ","",f3$gene_name); f3$gene_name = gsub("\"","",f3$gene_name)

names(f3) = c("chr","begin","end","gene_id","gene_type","gene_name")
f33 = subset(f3, f3$gene_type == "protein_coding")
f33 = f33[!duplicated(f33$gene_id), ]
#!!去掉_1
f33$gene_id = gsub('_1', '', f33$gene_id)

#变成MAGMA格式
gene_loc = f33[,c(4,1:3,6)]   

write.table(gene_loc, "loc_file/protein_coding_genes.gencode_v38lift37.chr1-22.loc", sep="\t", row.names=F, col.names=F, quote=F)
