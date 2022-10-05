#--------------------------------------------
Extract enhancer, and lifeover hg38 -> hg19
#--------------------------------------------
rm(list = ls())
setwd("./")
#hg38 -> hg19
f1 = read.table("loc_file/enhancer_region.lift_to_hg19.bed", header = F, check.names = F, sep = "\t")
f1$V1 = gsub('chr', '', f1$V1)
f1 = f1[,c(4,1:3)]
write.table(f1, "loc_file/enhancer.hg19.loc", sep="\t", row.names=F, col.names=F, quote=F)

#MAGMA annotate
#gene extenrd 2kb
/home/u/destop/Software/magma_v1.08b_static/magma --annotate window=2 --snp-loc snp.hg19.loc \
--gene-loc protein_coding_genes.gencode_v38lift37.chr1-22.loc  --out snp_to_gene.2kb
#enhancer
/home/u/destop/Software/magma_v1.08b_static/magma --annotate --snp-loc snp.hg19.loc \
--gene-loc enhancer.hg19.loc --out snp_to_enhancer

