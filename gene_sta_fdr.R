## figure2 ##
library(stringi)
library(stringr)
library(dplyr)
library(DOSE)
library(org.Hs.eg.db)
data_all<-read.table('FileS2.txt',sep = '\t',header = TRUE)
##
gene_list_us_enrichment<-data_all
disease_gene_list<-read.csv('TableS8 Disorder_gene.csv')
pd_gene<-disease_gene_list$Parkinson.s.Disease
pd_gene<-pd_gene[which(pd_gene!="")]
test = bitr(pd_gene, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENSEMBL",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
PD_data<-gene_list_us_enrichment[which(gene_list_us_enrichment$Disorder.Phenotype=='PD'),]
scz_data<-gene_list_us_enrichment[which(gene_list_us_enrichment$gwas_file=='SCZ'),]

getgene_number<-function(disorder_data){
  
 # disorder_data<-PD_data
  save_name<-disorder_data$Disorder.Phenotype%>%unique()
  disease_name<-disorder_data$Disorder.Phenotype%>%unique()
  gene_id<-'xxx'
  gene_symbol<-"xxx"
  p_value_d<-0
  state_data<-"statename"
  all_gene<-data.frame(disease_name,state_data,gene_id,gene_symbol,p_value_d)
  for (i in c(1:nrow(disorder_data))) {
   # i=1
    tem_gene<-disorder_data[i,4]
    tem_p<-disorder_data[i,6]
    tem_sym<-disorder_data[i,5]
    state_data<-disorder_data[i,]$StateName
    dis_name<-disorder_data[i,]$Disorder.Phenotype
    gene_id<-str_split(tem_gene,'[ ]')%>%unlist()
    p_value_d<-str_split(tem_p,'[ ]')%>%unlist()%>%as.numeric()
    gene_symbol<-str_split(tem_sym,'[ ]')%>%unlist()
    disease_name<-rep(dis_name,length(gene_id))
    state_data<-rep(state_data,length(p_value_d))
    tem_df<-data.frame(disease_name,state_data,gene_id,gene_symbol,p_value_d)
    all_gene<-rbind(all_gene,tem_df)
  }
  all_gene<-all_gene[-1,]
  ##
  colnames(all_gene)<-c('disorder_name','Tissue/Cell_name',"Ensembl_id",'Gene_SYmbol','P-value')
  write.table(all_gene,paste('./',save_name,".txt",sep = ""),quote = FALSE,row.names = FALSE)
}
name_id_dis<-unique(gene_list_us_enrichment$Disorder.Phenotype)
name_id_dis<-name_id_dis[-c(12,3,5,10,13,19)]
for (i in name_id_dis) {
  disorder_data<-gene_list_us_enrichment[which(gene_list_us_enrichment$Disorder.Phenotype==i),]
  getgene_number(disorder_data)
}



pd_inffer<-getgene_number(PD_data)
pd_inffer<-pd_inffer[!duplicated(pd_inffer$V1),]
rownames(pd_inffer)<-pd_inffer$V1
pd_unfind<-setdiff(pd_inffer$V1,pd_gene)
disease_unfind_set<-pd_inffer[pd_unfind,]
disease_unfind_set$V2<-as.numeric(disease_unfind_set$V2)
disease_unfind_set<-disease_unfind_set[order(disease_unfind_set$V2,)]
dis_name<-gene_list_us_enrichment$gwas_file%>%unique()
dis_name<-dis_name[-c(3,5,10,12,13,19)]
#### re organised the data ##
target_data<-data_all[which()]










