#expression data#
library(gdata)
library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(export)
library(ggalt)
library(gdata)
library(preprocessCore)
library(forecast)
library(tidyverse)
library(officer)
library(rvg)
require(ggsci)
require(ggpubr)
library(VennDiagram)
gene_list_us_enrichment<-read.table('state-signature genes_id.overlap.zscore.50%.fdr(1).txt',fill = TRUE,sep = '\t',header = TRUE  )
#rownames(gene_list_us_enrichment)<-gene_list_us_enrichment$state_name
#gene_list_us_fdr<-read.table('cell-signature genes_id.overlap.zscore.50%.fdr.txt',fill = TRUE,sep = '\t',header = TRUE  )
expression_data<-fread('./promoter_signal.matrix')
all_gene_id<-expression_data$promoter_id%>%as.matrix()
#
all_gene_id_ms<-apply(all_gene_id, 1, function(x){
  x<-str_split(x,'[|]')%>%unlist()
  return(x[1])
})
all_gene<-apply(all_gene_id, 1, function(x){
  x<-str_split(x,'[|]')%>%unlist()
  return(x[2])
})
#gene mapping #
gene_mapping<-data.frame(id=all_gene_id_ms,gene_sym=all_gene)
rownames(gene_mapping)<-gene_mapping[,1]
gene_mapping<-gene_mapping[which(!duplicated(gene_mapping$gene_sym)),]
rownames(gene_mapping)<-gene_mapping$id
expression_data<-expression_data[,-1]
expression_data<-t(scale(t(expression_data),center = FALSE))
expression_data<-as.matrix(expression_data)
rownames(expression_data)<-all_gene_id_ms

##
tissue_end<-c()
for (i in using_region_new) {
 # i=tissue_name[1]
  tem_data<-str_split(i,'[,]')%>%unlist()
  tissue_end<-c(tissue_end,tem_data[1])
}
## get the same brain region #
tissue_end<-as.matrix(table(tissue_end))
brain_region<-names(tissue_end[which(tissue_end>1),])
#
adult_region<-paste(tissue_end,', adult',sep = "")
new_born_region<-paste(brain_region,', newborn',sep = "")
###
using_region_new<-intersect(new_born_region,newborn_data)
using_region_adult<-intersect(adult_region,adult_data)
#
using_tissue<-c(using_region_new,using_region_adult)
using_tissue<-c(adult_data,newborn_data)
##
get_sample_id<-function(meta_data,using_tissue){
  re_matrix<-list()
 # rownames(re_matrix)<-using_tissue
  for (i in c(1:length(using_tissue))) {
    #i<-1
    re_matrix[[using_tissue[i]]]<-meta_data[which(meta_data$AggregatedBiosampleName==using_tissue[i]),2]
  }
  return(re_matrix)
}
tissue_using_id<-get_sample_id(meta_data,using_tissue)
##our expression data ##
disease_class<-cbind(annote_files_col,rownames(annote_files_col))
## combine the gene of every state #
tissue_type<-using_tissue
get_every_disorder_geneSet<-function(gene_list_us_enrichment,tissue_type){
  disorder_data<-gene_list_us_enrichment$gwas_file%>%unique()
  disorder_list<-list()
  for (i in c(1:length(disorder_data))) {
   # i=1
    tem_dis<-gene_list_us_enrichment[which(gene_list_us_enrichment$gwas_file==disorder_data[i]),]
    rownames(tem_dis)<-tem_dis$state_name
    tem_dis<-tem_dis[tissue_type,]
    gene_name<-c()
    for (j in c(1:nrow(tem_dis))) {
    # j=1
      gene_s<-str_split(tem_dis[j,4],'[ ]')%>%unlist()
      gene_name<-c(gene_name,gene_s)
      
    }
    gene_name<-unique(gene_name)
    if(length(gene_name)>1){
      disorder_list[[disorder_data[i]]]<-gene_name
    }
   
  }
  return(disorder_list)
}
get_different_intersection<-function(tissue1,tissue2){
  dis_list<-list()
  gwas_list<-gene_list_us_enrichment$gwas_file%>%unique() 
  save_ratio<-matrix(NA,length(gwas_list),1)
  rownames(save_ratio)<-gwas_list
  k=1
  for ( dis in gwas_list) {
    #dis<-gwas_list[1]
    inter_gene_dis<-get_every_disorder_geneSet(gene_list_us_enrichment,using_tissue[1])
    inter_gene_dis<-inter_gene_dis[[dis]]
  for (i in c(2:length(using_tissue))) {
     disorder_list_1<-get_every_disorder_geneSet(gene_list_us_enrichment,using_tissue[i])
     disorder_list_1<-disorder_list_1[[dis]]
     inter_gene_dis<-intersect(inter_gene_dis,disorder_list_1)
  }
    dis_list[[dis]]<-inter_gene_dis
    disorder_list_2<-get_every_disorder_geneSet(gene_list_us_enrichment,using_tissue)
    #disorder_list_2[[dis]]
    save_ratio[k,1]<-length(inter_gene_dis)/length(disorder_list_2[[dis]])
    k=k+1
  }
}
####



#
disease_age_group_expression<-function(disorder_type){
 # disorder_type<-'psychiatrics'
  disease_cla<-disease_class[which(disease_class==disorder_type),]%>%rownames()
  if(disorder_type=='Neurological'){
    disease_cla<-disease_cla[-4] 
  }
  for (i in c(1:length(disease_cla))) {
   # i=3
    print(i)
    tem_name_dis<-disease_cla[i]
    if(tem_name_dis=="AD_JENSE"){
      tem_name_dis<-'AD_Jansen'
    }else if(tem_name_dis=="Parkinson"){
      tem_name_dis<-"PD"
    }else if(tem_name_dis=="ANXIETY"){
      
      tem_name_dis<-"Anxiety"
    }
    ##
   # gene_stata_data<-gene_list_us_enrichment[which(gene_list_us_enrichment$gwas_file==tem_name_dis),]
    gene_stata_data<-disorder_list[[tem_name_dis]]
    if(length(gene_stata_data)>1){
    raw_expression<-expression_data[,1]
    tissue_vec<-vector()
    for (j in tissue_type) {
     # j=using_tissue[1]
      using_expression<-expression_data[gene_stata_data,meta_data[which(meta_data$AggregatedBiosampleName==j),2]]%>%as.matrix()
      if(ncol(using_expression)>1){
        using_expression<-apply(using_expression, 1, mean)
      }
      tissue_vec<-c(tissue_vec,rep(j,length(using_expression)))
      raw_expression<-cbind(raw_expression,using_expression)
    } 
    raw_expression<-raw_expression[,-1]
    colnames(raw_expression)<-tissue_type
###########
  ###
     new_born_exp<-unmatrix(raw_expression[,1:14],byrow=FALSE)
     adult_exp<-unmatrix(raw_expression[,15:28])
     all_exp<-c(new_born_exp,adult_exp)
    group_name<-c(rep('adult',length(new_born_exp)),rep('new_born',length(adult_exp)))
    tissue_vec<-apply(as.matrix(tissue_vec), 1, function(x){
     # x<-tissue_vec[1]
      x<-str_split(x,'[,]')%>%unlist()
      x<-x[1]
      return(x)
    })
    plot_df_data<-data.frame(expression_data=all_exp,tissue_group=tissue_vec,age_group=group_name)
    rownames(plot_df_data)<-NULL
  
   # ggplot(plot_df_data,aes(x=tissue_group,y=expression_data))+
   #  geom_boxplot(aes(fill=age_group),outlier.colour = NA)+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))+
   #  theme(panel.grid =element_blank())+
   #  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
   #  ylim(0,5)+
   #  theme(panel.border = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=12),
   #        axis.text.y = element_text(size=15)) +   ## 删去外层边框
   #  theme(axis.line = element_line(size=1, colour = "black"))+
   #  geom_signif(comparisons = list(tissue_type),map_signif_level=TRUE)
   # 
   
   p <- ggboxplot(plot_df_data, x = "tissue_group", y = "expression_data",
                  color = "age_group", palette = "jco",outlier.colour = NA)+ ylim(0,5)+
     theme(panel.border = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=12),
           axis.text.y = element_text(size=15))+
     ggtitle(tem_name_dis)+theme(plot.title = element_text(hjust = 0.5))
     
   p<-p + stat_compare_means(aes(group = age_group),label = "p.signif")
   print(p)
  graph2ppt(file="tissue_type_data_age_group", width=16, height=9,aspectr=sqrt(2), append=TRUE)
    }
  }
}
for (i in c('Neurological','psychiatrics','Trait')) {
  disease_age_group_expression(i)
}
disease_age_group_expression
###

age_data<-adult_data
disease_class_expression<-function(disease_type,age_data){
  
  disease_cla<-disease_class[which(disease_class==disease_type),]%>%rownames()
 # disease_cla<-rownames(disease_class)
  if(disease_type=='Neurological'){
    disease_cla<-disease_cla[-4]
  }
  return_data<-matrix(NA,1,3)
  for (i in c(1:length(disease_cla))) {
  #  i=1
    print(i)
    tem_name_dis<-disease_cla[i]
    if(tem_name_dis=="AD_JENSE"){
      tem_name_dis<-'AD_Jansen'
    }else if(tem_name_dis=="Parkinson"){
      tem_name_dis<-"PD"
    }else if(tem_name_dis=="ANXIETY"){
      
      tem_name_dis<-"Anxiety"
    }
    ## gene set plot #
    gene_stata_data<-gene_list_us_enrichment[which(gene_list_us_enrichment$gwas_file==tem_name_dis),]
    rownames(gene_stata_data)<-gene_stata_data$state_name
    gene_stata_data<-gene_stata_data[age_data,]
   #gene_stata_data<-gene_stata_data[tissue_type,]
    #### 针对每个state 取出对应的基因，然后得出对应的基因的表达值# 
    for (j in c(1:nrow(gene_stata_data))) {
    # print(j)
      #j=4
      tem_expression<-gene_stata_data[j,]
      gene_sen_te<-str_split(tem_expression$signature.genes,'[ ]')%>%unlist()
      if(gene_sen_te!=''){
        sample_id<-meta_data[which(meta_data$AggregatedBiosampleID==tem_expression$state_id),2]
        gene_stat_expression<-expression_data[gene_sen_te,sample_id]
        if(length(gene_sen_te)==1){
          gene_stat_expression<-mean(gene_stat_expression)
        }
        else if(length(sample_id)>1){
          gene_stat_expression<-apply(gene_stat_expression, 1, mean)
        }
        cell_name<-c(rep(tem_expression$state_name,length(gene_sen_te)))
        disease_name<-c(rep(disease_cla[i],length(gene_sen_te)))
        com_data<-cbind(gene_stat_expression,cell_name,disease_name)
        return_data<-rbind(return_data,com_data)
      }

    }

  }
  return_data<-return_data[-1,]
  re_df<-data.frame(expression_value=as.numeric(return_data[,1]),Tissue_name=return_data[,2],disorder_name=return_data[,3])
  re_df<-re_df[which(re_df$disorder_name!='ALS'),]
  p<-re_df%>%mutate(disorder_name = fct_reorder(disorder_name, expression_value, .fun='median'))%>%ggplot(aes(x=Tissue_name,y=expression_value))+
    geom_boxplot(aes(fill=disorder_name),outlier.colour = NA)+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))+
    theme(panel.grid =element_blank())+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=10),
          axis.text.y = element_text(size=10))+
    theme(panel.border = element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))+ylim(0,5)+
    ylab('Normalized expression')+xlab("Tissue name")
  print(p)
  graph2ppt(file="disorders_region_data.pptx", width=16, height=9,aspectr=sqrt(2), append=TRUE)
}
for (i in c('Neurological','psychiatrics','Trait')) {
  disease_class_expression(i,newborn_data)
  disease_class_expression(i,fetal_data)
  disease_class_expression(i,adult_data)
}

## enrichment value ##
disease_enrichment<-function(){
  
  disease_cla<-disease_class[which(disease_class=='Neurological'),]%>%rownames()
  disease_cla<-disease_cla[-4]
  return_data<-matrix(NA,1,3)
  for (i in c(1:length(disease_cla))) {
     #i=4
    print(i)
    tem_name_dis<-disease_cla[i]
    if(tem_name_dis=="AD_JENSE"){
      tem_name_dis<-'AD_Jansen'
    }else if(tem_name_dis=="Parkinson"){
      tem_name_dis<-"PD"
    }else if(tem_name_dis=="ANXIETY"){
      
      tem_name_dis<-"Anxiety"
    }
    gene_stata_data<-gene_list_us_enrichment[which(gene_list_us_fdr$gwas_file==tem_name_dis),]
    for (j in c(1:nrow(gene_stata_data))) {
      # print(j)
      #j=4
      tem_expression<-gene_stata_data[j,]
      gene_sen_te<-str_split(tem_expression$P,'[ ]')%>%unlist()
      if(gene_sen_te!=''){
        cell_name<-c(rep(tem_expression$state_name,length(gene_sen_te)))
        disease_name<-c(rep(disease_cla[i],length(gene_sen_te)))
        com_data<-cbind(as.numeric(gene_sen_te),cell_name,disease_name)
        return_data<-rbind(return_data,com_data)
      }
    }
    
  }
  return_data<-return_data[-1,]
  re_df<-data.frame(enrichment=as.numeric(return_data[,1]),cell_name=return_data[,2],disorder_name=return_data[,3])
  ##
  re_df$enrichment<-scale(re_df$enrichment,center = FALSE)
  ##
  ggplot(re_df,aes(x=disorder_name,y=enrichment))+
    geom_boxplot(aes(fill=cell_name))+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))+
    theme(panel.grid =element_blank())+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
    theme(panel.border = element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))
  graph2ppt(file="gene_enrichment_data", width=16, height=9,aspectr=sqrt(2), append=TRUE)
  
}
############ plot the disease intersection ####
disorder_list<-read.csv('disease_gene.csv')

disease_gene_all<-get_every_disorder_geneSet(gene_list_us_enrichment,using_tissue)
##ADHD
ADHD_gene<-disorder_list$ADHD
ASD_gene<-disorder_list$Autism.Spectrum.Disorders.ASD.
SCZ_gene<-disorder_list$Schizophrenia.
Bip_gene<-disorder_list$Bipolar
AD_gene<-disorder_list$Alzheimer..AD.
PD_gene<-disorder_list$Parkinson
####  
ADHD_gene<-ADHD_gene[which(ADHD_gene!="")]
ASD_gene<-ASD_gene[which(ASD_gene!="")]
SCZ_gene<-SCZ_gene[which(SCZ_gene!="")]
Bip_gene<-Bip_gene[which(Bip_gene!="")]
PD_gene<-PD_gene[which(PD_gene!="")]
AD_gene<-AD_gene[which(AD_gene!="")]
#out
our_gene_ASD<-gene_mapping[disease_gene_all[['ASD']],2]
our_gene_SCZ<-gene_mapping[disease_gene_all[['SCZ']],2]
our_gene_ADHD<-gene_mapping[disease_gene_all[['ADHD']],2]
our_gene_Bip<-gene_mapping[disease_gene_all[['BIP']],2]
our_gene_AD<-gene_mapping[disease_gene_all[['AD_Jansen']],2]
our_gene_PD<-gene_mapping[disease_gene_all[['PD']],2]
fisher_test_data<-function(genelist1,genelist2,disease_name){
  # genelist1<-our_gene_ASD
  # genelist2<-ASD_gene
  genelist1<-genelist1[which(!is.na(genelist1))]
  genelist2<-genelist2[which(!is.na(genelist2))]
  all_gene_list<-gene_mapping[,2]
  a<-intersect(genelist1,genelist2)%>%length()
  b<-setdiff(genelist1,a)%>%length()
  c<-setdiff(genelist2,a)%>%length()
  d<-setdiff(all_gene_list,c(genelist1,genelist2))%>%length()
  venn.diagram(
    x = list(genelist1,genelist2),
    category.names = c("out_res" , "Published_res"),
    filename = paste('IMG/',disease_name,'venn.png',sep = ""),
    output = TRUE ,
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    lwd = 1,
    col=c("#440154ff", '#21908dff'),
    fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
    cex = 0.5,
    fontfamily = "sans",
    cat.cex = 0.3,
    cat.fontfamily = "sans",
    cat.col = c("#440154ff", '#21908dff'),
  )
  fis_table<-rbind(c(a,b),c(c,d))
  print(fis_table)
  fisher.test(fis_table)
}
fisher_test_data(our_gene_SCZ,SCZ_gene,'SCZ')
fisher_test_data(our_gene_ADHD,ADHD_gene,'ADHD')
fisher_test_data(our_gene_Bip,Bip_gene,'BIP')
fisher_test_data(our_gene_ASD,ASD_gene,'ASD')
fisher_test_data(our_gene_AD,AD_gene,'AD')
fisher_test_data(our_gene_PD,PD_gene,'PD')

###plot the  frontal ##

###
all_disorder_type<-c('Neurological','psychiatrics','Trait')
plot_brain_region_diff<-function(sample_name){
  #sample_name<-test_data
  return_data<-data.frame(exp_value=0,disorder_name="sas",disorder_type='sadas')
  for (j in c(1:length(all_disorder_type))) {
   # i=3
    disease_type<-all_disorder_type[j]
    disease_cla<-disease_class[which(disease_class==disease_type),]%>%rownames()
    # disease_cla<-rownames(disease_class)
    if(disease_type=='Neurological'){
      disease_cla<-disease_cla[-c(3,4)]
    }
    for (i in c(1:length(disease_cla))) {
      #i=1
      print(i)
      tem_name_dis<-disease_cla[i]
      if(tem_name_dis=="AD_JENSE"){
        tem_name_dis<-'AD_Jansen'
      }else if(tem_name_dis=="Parkinson"){
        tem_name_dis<-"PD"
      }else if(tem_name_dis=="ANXIETY"){
        
        tem_name_dis<-"Anxiety"
      }
      ## gene set plot #
      gene_stata_data<-gene_list_us_enrichment[which(gene_list_us_enrichment$gwas_file==tem_name_dis),]
      rownames(gene_stata_data)<-gene_stata_data$state_name
      gene_stata_data<-gene_stata_data[sample_name,]
      #gene_stata_data<-gene_stata_data[tissue_type,]
      #### 针对每个state 取出对应的基因，然后得出对应的基因的表达值# 
      for (j in c(1:nrow(gene_stata_data))) {
        # print(j)
       # j=1
        tem_expression<-gene_stata_data[j,]
        gene_sen_te<-str_split(tem_expression$signature.genes,'[ ]')%>%unlist()
        if(gene_sen_te!=''){
          sample_id<-meta_data[which(meta_data$AggregatedBiosampleID==tem_expression$state_id),2]
          gene_stat_expression<-expression_data[gene_sen_te,sample_id]
          if(length(gene_sen_te)>1){
            #median  
            gene_stat_expression<-median(gene_stat_expression)
          }
          #xp_value=0,disorder_name="sas",disorder_type='sadas'
          com_data<-data.frame(exp_value=gene_stat_expression,disorder_name=tem_name_dis,disorder_type=disease_type)
          return_data<-rbind(return_data,com_data)
        }
      }
    }
  }
  return_data<-return_data[-1,] 
  id<-c(1:nrow(return_data))
  return_data<-cbind(return_data,id)
  label_data <- return_data
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  p <- ggplot(return_data, aes(x=as.factor(id), y=exp_value , fill=disorder_type )) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(stat="identity", alpha=0.5) +
    ylim(0,5)+
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5,color="black"),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +coord_polar()+
    geom_text(data=label_data, aes(x=id, y=exp_value+0.1, label=disorder_name, hjust=hjust), color="black", fontface="bold",alpha=0.7, size=3, angle= label_data$angle, inherit.aes = FALSE ) 
  
  print(p)
 # graph2ppt(file="special_region_data.pptx", width=16, height=9,aspectr=sqrt(2), append=TRUE)
}
###
test_data<-c('hippocampus, newborn','hippocampus, adult')


##cerebellum, adult"              "cerebellum, newborn" 
plot_brain_region_diff(c('cerebellum, adult, fetal','cerebellum, newborn'))
plot_brain_region_diff(c('occipital lobe, fetal','occipital cortex, newborn','occipital lobe, adult'))
plot_brain_region_diff(c('hippocampus, newborn','hippocampus, adult'))

plot_brain_region_diff(c('temporal lobe, fetal','temporal lobe, adult'))
plot_brain_region_diff(c('brain, fetal','brain, adult'))
plot_brain_region_diff(c('parietal lobe, adult','parietal lobe, fetal'))
plot_brain_region_diff(c('substantia nigra, newborn','substantia nigra, adult'))


### tissue expression data  from the new born and adult sample ##
ying_expression<-fread('./brain_development_RPKM.txt',header = FALSE)

ying_id<-ying_expression$V1%>%as.matrix()
ying_id<-apply(ying_id, 1, function(x){
  x<-str_split(x,'[|]')%>%unlist()
  return(x[2])
})

ying_expression<-ying_expression[,-1]
#
ying_new_data<-data.frame(data=ying_expression,symbol=ying_id)
ying_new_data<-ying_new_data[!duplicated(ying_new_data$symbol),]
rownames(ying_new_data)<-ying_new_data$symbol
ying_new_data<-ying_new_data[,-ncol(ying_new_data)]
sampple_infpr<-read.csv('RNA-sample.csv',header = TRUE)
## chosing the human  samples   ##
human_sample<-sampple_infpr[which(sampple_infpr$Species=='Human'),]##460  sample #
age_group_ying<-human_sample$Period
colnames(ying_new_data)<-sampple_infpr$Sample
ying_new_data<-ying_new_data[,human_sample$Sample]
ying_new_data_scale<-normalize.quantiles(as.matrix(ying_new_data))
rownames(ying_new_data_scale)<-rownames(ying_new_data)
colnames(ying_new_data_scale)<-colnames(ying_new_data)
###
pfc_exp<-human_sample[which(human_sample$Species=="Human"),1]

disorder_type<-'Trait'
plot_data<-function(ying_new_data_scale,disorder_type,disorder_gene_set){
  disease_cla<-disease_class[which(disease_class==disorder_type),]%>%rownames()
  using_diorder<-intersect(disease_cla,names(disorder_gene_set))
  if(disorder_type=='Neurological'){
    using_diorder<-c(using_diorder,'AD_Jansen','PD')
  }
  # return_data<-data.frame(mean_value=0,left_cond=0,
  #                         right_cond=0,period=2,disorder_type="test")
  return_data<-data.frame(mean_value=0,age_period=2,disorder_type="test")
  save_group_plot_data<-data.frame(Normalization_exp=0,age_period=2,disorder_type="test")
  for (i in c(1:length(using_diorder))) {
    #print(i)
   # i=4
    tem_gene<-disorder_gene_set[[using_diorder[i]]]
    tem_gene<-gene_mapping[tem_gene,2]
    temp_expression<-ying_new_data_scale[intersect(tem_gene,rownames(ying_new_data_scale)),]

    no_born<-human_sample[which(human_sample$Period==8|human_sample$Period==9),1]
    born_data<-human_sample[which(human_sample$Period>=13),1]
    
    #no_born<-human_sample[which(human_sample$Period<8),1]
    #born_data<-human_sample[which(human_sample$Period>=8),1]
    #no_born_exp<-apply(temp_expression[,no_born],1,mean)
    #born_exp<-apply(temp_expression[,born_data],1,mean)
    no_born_exp<-unmatrix(temp_expression[,intersect(no_born,pfc_exp)])
   born_exp<-unmatrix(temp_expression[,intersect(born_data,pfc_exp)])
    test_t<-t.test(no_born_exp,born_exp)
    if(test_t$p.value<=0.05){
      print(test_t$p.value)
      print(using_diorder[i])
    }
# different group data # 
    save_re<-cbind(c(no_born_exp,born_exp),c(rep('Prenatal',length(no_born_exp)),rep("Postnatal",length(born_exp))),rep(using_diorder[i],length(c(no_born_exp,born_exp))))
    colnames(save_re)<-colnames(save_group_plot_data)
    save_group_plot_data<-rbind(save_group_plot_data,save_re)
    ##divide the group
    for (j in unique(age_group_ying)) {
     #j=age_group_ying[1]
      age_tem_sam<-human_sample[which(human_sample$Period==j),1]
      if(length(intersect(pfc_exp,age_tem_sam))>0){
        age_tem_expression<-temp_expression[,intersect(pfc_exp,age_tem_sam)]%>%as.matrix()
        age_tem_expression_mean<-apply(age_tem_expression, 2, mean)
       # age_tem_expression_mean<-unmatrix(age_tem_expression,byrow = TRUE)
        raw_data<-data.frame(mean_value=age_tem_expression_mean,age_period=rep(j,length(age_tem_expression_mean)),
                             disorder_type=rep(using_diorder[i],length(age_tem_expression_mean)))
        return_data<-rbind(return_data,raw_data)
      }
    }
    
  }
  return_data<-return_data[-1,]
  save_group_plot_data<-save_group_plot_data[-1,]
  #rownames(save_group_plot_data)<-NULL
 # return_data$period<-levels(as.character(return_data$period))
  return_data$age_period <- factor(return_data$age_period,levels = unique(return_data$age_period))
  ggplot(data=return_data, 
         aes(x=age_period, y=mean_value, group=disorder_type, color=factor(disorder_type),
             fill=factor(disorder_type)))+ geom_smooth(method='loess')+
    xlab('periods') +ylab("Normalizated expression")
  
  #
  color1<-c('#DE47AB','#DE47AB','#72BE97','#72BE97','#E85726')
  chr1 <- '#DE47AB'
  chr2 <- '#72BE97'
  chr3 <- '#F7F797'
  chr4 <- '#7C749B'
  chr5 <- '#E85726'
  chr6 <- '#B395F8'
  chr7 <- '#DC8747'
  chr8 <- '#96D53D'
  # save_group_plot_data$Normalization_exp<-as.numeric(save_group_plot_data$Normalization_exp)
  # p<-ggplot(save_group_plot_data,aes(x=disorder_type,y=Normalization_exp,ylim()))+
  #   geom_boxplot(aes(fill=age_period))+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))+
  #   theme(panel.grid =element_blank())+
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
  #   theme(panel.border = element_blank()) +   ## 删去外层边框
  #   theme(axis.line = element_line(size=1, colour = "black"))
  # print(p)
  #graph2ppt(file="expression_data_age_group", width=16, height=9,aspectr=sqrt(2), append=TRUE)
  #####
 # graph2ppt(file="gene_enrichment_data", width=16, height=9,aspectr=sqrt(2), append=TRUE)
  # ggplot(return_data, aes(x=age_period, y=mean_value, group=disorder_type, color=factor(disorder_type))) +
  #   #geom_point(color="black") +
  #   geom_smooth(se=FALSE, linetype="dashed", size=0.5) +
  #   geom_xspline(spline_shape=0.4, size=0.5)
}
new_born_disorder<-get_every_disorder_geneSet(gene_list_us_enrichment,tissue_type[1:14])
adult_gene_set<-get_every_disorder_geneSet(gene_list_us_enrichment,tissue_type[14:28])
disorder_gene_set<-get_every_disorder_geneSet(gene_list_us_enrichment,tissue_type)
#fetal state
disorder_gene_set<-get_every_disorder_geneSet(gene_list_us_enrichment,c('brain, fetal','brain, fetal'))
#psychiatrics  Trait Neurological
plot_data(ying_new_data_scale,'Neurological',disorder_gene_set)
  ###########



library(tidyverse)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggpubr)
# create a dataset
plot_data<-exp_spec_m[,c(3,4,5)]
plot_data$Freq<-as.character(plot_data$Freq)
##
plot_data$Freq <- factor(plot_data$Freq, levels = c("0","1","2","3","4","5+"))
# Plot
ggboxplot(plot_data, x="Freq", y="value", palette = "jco") + ylim(0,40)+
  stat_compare_means(label="p.signif", method="t.test",
                     ref.group = "0.5")

ggboxplot(plot_data, x = "Freq", y = "value",
          color = "Freq", palette = "jco",outlier.shape =NA)+ylim(-12,-7)+
          stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                  # Pairwise comparison against all
##
ggboxplot(plot_data, x = "reg", y = "value",
          color = "reg", palette = "jco",outlier.shape =NA)+ylim(-12.5,-7)+
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                  # Pairwise comparison against all
##












    
    