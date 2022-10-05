##
library(UpSetR)
library(VennDiagram)
library(export)
library(ggplot2)
library(clusterProfiler)
all_files<-list.files('./')
tem_data<-read.table(all_files[1],sep ='\t',header = TRUE)
#####
for (j in c(2:length(all_files))) {
  
  tem_data_e<-read.table(all_files[j],sep ='\t',header = TRUE)
  tem_data<-rbind(tem_data,tem_data_e)
}
write.csv(tem_data,'disorder_csv.csv',row.names = FALSE,quote = FALSE)


col_name_class<-c("AD_Jansen","ALS","EPILEPSY","MS","PD",
                                   "ADHD","ASD","BIP","SCZ", 
                                "AUD","INTELLIGENCE", "INSOMNIA","Neuroticism","RISK_behavior")
########## neuro #
neuro_dis<-col_name_class[1:5]
psy_dis<-col_name_class[6:9]
cog<-col_name_class[10:14]
####
get_number_gene<-function(dis_clas){
 # dis_clas<-neuro_dis
  gene_set<-c()
  gene_set1<-c()
  re_data<-list()
  for (i in dis_clas) {
    #i=dis_clas[1]
    tem_gene<-read.table(paste(i,'.txt',sep = ''),fill = TRUE,header = TRUE,sep = '\t')
    ###
    gene_set<-c(gene_set,unique(tem_gene$Gene_SYmbol))
    gene_set1<-c(gene_set1,unique(tem_gene$Ensembl_id))
  }
  gene_share<-table(gene_set)%>%as.matrix()
  gene_share<-gene_share[which(gene_share>1),]
  gene_share1<-table(gene_set1)%>%as.matrix()
  gene_share1<-gene_share1[which(gene_share1>1),]
  print(length(gene_share))
  re_data[[1]]<-gene_set

  me_set_gene<-cbind(names(gene_share1),names(gene_share))%>%as.matrix()
  rownames(me_set_gene)<-me_set_gene[,1]
  re_data[[2]]<-me_set_gene
  return(re_data)
}
neuro_dis_count<-get_number_gene(neuro_dis)
psy_dis_count<-get_number_gene(psy_dis)
cog_count<-get_number_gene(cog)
###
write.csv(neuro_dis_count[[2]],'neo.csv',row.names = FALSE,quote = FALSE)
write.csv(psy_dis_count[[2]],'psy.csv',row.names = FALSE,quote = FALSE)
write.csv(cog_count[[2]],'cog.csv',row.names = FALSE,quote = FALSE)
##
gene_neo_d<-neuro_dis_count[[2]]%>%as.matrix()
gene_psy_d<-psy_dis_count[[2]]
gene_cog_d<-cog_count[[2]]
dd<-intersect(gene_neo_d[,1],gene_psy_d[,1])


###
# min=2  50 259 770 ##
### min=3  1 3 132 ##
value_class<-c(50,259) # Behavioral-cognitive 
class_name<-c('Neurological','Psychiatric') 
data_plot<-data.frame(group_name=class_name,number=value_class)
data_plot$group_name<-factor(data_plot$group_name,levels =class_name )
ggplot(data_plot,aes(x=group_name,y=number))+geom_bar(stat = 'identity')
##
##    ##
Reduce(intersect,list(neuro_dis_count[[2]],psy_dis_count[[2]],cog_count[[2]]))
upset_list <- list(neuro_dis_count[[2]],psy_dis_count[[2]],cog_count[[2]])   # 件
names(upset_list) <- c('Neurological','Psychiatric','Behavioral-cognitive ')    #
upset(fromList(upset_list),  # 
      nsets = 100,  
      main.bar.color = "brown",# 
      nintersects = 40, #
      order.by = "freq", #
      keep.order = F, # 
      mb.ratio = c(0.6,0.4), 
      sets.bar.color = ggpubr::get_palette('npg',dim(fromList(upset_list))[2]),#
      text.scale = 2 # 
)
##
write.table(neuro_dis_count[[2]],'Neurological.txt',row.names = FALSE,quote = FALSE,col.names = FALSE)
write.table(psy_dis_count[[2]],'Psychiatric.txt',row.names = FALSE,quote = FALSE,col.names = FALSE)
write.table(cog_count[[2]],'Behavioral.txt',row.names = FALSE,quote = FALSE,col.names = FALSE)
###
graph2ppt(file='up_set.pptx',width=8,height=6,append=TRUE)
####fisher test ##
all_gene_list<-read.table('E:/downLoad/NCBI37.3/NCBI37.3.gene.loc')
all_gene_data<-all_gene_list$V6
#####

fisher_test_f<-function(set1,set2){
  a<-intersect(set1,set2)%>%length()
  b<-setdiff(set1,set2)%>%length()
  c<-setdiff(set2,set1)%>%length()
  d<-setdiff(all_gene_data,union(set1,set2))%>%length()
  fisher_ma<-rbind(c(a,b),c(c,d))
  print(fisher_ma)
  print(fisher.test(fisher_ma))
}
fisher_test_f(neuro_dis_count[[2]],psy_dis_count[[2]])
fisher_test_f(neuro_dis_count[[2]],cog_count[[2]])
fisher_test_f(cog_count[[2]],psy_dis_count[[2]])
#### share pathway analysis #
Neurological_path<-read.table('Neurological_path.txt',fill = TRUE,header = TRUE,sep = '\t')
Psychiatric_path<-read.table('Psychiatric_path.txt',fill = TRUE,header = TRUE,sep = '\t')
Behavioral_path<-read.table('Behavioral_path.txt',fill = TRUE,header = TRUE,sep = '\t')
########################## Analysis the result ######################
Neurological_bp<-Neurological_path[which(Neurological_path$Category=='GO: Biological Process'),]
Psychiatric_bp<-Psychiatric_path[which(Psychiatric_path$Category=='GO: Biological Process'),]
Behavioral_bp<-Behavioral_path[which(Behavioral_path$Category=='GO: Biological Process'),]
##########################
intersect(Neurological_bp[1:,]$Name,Psychiatric_bp[1:80,]$Name)
intersect(Behavioral_bp[1:3,]$Name,Psychiatric_bp[1:3,]$Name)
intersect(Behavioral_bp[1:3,]$Name,Neurological_bp[1:3,]$Name)
########
Neurological_kegg<-Neurological_path[which(Neurological_path$Source=='BioSystems: KEGG'),]
Psychiatric_kegg<-Psychiatric_path[which(Psychiatric_path$Source=='BioSystems: KEGG'),]
Behavioral_kegg<-Behavioral_path[which(Behavioral_path$Source=='BioSystems: KEGG'),]
##
intersect(Neurological_kegg[1:3,]$Name,Psychiatric_kegg[1:3,]$Name)
intersect(Behavioral_kegg[1:3,]$Name,Psychiatric_kegg[1:3,]$Name)
intersect(Behavioral_kegg[1:3,]$Name,Neurological_kegg[1:3,]$Name)
###### manager the pathway data ##
bp_data<-rbind(Neurological_bp[1:5,],Psychiatric_bp[1:5,],Behavioral_bp[1:5,])
##
write.csv(bp_data,'bp_data.csv')
group_bp_data<-c(rep('Neurological',3),rep('Psychiatric',3),rep('Behavioral-cognitive',3))
kegg_data<-rbind(Neurological_kegg[1:3,],Psychiatric_kegg[1:3,],Behavioral_kegg[1:3,])
path_name<-c(bp_data$Name,kegg_data$Name)
gene_ratio<-c(bp_data$Hit.Count.in.Query.List/as.numeric(bp_data$Hit.Count.in.Genome),
              kegg_data$Hit.Count.in.Query.List/as.numeric(kegg_data$Hit.Count.in.Genome))
fdr_bh<-c(bp_data$q.value.FDR.B.H,kegg_data$q.value.FDR.B.H)
group_sca<-rep(group_bp_data,2)
###
plot_df<-data.frame(path_name=path_name,group=group_sca,gene_ratio=gene_ratio,FDR=fdr_bh)
plot_df$path_name<-factor(plot_df$path_name,levels = rev(unique(plot_df$path_name)))
###
## plot  data
ggplot(plot_df, aes(group, path_name)) +
  geom_point(aes(color=fdr_bh,size = 3))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,vjust=0.5))+scale_color_gradient(low='#6699CC',high='#CC3333')
##   
graph2ppt(file='path.pptx',width=8,height=6,append=TRUE)

