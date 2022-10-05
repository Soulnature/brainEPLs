# mental disorder enrichment plot #
library(export)
library(stringi)
library(stringr)
library(gdata)
library(Matrix)
library(gdata)
library(RColorBrewer) 
disorders_num<-c('ASD','ALCDEP','BIP','SCZ','PTSD','AUD','OCD','ALCDEP','EPILEPSY','Parkinson','TS','AD_JANSEN','AD_KUNKLE','MS',
   'ANXIETY','ASB','INSOMNIA','INTELLIGENCE','RISK_behavior','Neuroticism','ALS')
all_filr1<-list.files('./disorder_res_LDSC/')
all_filr<-c()
for (i in all_filr1) {
  wheher<-regexpr('results',i)
  if(wheher>0){
    all_filr<-c(all_filr,i)
  }
}
library(pheatmap)
all_name<-vector()

for (i in all_filr) {
  #i=all_filr[1]
  temp<-str_split(i,'_')%>%unlist()
  all_name<-c(all_name,temp[1])
}
all_name<-unique(all_name)
all_name[1]<-'AD_JENSE'
all_name[17]<-'RISK_behavior'
all_name<-c(all_name,'AD_KUNKLE')
all_name<-all_name[-c(20,16,3)]
## 55  states #
states<-c(31:85)
enhancer_data<-matrix(NA, nrow = length(all_name), ncol = 55)
example_data<-read.table(paste(dis_name,'_State_0',as.character(states[1]),'..results',sep = ""),header = TRUE)
example_data<-example_data[1,c(2:ncol(example_data))]
state_cell=states[1]
example_data<-cbind(dis_name,state_cell,example_data)
#promoter_data<-matrix(NA, nrow = length(all_name), ncol = 55)
rownames(enhancer_data)<-all_name
#rownames(promoter_data)<-all_name
colnames(enhancer_data)<-as.character(c(31:85))
#colnames(promoter_data)<-as.character(c(31:85))
for(j in c(1:length(all_name))){
  dis_name=all_name[j]
  for (i in c(1:55)) {
   # dis_name<-all_name[1]
    #i=1
    state_cell<-nama_tissue[i]
    tem_files<-paste(dis_name,'_State_0',as.character(states[i]),'..results',sep = "")
    read_data<-read.table(tem_files,header = TRUE)
    read_data<-cbind(dis_name,state_cell,read_data[1,2:ncol(read_data)])
    example_data<-rbind(example_data,read_data)
   # enhancer_data[j,i]<-read_data$Enrichment_p[1]
      #pnorm(read_data$Coefficient_z.score[1])
  }
}
example_data<-example_data[-1,]
write.csv(example_data,'disease_ldsc.csv',quote = FALSE,row.names = FALSE)
##state name load #
meta_data<-read.csv('./metedata.csv')
state_name_d<-meta_data[which(meta_data$BiosampleGroup=='Brain'),]
nama_tissue<-state_name_d$AggregatedBiosampleName%>%unique()
colnames(enhancer_data)<-nama_tissue
###
enhancer_one<-unmatrix(enhancer_data,byrow=T)
enhancer_p_adjust<-p.adjust(enhancer_one,method='BH')
enhancer_p_adjust_re<- array(enhancer_p_adjust[1:length(enhancer_p_adjust)], dim=c(55,20))
#### result  ##
rownames(enhancer_p_adjust_re)<-nama_tissue
colnames(enhancer_p_adjust_re)<-all_name
enhancer_data_log<--log10(enhancer_p_adjust_re)
######  staistic  the number  of significant p  #####
disorder_tissue_name<-function(enhancer_data_log)
{
  names_disorder<-colnames(enhancer_data_log)
  names_tissue<-rownames(enhancer_data_log)
  vec_num_tissye<-c()
  vec_num_disorder<-c()
  for (i in c(1:nrow(enhancer_data_log))) {
    #i=1
    temp_num<-enhancer_data_log[i,]
    tem_tis_num<-temp_num[which(temp_num>1.3)]%>%length()
    vec_num_tissye<-c(vec_num_tissye,tem_tis_num)
  }
  for (i in c(1:ncol(enhancer_data_log))) {
    #i=1
    temp_num<-enhancer_data_log[,i]
    tem_dis_num<-temp_num[which(temp_num>1.3)]%>%length()
    vec_num_disorder<-c(vec_num_disorder,tem_dis_num)
  }
  disorder_df<-data.frame(region_name=names_disorder,significant_disoerser_num=vec_num_disorder)
  Tissue_df<-data.frame(region_name=names_tissue,significant_tissue_num=vec_num_tissye)
  return(list(disorder_df,Tissue_df))
}
library(ggplot2)
disorder_tissue_df<-disorder_tissue_name(enhancer_data_log)
####
disorder_df<-disorder_tissue_df[[1]]
Tissue_df<-disorder_tissue_df[[2]]
##plot the number
ggplot(data = Tissue_df,mapping = aes(x = region_name, y = significant_tissue_num)) + geom_bar(stat = 'identity')+theme(axis.text.x = element_text(size = 8, vjust =1, hjust = 1, angle = 45))
#+theme(plot.title = element_text(hjust = 0.5,size=18)) 
##
disorder_name<-disorder_df[which(disorder_df$significant_disoerser_num>=1),1]
tissue_name<-Tissue_df[which(Tissue_df$significant_tissue_num>=1),1]
chosing_data<-enhancer_data_log[as.character(tissue_name),as.character(disorder_name)]
# split the data into tissue and cell data #
##cell data

cell_data<-enhancer_data_log[1:6,]
col_name_class<-c("AD_JENSE","Parkinson","ALS","EPILEPSY","MS","TS",
                  "ADHD","ANXIETY","ASD","BIP","SCZ", "OCD",
                  "Neuroticism","RISK_behavior", "INTELLIGENCE","AUD", "INSOMNIA")
cell_data<-cell_data[,col_name_class]
##set the trait #
annote_files_col<-c('Neurological','Neurological','Neurological','Neurological','Neurological','Neurological',
                    'psychiatrics','psychiatrics','psychiatrics','psychiatrics','psychiatrics','psychiatrics',
                    'Trait','Trait','Trait','Trait','Trait')
annote_files_col<-data.frame(Disorder_type=annote_files_col)
rownames(annote_files_col)<-as.character(disorder_name)

annoye_row_cell<-c('Astrocyte','Astrocyte','Meningeal Cells','Neural stem cells','Neurons')
annoye_row_cell<-data.frame(cell_type=annoye_row_cell)
rownames(annoye_row_cell)<-as.character(rownames(cell_data))
## plot the enrichment of cell data # 
#cell_data<-cell_data[,-2]
pheatmap(cell_data,cluster_rows = FALSE,cluster_cols =FALSE,scale='column',border_color = F, 
         fontsize_number = 15,
         color=colorRampPalette(c(brewer.pal(9,"Blues")[8] , "white", brewer.pal(9,"RdBu")[2]))(30),
         display_numbers = matrix(ifelse(cell_data>1.30103,"×",""),nrow(cell_data)),
         annotation_col =annote_files_col,annotation_row = annoye_row_cell)
#set 
graph2ppt(file="cell_enrichment_2.pptx", width=16, height=9)
tissue_data<-chosing_data[c(using_region_new,using_region_adult),disorder_name]
#annote the rownames and colnames # 
get_row_annote_files<-function(tissue_data){
  tissue_name<-rownames(tissue_data)
  adult_data<-vector()
  newborn_data<-vector()
  fetal_data<-vector()
  ###
  annote_files_row<-vector()
  for (i in tissue_name) {
    # i=names_tissue[7]
    if(regexpr("adult",i)[1]>0){
      adult_data<-c(adult_data,i)
      annote_files_row<-c(annote_files_row,'adult')
    }else if(regexpr("newborn",i)[1]>0){
      annote_files_row<-c(annote_files_row,'newborn')
      newborn_data<-c(newborn_data,i)
    } else if(regexpr("fetal",i)[1]>0){
      annote_files_row<-c(annote_files_row,'fetal')
      fetal_data<-c(fetal_data,i)
    }
  }
  annote_files_row<-data.frame(Age_group=annote_files_row)
  rownames(annote_files_row)<-tissue_name
  return(annote_files_row)
}
## col annote #
annote_files_col<-c('Neurological','psychiatrics','Neurological','psychiatrics','psychiatrics','psychiatrics','psychiatrics','Neurological','Trait','Trait','Neurological','psychiatrics','psychiatrics','Neurological',
                    'Trait','psychiatrics','Neurological','Neurological')
annote_files_col<-data.frame(Disorder_type=annote_files_col)
rownames(annote_files_col)<-as.character(disorder_name)


##
pheatmap(tissue_data,cluster_rows = FALSE,cluster_cols = FALSE,scale='column',border_color = F, 
         fontsize_number = 15,
         color=colorRampPalette(c(brewer.pal(9,"Blues")[8] , "white", brewer.pal(9,"RdBu")[2]))(30),
         display_numbers = matrix(ifelse(tissue_data>1.4,"×",""),nrow(tissue_data)),
         annotation_row  =annote_files_row,annotation_col =annote_files_col )
#pheatmap(chosing_data,cluster_rows = TRUE,cluster_cols = TRUE,scale='column',clustering_method = "average",display_numbers=T)
## distinguish the age 

##save imge
ggsave("figure3.tiff", dpi=300)       #这里dpi调整分辨率，可以72，96，300或600，默认是300 
tiff("figure3.tiff", width=1920, height=1080)

tissue_name<-tissue_name[-c()]

###### plot the boxplot ###
get_data<-function(pfdata){
  disorder_name<-colnames(pfdata)
  scale_enhancer_tissue<-scale(pfdata)
  adult_disorder<-scale_enhancer_tissue[intersect(tissue_name,using_region_adult),as.character(disorder_name)]
  newborn_disorder<-scale_enhancer_tissue[intersect(tissue_name,using_region_new),as.character(disorder_name)]
 # fetal_disorder<-scale_enhancer_tissue[intersect(tissue_name,fetal_data),as.character(disorder_name)]
  ###
  adult_disorder_one<-unmatrix(adult_disorder,byrow =T)
  #fetal_disorder_one<-unmatrix(fetal_disorder,byrow =T)
  new_born_disorder_one<-unmatrix(newborn_disorder,byrow =T)
  disorder_len_adult<-rep(as.character(disorder_name),nrow(adult_disorder))
  disorder_len_new_born<-rep(as.character(disorder_name),nrow(newborn_disorder))
 # disorder_len_fetal<-rep(as.character(disorder_name),nrow(fetal_disorder))
  #age_range<-c(rep("adult",length(disorder_len_adult)),rep('new_born',length(disorder_len_new_born)),rep('feral',length(disorder_len_fetal)))
  age_range<-c(rep("adult",length(disorder_len_adult)),rep('new_born',length(disorder_len_new_born)))
  #disorder_enrichment<-c(adult_disorder_one,new_born_disorder_one,fetal_disorder_one)
 # disorder_group<-c(disorder_len_new_born,disorder_len_adult,disorder_len_fetal)
  disorder_enrichment<-c(adult_disorder_one,new_born_disorder_one)
  disorder_group<-c(disorder_len_new_born,disorder_len_adult)
  plot_box_plt<-data.frame(disorder_enrichment=disorder_enrichment,disorder_group=disorder_group,age_range=age_range)
  rownames(plot_box_plt)<-NULL
  ### 
  plot_box_plt$disorder_group <- factor(plot_box_plt$disorder_group,levels = col_name_class_1)
  ggplot(plot_box_plt,aes(x=disorder_group,y=disorder_enrichment))+
    geom_boxplot(aes(fill=age_range))+theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))+
    theme(panel.grid =element_blank())+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
    theme(panel.border = element_blank()) +   ## 删去外层边框
    theme(axis.line = element_line(size=1, colour = "black"))
  #
  graph2ppt(file="cell_enrichment", width=16, height=9)
  
  plot_dataF<-data.frame(region_name=names_tissue,significant_tissue_num=vec_num)
  ggplot(data = plot_dataF,mapping = aes(x = region_name, y = significant_tissue_num)) + geom_bar(stat = 'identity')+theme(axis.text.x = element_text(size = 8, vjust =1, hjust = 1, angle = 45))
  #+theme(plot.title = element_text(hjust = 0.5,size=18)) 
  ## t test on the different group data ## 
  for (i in c(1:ncol(adult_disorder))) {
  #  i=2
    data_test<-wilcox.test(adult_disorder[,i],newborn_disorder[,i])
  #  print(data_test)
    if(data_test$p.value<=0.05){
      print(colnames(adult_disorder)[i])
      print(data_test)
    }
  }
}
get_data(tissue_data[,-2])
