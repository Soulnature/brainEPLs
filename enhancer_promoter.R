# combine the enhancer and promoter #
promoter_dir<-list.files('./promoter/')
enhancer_dir<-list.files('./enhancer_new/')
root_en='./enhancer_new/'
root_pro='./promoter/'
out_dir<-'./enhancer_promoter_end/'
for (i in c(1:length(enhancer_dir))) {
  #i=1
  temp_name=str_split(enhancer_dir[i],'enhancer')%>%unlist()
  data2<-read.table(paste(root_en,enhancer_dir[i],sep = ""),header = FALSE)
  data1<-read.table(paste(root_pro,promoter_dir[i],sep = ""),header = FALSE)
  combind_data<-rbind(data1,data2)
  write.table(combind_data,paste(out_dir,temp_name[1],'enhancer_promoter.bed',sep = ""),row.names = FALSE,col.names = FALSE,quote = FALSE,sep = '\t')
}
