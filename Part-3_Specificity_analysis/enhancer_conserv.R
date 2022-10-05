setwd('/Users/songlt/Desktop/fantom/enhancer_cat')

freq_rare <- c()
for(enhan_cat in dir('/Users/songlt/Desktop/fantom/enhancer_cat',pattern = '_gnomad_daf.txt')){
  #enhan_cat <- sub('.gnomad.txt','',enhan_cat)
  print(enhan_cat)
  enha <- read.table(paste0(enhan_cat))
  enha <- subset(enha, !V1%in%c('chrX','chrY'))
  enha$pos <- paste(enha$V1,enha$V2)
  #enha$maf <- as.numeric(gsub(pattern ='A.*=' ,replacement ='', enha$V10))
  enha$maf <-  as.numeric(enha$V4)
  #enha1  <- enha%>%group_by(pos)%>%top_n(1,maf)
  enha1 <- enha
  #n_rare <- table(enha1$maf>0.001 & enha1$maf<0.01)[2]
  #n_nonrare <- table(enha1$maf>0.01)[2]
  
  n_rare <- table(enha1$maf<0.0001)[2]
  n_nonrare <- table(enha1$maf>=0.0001)[2]
  
  
  
  # enha_region <- read.table(paste0(enhan_cat,'.sort.bed'))
  # enha_region <- subset(enha_region,!V1%in%c('chrX','chrY'))
  # print(enhan_cat)
  # 
  # enha_region$length <- enha_region$V3-enha_region$V2+1
  # n_allele <- sum(enha_region$length)
  
  freq_rare <- rbind(freq_rare, c(n_rare,n_nonrare))
}

freq_rare <- as.data.frame(freq_rare)
colnames(freq_rare) <- c('n_rare','n_nonrare')

freq_rare$cat <-  str_split_fixed(dir('/Users/songlt/Desktop/fantom/enhancer_cat',pattern = '_gnomad_daf.txt'),'_',2)[,1]
freq_rare$cat[freq_rare$cat=='Non'] <- 'Non_specific'

freq_rare$data <-  str_split_fixed(dir('/Users/songlt/Desktop/fantom/enhancer_cat',pattern = '_gnomad_daf.txt'),'\\.',3)[,2]
#freq_rare$dataset <- '1000G'
freq_rare$dataset[freq_rare$data=='_gnomad_daf'] <- 'gnomAD_DAF'


freq_rare$freq1 <- freq_rare$n_rare/(freq_rare$n_rare+freq_rare$n_nonrare)
#freq_rare$freq2 <- 1-freq_rare$n_nonrare/freq_rare$n_allele


p1 <- ggplot(subset(freq_rare,!cat%in%c('adult','fetal','newborn') ) ,aes(x=reorder(cat,freq1),y=freq1,fill=dataset))+
  geom_bar(stat ='identity',position = 'dodge',width = 0.7)+
  coord_cartesian(ylim=(c(0.7,0.75)))+
  theme_bw()+xlab('')+ylab('rare DAF = r/(r+c)')+scale_fill_jco()+theme(legend.position = '',axis.text.x =element_text(angle = 0, hjust = 0.5,vjust = 1))

p2 <- ggplot(subset(freq_rare,cat%in%c('adult','fetal','newborn') ) ,aes(x=reorder(cat,freq1),y=freq1,fill=dataset))+
  geom_bar(stat ='identity',position = 'dodge',width = 0.7)+
  coord_cartesian(ylim=(c(0.7,0.75)))+
  theme_bw()+xlab('')+ylab('rare DAF = r/(r+c)')+scale_fill_jco()+
  theme(axis.text.x =element_text(angle = 0, hjust = 0.5,vjust = 1),legend.title = element_blank(),legend.position = '')#c(0.8,0.9))

plot_grid(p1,p2)
dev.print(pdf, file='~/Desktop/fantom/figures/gnomad_conservation.pdf')



enha_phastCon <- c()
for(enhan_cat in dir('/Users/songlt/Desktop/fantom/enhancer_cat',pattern = 'phastCons.bed')){
  #enhan_cat <- sub('.gnomad.txt','',enhan_cat)
  print(enhan_cat)
  enha <- read.table(paste0(enhan_cat))
  enha$cat <- str_split_fixed(enhan_cat,'_',2)[,1]
  enha$cat[enha$cat=='Non'] <- 'Non_specific'
  enha_phastCon <- rbind(enha_phastCon,enha)
}

library(ggplot2)
library(ggsci)
library(ggpubr)

enha_phastCon$cat <- factor(enha_phastCon$cat, levels = c("adult", "fetal",'newborn',"Non_specific", "NonBrain","Brain" ))

p2_phastCon <- ggplot(subset(enha_phastCon,cat%in%c('adult','fetal','newborn') ) ,aes(x=cat,y=V6,fill=cat,color=cat))+
  geom_boxplot(width=0.05,outlier.size =0.05,alpha=0.2)+
  geom_violin(alpha=0.2,width=0.5)+
  scale_fill_jco()+theme_bw()+#coord_cartesian(ylim=(c(0.7,1)))+
  theme_bw()+xlab('')+ylab('PhastCons score')+scale_fill_jco()+scale_color_jco()+
  theme(axis.text.x =element_text(angle = 0, hjust = 0.5,vjust = 1),legend.title = element_blank(),legend.position = '')+
  geom_signif(comparisons = list(c("adult", "fetal"),
                                 c("adult", "newborn"),
                                 c("fetal", "newborn")),
              step_increase=0.05)

p1_phastCon <- ggplot(subset(enha_phastCon,!cat%in%c('adult','fetal','newborn') ) ,aes(x=reorder(cat,V6),y=V6,fill=cat,color=cat))+
  geom_boxplot(width=0.05,outlier.size =0.05,alpha=0.2)+
  geom_violin(alpha=0.2,width=0.5)+
  theme_bw()+#coord_cartesian(ylim=(c(0.7,1)))+
  theme_bw()+xlab('')+ylab('PhastCons score')+scale_fill_jco()+scale_color_jco()+
  theme(axis.text.x =element_text(angle = 0, hjust = 0.5,vjust = 1),legend.position = '')+
  geom_signif(comparisons = list(c("Non_specific", "NonBrain"),
                                 c("Non_specific", "Brain"),
                                 c("NonBrain", "Brain")),
              step_increase=0.05)

plot_grid(p1_phastCon,p2_phastCon)



