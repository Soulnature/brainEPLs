# to identify brain-elevated and brain period-elevated enhancers 
# remove cancer
# for each stage

library(DESeq2)
library("migest")
library(circlize)
library(VennDiagram)
library(ggplot2)
library(ggpubr)

# enhancer category

enhancer <- (read.table('~/Desktop/fantom/enhancer_count.txt',  header = T,check.names = F))

meta <- read.table('/Users/songlt/Desktop/fantom/metadata.FANTOM/metadata.FANTOM.biosample.txt', header = T,sep='\t',quote = '', row.names = 2)
meta$category <- 'Brain'
meta$category[meta$c!='Brain'] <- 'Non_Brain'
meta$category1 <- meta$category
meta$category1[(grepl('fetal',meta$BiosampleName) | meta$Lifestage%in%c('embryonic') )& meta$category=='Brain'] <- 'fetal'
meta$category1[(grepl('adult',meta$BiosampleName) | meta$Lifestage%in%c('adult') )& meta$category=='Brain'] <- 'adult'
meta$category1[(grepl('newborn',meta$BiosampleName) | meta$Lifestage%in%c('newborn') )& meta$category=='Brain'] <- 'newborn'
meta$category1[meta$category1=='Brain'] <- 'unknown'

#meta$category1[(grepl('fetal|newborn',meta$BiosampleName) | meta$Lifestage%in%c('embryonic','newborn') )& meta$category=='Brain'] <- 'fetal'
#meta$category1[!(grepl('fetal|newborn',meta$BiosampleName) | meta$Lifestage%in%c('embryonic','newborn') ) & meta$category=='Brain'] <- 'adult'

meta$category <- factor(meta$category, levels=c('Non_Brain','Brain'))

# total counts for each group
total_counts_sample <- as.data.frame(colSums(enhancer))
total_counts_sample$category1 <- meta[rownames(total_counts_sample),'category1']
total_counts_sample$category <- meta[rownames(total_counts_sample),'category']
total_counts_sample$c <- meta[rownames(total_counts_sample),'c']


ggplot(subset(total_counts_sample,category1%in%c('adult','fetal','newborn')), 
       aes(x=factor(category1,levels = c('adult','newborn','fetal')),y=`colSums(enhancer)`, color=category1))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')+xlab('')+
  geom_signif(map_signif_level=F,test = "wilcox.test",step_increase = 0.1,
              comparisons = list(c("adult", "fetal"),c("adult", "newborn"),c("newborn", "fetal")))+
  scale_color_jco()

dev.print(pdf, file='~/Desktop/fantom/figures/enhancer_counts_brain.pdf')

ggplot(total_counts_sample, aes(x=reorder(category,`colSums(enhancer)`),y=`colSums(enhancer)`, color=category))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')+xlab('')+
  geom_signif(map_signif_level=F,test = "wilcox.test",step_increase = 0.1,
              comparisons = list(c("Non_Brain","Brain")))+
  scale_color_jco()

dev.print(pdf, file='~/Desktop/fantom/figures/enhancer_counts.pdf')


ggplot(total_counts_sample, aes(x=reorder(c,`colSums(enhancer)`),y=`colSums(enhancer)`, color=category))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')+xlab('')+
  scale_color_jco()


# normalize to TPM
#enhancer <- enhancer[rowSums(enhancer== 0) < ncol(enhancer)*0.7,]
#enhancer <- round(t(t(enhancer)/(colSums(enhancer)/1000000)))
enhancer <- enhancer[rowSums(enhancer!= 0) > 10,]

# count matrix
DEseq2_enhancer <- function(){
  
  meta <- subset(meta,c!='Cancer')
  
  enhancer <- enhancer[,rownames(meta)]
  
  # design= ~ batch + condition
  dds <- DESeqDataSetFromMatrix(countData = enhancer,
                                colData = meta,
                                design = ~ category)
  dds <- estimateSizeFactors( dds )
  
  ## DEseq
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  
  brain_meta <- subset(meta,category1%in%c('fetal','adult'))
  
  brain_meta$category1 <- factor(brain_meta$category1, levels=c('fetal','adult'))
  
  enhancer <- enhancer[,rownames(brain_meta)]
  
  # design= ~ batch + condition
  dds_brain <- DESeqDataSetFromMatrix(countData = enhancer,
                                      colData = brain_meta,
                                      design = ~ category1)
  
  dds_brain <- estimateSizeFactors( dds_brain )
  
  ## DEseq
  dds_brain <- DESeq(dds_brain)
  res_brain <- results(dds_brain)
  res_brain <- as.data.frame(res_brain)
  
  return(list(res=res, res_brain=res_brain))
}

DEseq2_enhancer_period <- function(period){
  
  brain_meta <- subset(meta,category1%in%c('fetal','adult','newborn'))
  
  brain_meta$category2 <- 'control'
  brain_meta$category2[brain_meta$category1==period] <- 'case'
  brain_meta$category2 <- factor(brain_meta$category2, levels=c('control','case'))
  
  enhancer <- enhancer[,rownames(brain_meta)]
  
  # design= ~ batch + condition
  dds_brain <- DESeqDataSetFromMatrix(countData = enhancer,
                                      colData = brain_meta,
                                      design = ~ category2)
  
  dds_brain <- estimateSizeFactors( dds_brain )
  
  ## DEseq
  dds_brain <- DESeq(dds_brain)
  res_brain <- results(dds_brain)
  res_brain <- as.data.frame(res_brain)
  
  res_brain$period <- period
  
  return(res_brain)
}

de_res <- DEseq2_enhancer()
res_fetal <- DEseq2_enhancer_period('fetal')
res_adult <- DEseq2_enhancer_period('adult')
res_newborn <- DEseq2_enhancer_period('newborn')

save(de_res,res_fetal, res_adult,res_newborn, file='~/Desktop/fantom/enhancer.DE.RData')
load('~/Desktop/fantom/enhancer.DE.RData' )
#load('~/Desktop/fantom/enhancer.DE.TPM.RData' )


res_all <- de_res$res
res_all$category <- 'Non_specific'
res_all$category[res_all$log2FoldChange>log2(2) & res_all$padj < 0.05 ] <- 'Brain'
res_all$category[res_all$log2FoldChange< -log2(2) & res_all$padj < 0.05 ] <- 'Non_Brain'
table(res_all$category)

write.table(subset(res_all, category=='Brain'), file='/Users/songlt/Desktop/fantom/enhancer_cat/Brain_elevated.txt',quote = F,row.names = T)
write.table(subset(res_all, category=='Non_Brain'), file='/Users/songlt/Desktop/fantom/enhancer_cat/NonBrain_elevated.txt',quote = F,row.names = T)
write.table(subset(res_all, category=='Non_specific'), file='/Users/songlt/Desktop/fantom/enhancer_cat/Non_specific.txt',quote = F,row.names = T)

res_brain <- de_res$res_brain
res_brain$category1 <- 'Non_specific_brain'
res_brain$category1[res_brain$log2FoldChange>log2(2) & res_brain$padj < 0.05 ] <- 'Adult_brain'
res_brain$category1[res_brain$log2FoldChange< -log2(2) & res_brain$padj < 0.05 ] <- 'Fetal_brain'
table(res_brain$category1)

enhancer_conservation <- read.table('/Users/songlt/Desktop/fantom/enhancer_region_conservation.txt', row.names = 1)

res_all$conservation<- enhancer_conservation[rownames(res_all),'V2']
res_all$category1 <- res_brain$category1
#write.table(rownames(subset(res_all,category=='Brain')), file='~/Desktop/fantom/enhancer_brain_exp.txt',quote = F, sep='\t',row.names = F)

# overlap between brain elevated enhancers and brain specific enhancers
load('~/Desktop/fantom/specif_enhancer_0.5.RData')
#load('~/Desktop/fantom/specif_enhancer_0.3.RData')
intersect(rownames(res_all)[ res_all$category=='Brain'], rownames(specif_enhancer)[specif_enhancer$source=='Brain'])
venn <- venn.diagram(list(expression_elevated=rownames(res_all)[ res_all$category=='Brain'],
                          regulation_specific=rownames(specif_enhancer)[specif_enhancer$source=='Brain']),
                     fill =  pal_jco("default")(2)[c(1,2)], alpha = 0.7, filename = NULL,col='white',
                     cex=1,cat.cex = 1)

plot_grid(venn,scale = 0.5)
dev.print(pdf, file='~/Desktop/fantom/figures/venn_specific.pdf')

# phyper test
enhancer_p <- (read.table('~/Desktop/fantom/enhancer_tpm.txt',  header = T,check.names = F))
bg <- unique(as.character(subset(melt(t(enhancer_p[,rownames(meta)[meta$c=='Brain']])),value>0.1)[,2]))
bg <- rownames(enhancer)
expression_elevated=rownames(res_all)[ res_all$category=='Brain']
regulation_specific=rownames(specif_enhancer)[specif_enhancer$source=='Brain']

venn <- venn.diagram(list(expression_elevated=expression_elevated,
                          regulation_specific=regulation_specific,
                          active_enhancer=bg),
                     fill =  pal_jco("default")(3)[c(1:3)], alpha = 0.7, filename = NULL,col='white',
                     cex=1,cat.cex = 1)

plot_grid(venn,scale = 0.5)
dev.print(pdf, file='~/Desktop/fantom/figures/venn_specific_exp.pdf')

hp <- phyper(length(intersect(expression_elevated,regulation_specific))-1,
             length(regulation_specific),
             length(setdiff(bg,regulation_specific)),
             length(expression_elevated),lower.tail = F)


### brain specific；non-brain specific; non-specific conservation
d0 <- as.data.frame(xtabs(~category+conservation, res_all))
d0$prop <- round(d0$Freq/as.numeric(table(res_all$category)),4)
d1 <- as.data.frame(cbind(category=c(levels(d0$category), levels(d0$conservation)) ,order=1:length(c(levels(d0$category), levels(d0$conservation))), color=pal_jco("default")(length(c(levels(d0$category), levels(d0$conservation)))) ))
d1$num <- c(table(res_all$category),table(res_all$conservation))

brain_conservation_prop <- ggplot(d0, aes(x=category,y=prop, fill=conservation))+
  geom_bar(stat = 'identity',position = 'stack')+theme_bw()+ylab('proportion')+
  scale_fill_manual(values = pal_jco("default")(8)[4:8])+
  theme(text = element_text(size = 14))

fisher.test( matrix(c(d1[1:2,4] - d0[c(13,14),3], d0[c(13,14),3]),nrow=2))

ggsave(brain_conservation_prop, file='~/Desktop/fantom/figures/brain_conservation_prop.pdf',width =4 ,height =4 )

circos.par(start.degree = 90)
chordDiagram(x = d0[,1:3], 
             directional = 2, 
             order = d1$category,
             grid.col = d1$color, 
             annotationTrack = "grid",
             transparency = 0.25,
             annotationTrackHeight = mm_h(c(3, 2)),
             diffHeight  = 0,
             big.gap = 5

)

circos.track(track.index = 1, bg.border = NA, 
             panel.fun = function(x, y) {
               xlim = get.cell.meta.data("xlim")
               sector.index = get.cell.meta.data("sector.index")
               category = d1 %>% filter(category == sector.index) %>% pull(category)
               numx = as.character(d1 %>% filter(category == sector.index) %>% pull(num))
               circos.text(x = mean(xlim), y = 2,labels = category, facing = "bending", cex =0.6)
               circos.text(x = mean(xlim), y = 0.5,labels = numx, facing = "bending", cex =0.6)
             })
circos.clear()
dev.print(pdf, file='~/Desktop/fantom/figures/brain_conservation.pdf')


# conservation of enhancers elevated in different stages 

res_adult$conservation <- enhancer_conservation[rownames(res_adult),]
res_fetal$conservation <- enhancer_conservation[rownames(res_fetal),]
res_newborn$conservation <- enhancer_conservation[rownames(res_newborn),]


res_adult_brain <- subset(res_adult, log2FoldChange>1 & padj < 0.05)
res_fetal_brain <- subset(res_fetal, log2FoldChange>1 & padj < 0.05)
res_newborn_brain <- subset(res_newborn, log2FoldChange>1 & padj < 0.05)

write.table((res_adult_brain), file='/Users/songlt/Desktop/fantom/enhancer_cat/adult_elevated.txt',quote = F,row.names = T)
write.table((res_fetal_brain), file='/Users/songlt/Desktop/fantom/enhancer_cat/fetal_elevated.txt',quote = F,row.names = T)
write.table((res_newborn_brain), file='/Users/songlt/Desktop/fantom/enhancer_cat/newborn_elevated.txt',quote = F,row.names = T)


res_period <- rbind(res_adult_brain,res_fetal_brain,res_newborn_brain)
ggplot(res_period,aes(x=period,fill=conservation))+
  geom_bar()+theme_bw()+ylab('number of enhancers')+
  scale_fill_manual(values=pal_jco("default")(8)[4:8])
dev.print(pdf, file='~/Desktop/fantom/figures/brain_period_n.pdf')

brain_period_p <- ggplot(res_period,aes(x=period,fill=conservation))+
  geom_bar(position = position_fill() )+theme_bw()+ylab('proportion')+
  scale_fill_manual(values=pal_jco("default")(8)[4:8])+theme(text = element_text(size = 14))
#dev.print(pdf, file='~/Desktop/fantom/figures/brain_period_p.pdf')
ggsave(brain_period_p, file='~/Desktop/fantom/figures/brain_period_p.pdf',width =4 ,height =4 )

fisher.test(matrix(c( table(res_period$period) -as.data.frame(xtabs(~period+conservation,res_period))[13:15, 3],
  as.data.frame(xtabs(~period+conservation,res_period))[13:15, 3]),nrow=3)[1:3,])

library(Vennerable)
enhancer_venn <- Venn(list(fetal=rownames(res_fetal_brain),
                          adult=rownames(res_adult_brain),
                          newborn=rownames(res_newborn_brain),
                          brain=rownames(subset(res_all, log2FoldChange>1 & padj < 0.05))))

plot(enhancer_venn, type = "ChowRuskey", show = list(SetLabels = F),doWeights = T)
#plot(enhancer_venn, type = "ChowRuskey", show = list(SetLabels = T),doWeights = T)

plot_grid(venn,scale = 0.5)
dev.print(pdf, file='~/Desktop/fantom/figures/venn_period_enhancer.pdf')


### TF enrichment
write.table(rownames(res_fetal_brain), file='~/Desktop/fantom/enhancer_fetal_brain.bed', quote = F,row.names = F)
write.table(rownames(res_adult_brain), file='~/Desktop/fantom/enhancer_adult_brain.bed', quote = F,row.names = F)
write.table(rownames(res_newborn_brain), file='~/Desktop/fantom/enhancer_newborn_brain.bed', quote = F,row.names = F)
write.table(rownames(subset(res_all, log2FoldChange>1 & padj < 0.05)), file='~/Desktop/fantom/brain_enhancer.bed', quote = F,row.names = F)


## conservation analysis using gnomAD data and phastCon score
## gene_proc: /home1/GENE_proc/SONGLITING/FANTOM/codes/tabix.sh
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
  geom_bar(stat ='identity',position = 'dodge',width = 0.7,alpha=0.7)+
  coord_cartesian(ylim=(c(0.7,1)))+
  theme_bw()+xlab('')+ylab('rare DAF = r/(r+c)')+scale_fill_jco()+theme(legend.position = '',axis.text.x =element_text(angle = 0, hjust = 0.5,vjust = 1))

p2 <- ggplot(subset(freq_rare,cat%in%c('adult','fetal','newborn') ) ,aes(x=reorder(cat,freq1),y=freq1,fill=dataset))+
  geom_bar(stat ='identity',position = 'dodge',width = 0.7,alpha=0.7)+
  coord_cartesian(ylim=(c(0.7,1)))+
  theme_bw()+xlab('')+ylab('rare DAF = r/(r+c)')+scale_fill_jco()+
  theme(axis.text.x =element_text(angle = 0, hjust = 0.5,vjust = 1),legend.title = element_blank(),legend.position = c(0.8,0.9))

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




