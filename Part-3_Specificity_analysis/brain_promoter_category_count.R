## Differrntially expressed promoter

library(DESeq2)
library("migest")
library(circlize)
library(VennDiagram)
library(ggsci)
library(ggpubr)

# promoter category

promoter <- (read.table('~/Desktop/fantom/promoter_count.txt',  header = T,check.names = F))

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
total_counts_sample <- as.data.frame(colSums(promoter))
total_counts_sample$category1 <- meta[rownames(total_counts_sample),'category1']
total_counts_sample$category <- meta[rownames(total_counts_sample),'category']
total_counts_sample$c <- meta[rownames(total_counts_sample),'c']

ggplot(subset(total_counts_sample,category1%in%c('adult','fetal','newborn')), 
       aes(x=factor(category1,levels = c('adult','newborn','fetal')),y=`colSums(promoter)`, color=category1))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')+xlab('')+
  geom_signif(map_signif_level=F,test = "wilcox.test",step_increase = 0.1,
              comparisons = list(c("adult", "fetal"),c("adult", "newborn"),c("newborn", "fetal")))+
  scale_color_jco()

dev.print(pdf, file='~/Desktop/fantom/figures/promoter_counts_brain.pdf')

ggplot(total_counts_sample, aes(x=reorder(category,`colSums(promoter)`),y=`colSums(promoter)`, color=category))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')+xlab('')+
  geom_signif(map_signif_level=F,test = "wilcox.test",step_increase = 0.1,
              comparisons = list(c("Non_Brain","Brain")))+
  scale_color_jco()

dev.print(pdf, file='~/Desktop/fantom/figures/promoter_counts.pdf')


ggplot(total_counts_sample, aes(x=reorder(c,`colSums(promoter)`),y=`colSums(promoter)`, color=category))+
  geom_boxplot()+theme_bw()+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')+xlab('')+
  scale_color_jco()

# count matrix
#promoter <- promoter[rowSums(promoter== 0) < ncol(promoter)*0.7,]
promoter <- promoter[rowSums(promoter!= 0) > 10,]

# count matrix

DEseq2_promoter <- function(){
  
  meta <- subset(meta,c!='Cancer')
  
  promoter <- promoter[,rownames(meta)]
  
  # design= ~ batch + condition
  dds <- DESeqDataSetFromMatrix(countData = promoter,
                                colData = meta,
                                design = ~ category)
  dds <- estimateSizeFactors( dds )
  
  ## DEseq
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  
  brain_meta <- subset(meta,category1%in%c('fetal','adult'))
  
  brain_meta$category1 <- factor(brain_meta$category1, levels=c('fetal','adult'))
  
  promoter <- promoter[,rownames(brain_meta)]
  
  # design= ~ batch + condition
  dds_brain <- DESeqDataSetFromMatrix(countData = promoter,
                                      colData = brain_meta,
                                      design = ~ category1)
  
  dds_brain <- estimateSizeFactors( dds_brain )
  
  ## DEseq
  dds_brain <- DESeq(dds_brain)
  res_brain <- results(dds_brain)
  res_brain <- as.data.frame(res_brain)
  
  return(list(res=res, res_brain=res_brain))
}

DEseq2_promoter_period <- function(period){
  
  brain_meta <- subset(meta,category1%in%c('fetal','adult','newborn'))
  
  brain_meta$category2 <- 'control'
  brain_meta$category2[brain_meta$category1==period] <- 'case'
  brain_meta$category2 <- factor(brain_meta$category2, levels=c('control','case'))
  
  promoter <- promoter[,rownames(brain_meta)]
  
  # design= ~ batch + condition
  dds_brain <- DESeqDataSetFromMatrix(countData = promoter,
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

promoter_de_res <- DEseq2_promoter()
promoter_res_fetal <- DEseq2_promoter_period('fetal')
promoter_res_adult <- DEseq2_promoter_period('adult')
promoter_res_newborn <- DEseq2_promoter_period('newborn')


save(promoter_de_res,promoter_res_fetal, promoter_res_adult,promoter_res_newborn, file='~/Desktop/fantom/promoter.DE.RData')

load('~/Desktop/fantom/promoter.DE.RData')

promoter_res_all <- promoter_de_res$res
promoter_res_all$category <- 'Non_specific'
promoter_res_all$category[promoter_res_all$log2FoldChange>log2(2) & promoter_res_all$padj < 0.05 ] <- 'Brain'
promoter_res_all$category[promoter_res_all$log2FoldChange< -log2(2) & promoter_res_all$padj < 0.05 ] <- 'Non_Brain'
table(promoter_res_all$category)

write.table((promoter_res_all[promoter_res_all$category=='Brain',]), file='/Users/songlt/Desktop/fantom/supple_fig_tables/promoter_cat/Brain_elevated.txt',quote = F,row.names = T)
write.table((promoter_res_all[promoter_res_all$category=='Non_Brain',]), file='/Users/songlt/Desktop/fantom/supple_fig_tables/promoter_cat/Non_Brain_elevated.txt',quote = F,row.names = T)


promoter_res_brain <- promoter_de_res$res_brain
promoter_res_brain$category1 <- 'Non_specific_brain'
promoter_res_brain$category1[promoter_res_brain$log2FoldChange>log2(2) & promoter_res_brain$padj < 0.05 ] <- 'Adult_brain'
promoter_res_brain$category1[promoter_res_brain$log2FoldChange< -log2(2) & promoter_res_brain$padj < 0.05 ] <- 'Fetal_brain'
table(promoter_res_brain$category1)

promoter_res_all$category1 <- promoter_res_brain$category1
write.table(rownames(subset(promoter_res_all,category=='Brain')), file='~/Desktop/fantom/promoter_brain_exp.txt',quote = F, sep='\t',row.names = F)
save(promoter_res_all, file='~/Desktop/fantom/promoter_res_all.RData')



promoter_adult_brain <- subset(promoter_res_adult, log2FoldChange>1 & padj < 0.05)
promoter_fetal_brain <- subset(promoter_res_fetal, log2FoldChange>1 & padj < 0.05)
promoter_newborn_brain <- subset(promoter_res_newborn, log2FoldChange>1 & padj < 0.05)

write.table((promoter_adult_brain), file='/Users/songlt/Desktop/fantom/supple_fig_tables/promoter_cat/adult_elevated.txt',quote = F,row.names = T)
write.table((promoter_fetal_brain), file='/Users/songlt/Desktop/fantom/supple_fig_tables/promoter_cat/fetal_elevated.txt',quote = F,row.names = T)
write.table((promoter_newborn_brain), file='/Users/songlt/Desktop/fantom/supple_fig_tables/promoter_cat/newborn_elevated.txt',quote = F,row.names = T)


promoter_period <- rbind(promoter_adult_brain,promoter_fetal_brain,promoter_newborn_brain)
brain_promoter_n <- ggplot(promoter_period,aes(x=period,fill=period))+
  geom_bar()+theme_bw()+ylab('number of promoters')+
  scale_fill_jco()+theme(text = element_text(size = 14))+
  geom_text(aes(label=after_stat(count)), vjust=0, stat = "count")

ggsave(brain_promoter_n, file='~/Desktop/fantom/figures/brain_promoter_n.pdf',width =4 ,height =2.5 )
  
dev.print(pdf, file='~/Desktop/fantom/figures/brain_promoter_n.pdf')

brain_period_promoter_n <- ggplot(promoter_res_all,aes(x=category,fill=category))+
  geom_bar()+theme_bw()+ylab('number of promoters')+
  scale_fill_jco()+theme(axis.text.x  = element_text(size = 8))+
  geom_text(aes(label=after_stat(count)), vjust=0, stat = "count")

ggsave(brain_period_promoter_n, file='~/Desktop/fantom/figures/brain_period_promoter_n.pdf',width =4 ,height =2.3 )


# oevrlap

library(Vennerable)
promoter_venn <- Venn(list(fetal=rownames(promoter_fetal_brain),
                 adult=rownames(promoter_adult_brain),
                 newborn=rownames(promoter_newborn_brain),
                 brain=rownames(subset(promoter_res_all, log2FoldChange>1 & padj < 0.05))))

plot(promoter_venn, type = "ChowRuskey", show = list(SetLabels = T),doWeights = T)
plot(promoter_venn, type = "ChowRuskey", show = list(SetLabels = F),doWeights = T)

dev.print(pdf, file='~/Desktop/fantom/figures/venn_period_promoter.pdf')

# 

load('~/Desktop/fantom/specif_promoter_0.5.RData')
brain_specific_promoters <- rownames(subset(specif_promoter, source=='Brain'))# regulating

# the genes specifically expressed in brain group
brain_exp_genes <- rownames(subset(promoter_res_all, category=='Brain'))

# specific promoters
load('~/Desktop/fantom/specif_promoter_0.5.RData')
load('~/Desktop/fantom/specif_enhancer_0.5.RData')
states <- read.table('/Users/songlt/Desktop/fantom/link.single_sample/stats.txt',sep='\t')
states$edge <- paste(states$V1,states$V2,sep='--')
states$class <- meta[states$V4,'c']
states$BiosampleType <- meta[states$V4,'BiosampleType']

brain_specific_enhancers <- rownames(subset(specif_enhancer, source=='Brain'))# regulating
genes_regulated_by_specific_enhancer <- unique(subset(states, class=='Brain' & V2%in%brain_specific_enhancers)[,'V1'])
brain_specific_genes <- rownames(subset(specif_promoter, source=='Brain'))

intersect(brain_exp_genes, genes_regulated_by_specific_enhancer)

venn <- venn.diagram(list(elevated_genes=brain_exp_genes,
                          genes_regulated_by_specific_enhancers=genes_regulated_by_specific_enhancer,
                          specifically_regulated_genes=brain_specific_genes),
                     fill =  pal_jco("default")(3), alpha = 0.7, filename = NULL,col='white',
                     cex=1,cat.cex = 1)

plot_grid(venn,scale = 0.5)
dev.print(pdf, file='~/Desktop/fantom/figures/venn_specific_exp.pdf')

