library(gProfileR)
library(dplyr)

load('~/Desktop/fantom/specif_promoter_0.5.RData')
load('~/Desktop/fantom/specif_enhancer_0.5.RData')
states <- read.table('/Users/songlt/Desktop/fantom/link.single_sample/stats.txt',sep='\t')
states$edge <- paste(states$V1,states$V2,sep='--')
states$class <- meta[states$V4,'c']
states$BiosampleType <- meta[states$V4,'BiosampleType']

brain_specific_enhancers <- rownames(subset(specif_enhancer, source=='Brain'))# regulating
genes_regulated_by_specific_promoters <- unique(subset(states, class=='Brain' & V2%in%brain_specific_enhancers)[,'V1'])
brain_specific_genes <- rownames(subset(specif_promoter, source=='Brain'))

write.table(brain_specific_enhancers, file='/Users/songlt/Desktop/fantom/enhancer_cat/Brain_specific.txt',quote = F,row.names = F)
write.table(brain_specific_genes, file='/Users/songlt/Desktop/fantom/promoter_cat/Brain_specific.txt',quote = F,row.names = F)


# pathway enrichment analysis for different categories of enhancers
## genes regulated by brain elevated enhancers 
brain_ExpEnhancer_genes <- unique(subset(states,V2%in% rownames(res_all[res_all$category=='Brain',]) & class=='Brain')[,'V1'])
brain_ExpEnhancer_genes_p <- gProfileR::gprofiler(str_split_fixed(brain_ExpEnhancer_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
brain_ExpEnhancer_genes_p$group <- 'brain_EleEnhancer_genes'

## genes regulated by adult brain elevated enhancers
adultBrain_ExpEnhancer_genes <- unique(subset(states,V2%in%rownames(res_adult_brain) & V4%in%rownames(meta[meta$category1=='adult',]))[,'V1'])
adultBrain_ExpEnhancer_genes_p <- gProfileR::gprofiler(str_split_fixed(adultBrain_ExpEnhancer_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
adultBrain_ExpEnhancer_genes_p$group <- 'adultBrain_EleEnhancer_genes'

## genes regulated by fetal brain elevated enhancers
fetalBrain_ExpEnhancer_genes <- unique(subset(states,V2%in% rownames(res_fetal_brain) & V4%in%rownames(meta[meta$category1=='fetal',]))[,'V1'])
fetalBrain_ExpEnhancer_genes_p <- gProfileR::gprofiler(str_split_fixed(fetalBrain_ExpEnhancer_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
fetalBrain_ExpEnhancer_genes_p$group <- 'fetalBrain_EleEnhancer_genes'

## genes regulated by newborn brain elevated enhancers
newbornBrain_ExpEnhancer_genes <- unique(subset(states,V2%in% rownames(res_newborn_brain) & V4%in%rownames(meta[meta$category1=='newborn',]))[,'V1'])
newbornBrain_ExpEnhancer_genes_p <- gProfileR::gprofiler(str_split_fixed(newbornBrain_ExpEnhancer_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
newbornBrain_ExpEnhancer_genes_p$group <- 'newbornBrain_EleEnhancer_genes'

# genes specifically regulated in brain
brain_regulated_genes <- unique(rownames(subset(specif_promoter, source=='Brain')))
brain_regulated_genes_p <- gProfileR::gprofiler(str_split_fixed(brain_regulated_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
brain_regulated_genes_p$group <- 'brain_regulated_genes'

# genes elevated in brain
load('~/Desktop/fantom/promoter_res_all.RData')
brain_exp_genes <- rownames(subset(promoter_res_all, category=='Brain'))
brain_exp_genes_p <- gProfileR::gprofiler(str_split_fixed(brain_exp_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
brain_exp_genes_p$group <- 'brain_elevated_genes'

# genes elevated in adult brain
adultBrain_exp_genes <- rownames(promoter_adult_brain)
adultBrain_exp_genes_p <- gProfileR::gprofiler(str_split_fixed(adultBrain_exp_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
adultBrain_exp_genes_p$group <- 'adultBrain_elevated_genes'

# genes elevated in fetal brain
fetalBrain_exp_genes <- rownames(promoter_fetal_brain)
fetalBrain_exp_genes_p <- gProfileR::gprofiler(str_split_fixed(fetalBrain_exp_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
fetalBrain_exp_genes_p$group <- 'fetalBrain_elevated_genes'

# genes elevated in newborn brain
newbornBrain_exp_genes <- rownames(promoter_newborn_brain)
newbornBrain_exp_genes_p <- gProfileR::gprofiler(str_split_fixed(newbornBrain_exp_genes,'\\|',2 )[,2], src_filter=c('GO:BP'), max_set_size = 500)
newbornBrain_exp_genes_p$group <- 'newbornBrain_elevated_genes'


# top 5 enriched pathway in each category

all_pat <- rbind(brain_ExpEnhancer_genes_p,adultBrain_ExpEnhancer_genes_p,newbornBrain_ExpEnhancer_genes_p,fetalBrain_ExpEnhancer_genes_p,
                 brain_regulated_genes_p,brain_exp_genes_p,adultBrain_exp_genes_p,fetalBrain_exp_genes_p,newbornBrain_exp_genes_p)

top_path <- all_pat%>%group_by(group)%>%top_n(5,(-p.value))
top_path <- subset(all_pat, term.name%in%top_path$term.name)

#top_path$group <- factor(top_path$group, levels = c('brain_EleEnhancer_genes','adultBrain_EleEnhancer_genes','fetalBrain_EleEnhancer_genes','newbornBrain_EleEnhancer_genes',
#                                                    'brain_regulated_genes','brain_elevated_genes','adultBrain_elevated_genes','fetalBrain_elevated_genes','newbornBrain_elevated_genes'))

top_path$group <- factor(top_path$group, levels = c('brain_EleEnhancer_genes','brain_regulated_genes','brain_elevated_genes','adultBrain_EleEnhancer_genes','newbornBrain_EleEnhancer_genes','fetalBrain_EleEnhancer_genes',
                                                    'adultBrain_elevated_genes','fetalBrain_elevated_genes','newbornBrain_elevated_genes'))

top_path$term.name <- factor(top_path$term.name,levels = unique(top_path$term.name))

ggplot(top_path, aes(x=group,y=term.name,fill=-log10(p.value)))+
  geom_tile()+theme_bw()+
  #scale_fill_gradient2(low='#0073C2FF',mid='white',high="#EFC000FF",name='-log10(FDR)') +
  #scale_fill_gradient(low='#0073C2FF',high="#EFC000FF",name='corr') +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 70) )+  
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust =1 ,size=10, lineheight = 0.75),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'), 
        axis.text.y = element_text(size=8),
        text = element_text(size = 8)) +   xlab('')+ylab('')

dev.print(pdf, file='~/Desktop/fantom/figures/brain_path.pdf')

top_pat <- all_pat%>%group_by(group)%>%top_n(5,-p.value)
top_pat$group <- factor(top_pat$group, levels = c('brain_exp_genes','adultBrain_exp_genes',
                                                  'fetalBrain_exp_genes','brain_regulated_genes','brain_ExpEnhancer_genes','adultBrain_ExpEnhancer_genes','fetalBrain_ExpEnhancer_genes'))
ggplot(top_pat, aes(x=term.name,y=-log10(p.value)))+
  geom_bar(stat = 'identity')+facet_grid(.~group, scales='free_x')+
  theme_bw()+  
  theme(axis.text.x =element_text(angle = 90, hjust = 1,vjust = 1), legend.position = '')+xlab('')


dev.print(pdf, file='~/Desktop/fantom/figures/brain_path.pdf')



