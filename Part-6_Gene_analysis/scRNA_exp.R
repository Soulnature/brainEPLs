
magma <- read.table('~/Desktop/fantom/state-signature genes_id.overlap.zscore.50%.fdr.txt',sep='\t',header = T)
magma <-subset(magma, signature.genes!='' & state_name%in%c('Astrocyte - cerebellum','Astrocyte - cerebral cortex','Neurons',"Oligodendrocyte - precursors"))

colnames(promoter_m) -> ensg
cbind(str_split_fixed(ensg,'\\|',2)[,1],str_split_fixed(ensg,'\\|',2)[,2]) -> gene_info

magma$gene <- ''
for(i in 1:nrow(magma)){
  magma$gene[i] <- paste(gene_info[unlist(strsplit( magma$signature.genes[i],' ')),2],collapse = ' ')
}

write.table(magma, file='~/Desktop/fantom/state-signature_cell.txt',sep='\t',row.names = F)
magma <- read.table('~/Desktop/fantom/state-signature_cell.txt',sep='\t',header = T)


library(Vennerable)
neurodegene_venn <- Venn(list(#AD_Jansen=unlist(strsplit( subset(magma,gwas_file=='AD_Jansen'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
  #ALS=unlist(strsplit( subset(magma,gwas_file=='ALS'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
  MS=unlist(strsplit( subset(magma,gwas_file=='MS'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
  PD=unlist(strsplit( subset(magma,gwas_file=='PD'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
  AD_Kunkle=unlist(strsplit( subset(magma,gwas_file=='AD_Kunkle'&state_name=='Astrocyte - cerebellum')[,'gene'],' '))
))
plot(neurodegene_venn, type = "ChowRuskey", show = list(SetLabels = T),doWeights = T)

neuropsych_venn <- Venn(list(ADHD=unlist(strsplit( subset(magma,gwas_file=='ADHD'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
                             #ASD=unlist(strsplit( subset(magma,gwas_file=='ASD'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
                             BIP=unlist(strsplit( subset(magma,gwas_file=='BIP'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
                             MDD=unlist(strsplit( subset(magma,gwas_file=='MDD'&state_name=='Astrocyte - cerebellum')[,'gene'],' ')),
                             SCZ=unlist(strsplit( subset(magma,gwas_file=='SCZ'&state_name=='Astrocyte - cerebellum')[,'gene'],' '))
))
plot(neuropsych_venn, type = "ChowRuskey", show = list(SetLabels = T),doWeights = T)


load('~/Desktop/fantom/h.cs.RData')

load('~/Desktop/fantom/sub_h.cs.RData')


#Astrocyte_cortex
disease_sub <- c("AD_Jansen" ,   "ADHD"  , "ALS"   ,   "ANXIETY" ,   
                 "ASD"  ,"AUD"  ,    "BIP" ,  "EPILEPSY"  ,"INSOMNIA" ,"INTELLIGENCE" , 
                 "MS"   , "Neuroticism" ,  "OCD" ,   "PD"    ,"RISK_behavior", "SCZ"    ,"TS")

h.cs_Astrocyte_cortex <- subset(h.cs_Astrocyte_cortex, disease%in%disease_sub)

h.cs_Astrocyte_cortex$disease_cat <- car::recode( h.cs_Astrocyte_cortex$disease, "c('AD_Jansen','MS','PD','ALS','AD_Kunkle')='Neurological'; 
                                                     c('MDD','SCZ','BIP','ASD','ADHD')='Psychiatric';else='Trait'")

sub_h.cs_Astrocyte_cortex$disease_cat <- car::recode( sub_h.cs_Astrocyte_cortex$disease, "c('AD_Jansen','MS','PD','ALS','AD_Kunkle')='Neurological'; 
                                                     c('MDD','SCZ','BIP','ASD','ADHD')='Psychiatric';else='Trait'")

h.cs_Astrocyte_cortex$period_cat <- car::recode( h.cs_Astrocyte_cortex$disease, "c('P1','P2','P3','P4','P5','P6','P7')='prenatal';else='postnatal'")
sub_h.cs_Astrocyte_cortex$period_cat <- car::recode( sub_h.cs_Astrocyte_cortex$disease, "c('P1','P2','P3','P4','P5','P6','P7')='prenatal';else='postnatal'")

h.cs_mean_Astrocyte_cortex <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = h.cs_Astrocyte_cortex, mean)
h.cs_median_Astrocyte_cortex <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = h.cs_Astrocyte_cortex, median)
sub_h.cs_mean_Astrocyte_cortex <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = sub_h.cs_Astrocyte_cortex, mean)
sub_h.cs_median_Astrocyte_cortex <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = sub_h.cs_Astrocyte_cortex, median)

#neuron
h.cs_neuron <- subset(h.cs_neuron, disease%in%disease_sub)

h.cs_neuron$disease_cat <- car::recode( h.cs_neuron$disease, "c('AD_Jansen','MS','PD','ALS','AD_Kunkle')='Neurological'; 
                                                     c('MDD','SCZ','BIP','ASD','ADHD')='Psychiatric';else='Trait'")

sub_h.cs_neuron$disease_cat <- car::recode( sub_h.cs_neuron$disease, "c('AD_Jansen','MS','PD','ALS','AD_Kunkle')='Neurological'; 
                                                     c('MDD','SCZ','BIP','ASD','ADHD')='Psychiatric';else='Trait'")

h.cs_neuron$period_cat <- car::recode( h.cs_neuron$disease, "c('P1','P2','P3','P4','P5','P6','P7')='prenatal';else='postnatal'")
sub_h.cs_neuron$period_cat <- car::recode( sub_h.cs_neuron$disease, "c('P1','P2','P3','P4','P5','P6','P7')='prenatal';else='postnatal'")

h.cs_mean_neuron <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = h.cs_neuron, mean)
h.cs_median_neuron <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = h.cs_neuron, median)
sub_h.cs_mean_neuron <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = sub_h.cs_neuron, mean)
sub_h.cs_median_neuron <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = sub_h.cs_neuron, median)


#Astrocyte_OPC
h.cs_OPC <- subset(h.cs_OPC, disease%in%disease_sub)
h.cs_OPC$disease_cat <- car::recode( h.cs_OPC$disease, "c('AD_Jansen','MS','PD','ALS','AD_Kunkle')='Neurological'; 
                                                     c('MDD','SCZ','BIP','ASD','ADHD')='Psychiatric';else='Trait'")

sub_h.cs_OPC$disease_cat <- car::recode( sub_h.cs_OPC$disease, "c('AD_Jansen','MS','PD','ALS','AD_Kunkle')='Neurological'; 
                                                     c('MDD','SCZ','BIP','ASD','ADHD')='Psychiatric';else='Trait'")

h.cs_OPC$period_cat <- car::recode( h.cs_OPC$disease, "c('P1','P2','P3','P4','P5','P6','P7')='prenatal';else='postnatal'")
sub_h.cs_OPC$period_cat <- car::recode( sub_h.cs_OPC$disease, "c('P1','P2','P3','P4','P5','P6','P7')='prenatal';else='postnatal'")


h.cs_mean_OPC <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = h.cs_OPC, mean)
h.cs_median_OPC <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = h.cs_OPC, median)
sub_h.cs_mean_OPC <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = sub_h.cs_OPC, mean)
sub_h.cs_median_OPC <- aggregate(value ~ Gene+disease+Area+Period+ disease_cat+period_cat, data = sub_h.cs_OPC, median)


### BOXPLOT 出生前（P6: 19<= Age < 24 PCW），出生后 40<= Age < 60 Years 全部细胞，均值

ggboxplot(subset(h.cs_mean_neuron,Period%in%c('P6','P14') &Area=='PFC' ), x = "disease", y = "value",
          color = "Period", palette = "jco",size=0.1,
          add = "boxplot")+
  facet_grid(Area~disease_cat,scales = 'free_x')+xlab('')+ylab('Normalized Expression')+
  stat_compare_means(aes(group = Period), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 44,vjust = 0.5,hjust = 0.5))+
  ylim(c(0,1))

dev.print(pdf, file='~/Desktop/fantom/figures/h.cs_mean_neuron.box.pdf')

ggboxplot(subset(h.cs_mean_OPC,Period%in%c('P6','P14') &Area=='PFC' ), x = "disease", y = "value",
          color = "Period", palette = "jco",size=0.1,
          add = "boxplot")+
  facet_grid(Area~disease_cat,scales = 'free_x')+xlab('')+ylab('Normalized Expression')+
  stat_compare_means(aes(group = Period), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 44,vjust = 0.5,hjust = 0.5))+
  ylim(c(0,1))

dev.print(pdf, file='~/Desktop/fantom/figures/h.cs_mean_OPC.box.pdf')


ggboxplot(subset(h.cs_mean_Astrocyte_cortex,Period%in%c('P6','P14') &Area=='PFC' ), x = "disease", y = "value",
          color = "Period", palette = "jco",size=0.1,
          add = "boxplot")+
  facet_grid(Area~disease_cat,scales = 'free_x')+xlab('')+ylab('Normalized Expression')+
  stat_compare_means(aes(group = Period), label = "p.signif")+
  theme(axis.text.x = element_text(angle = 44,vjust = 1,hjust = 1))+
  ylim(c(0,1))

dev.print(pdf, file='~/Desktop/fantom/figures/h.cs_mean_Astrocyte_cortex.box.pdf')


p1 <- ggplot(list(h.cs_mean_neuron,h.cs_median_neuron,
                  sub_h.cs_mean_neuron,sub_h.cs_median_neuron)[[1]],
             aes(x=Period, y=value, group=disease,fill=disease,color=disease)) +xlab('')+ylab('Normalized Expression')+
  geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=T,alpha = 0.15) + 
  facet_wrap(.~disease_cat+Area)+
  xlab('Period') +
  ylab('Expression Value')+
  ggtitle('Neuron')+
  theme_light()+
  theme(legend.key.size = unit(0.5, "cm")) 

p1 <- ggplot(list(h.cs_mean_Astrocyte_cortex,h.cs_median_Astrocyte_cortex,
                  sub_h.cs_mean_Astrocyte_cortex,sub_h.cs_median_Astrocyte_cortex)[[1]],
             aes(x=Period, y=value, group=disease,fill=disease,color=disease)) +xlab('')+ylab('Normalized Expression')+
  geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=F,alpha = 0.15) + 
  facet_wrap(.~disease_cat)+
  xlab('Period') +
  ylab('Expression Value')+
  ggtitle('Astrocyte_cortex')+
  theme_light()+
  theme(legend.key.size = unit(0.5, "cm")) 

ggplot(subset(list(h.cs_mean_Astrocyte_cortex,h.cs_mean_OPC)[[2]], Area=='PFC'),
       aes(x=Period, y=value, group=disease,fill=disease,color=disease)) +
  geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=T,alpha = 0.15) + 
  facet_grid(Area~disease_cat)+
  xlab('Period') +
  ylab('Expression Value')+
  theme_light()+
  theme(legend.key.size = unit(0.5, "cm")) 

dev.print(pdf, file='~/Desktop/fantom/figures/h.cs_mean_neuron_smooth.pdf')
dev.print(pdf, file='~/Desktop/fantom/figures/h.cs_astro_smooth.pdf')
dev.print(pdf, file='~/Desktop/fantom/figures/h.cs_mean_OPC_smooth.pdf')








