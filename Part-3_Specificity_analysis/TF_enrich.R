### TF enrichment
# known
library(viridis)
library(gplots)
library(scales)
library(dplyr)
library(reshape2)
library(stringr)
library(matrixStats)
library(ggsci)
library(ggplot2)

setwd('/Users/songlt/Desktop/fantom/TF_motif/')
know_tfs <- c()
for(cat in dir('./')){
  print(cat)
  know_tf <- read.table(paste0('./',cat,'/knownResults.txt'),sep='\t',quote = '',comment.char = '',header = F)
  know_tf$cat <- cat
  know_tfs <- rbind(know_tfs,know_tf)
}

# heatmap

know_tfs <- know_tfs[know_tfs$V1!='Motif Name',]
know_tfs$TF <- toupper(str_split_fixed(know_tfs$V1,'\\(',2)[,1])
know_tfs$log10p <- -as.numeric((know_tfs$V4))
know_tfs$cat_tf <- paste(know_tfs$cat, know_tfs$TF)
know_tfs <- know_tfs[!duplicated(know_tfs$cat_tf),]


# heatmap matrix
melt_know_tfs <- reshape2::acast((know_tfs[,c('TF','cat','log10p')]),TF~cat)
m_know_tfs <- apply(melt_know_tfs, 1, as.numeric)
rownames(m_know_tfs) <- colnames(melt_know_tfs)
m_know_tfs[m_know_tfs>50] <- 50

melt_know_tfs_FDR <- reshape2::acast((know_tfs[,c('TF','cat','V5')]),TF~cat)
m_know_tfs_FDR <- apply(melt_know_tfs_FDR, 1, as.numeric)
rownames(m_know_tfs_FDR) <- colnames(melt_know_tfs_FDR)

# global: subset TF significant enriched in at least 1 category
sub_know_tfs <- m_know_tfs[,colMins(m_know_tfs_FDR) < 0.05]

rownames(sub_know_tfs) <- car::recode(rownames(sub_know_tfs),"'enhancer_brain_specific'='brain_specific'; 'enhancer_brain'='brain_elevated';
                             'enhancer_brain_adult'='adult_elevated';'enhancer_brain_fetal'='fetal_elevated';'enhancer_brain_newborn'='newborn_elevated'")
rownames(sub_know_tfs)<- factor(rownames(sub_know_tfs), levels = c('brain_specific','brain_elevated','adult_elevated','fetal_elevated','newborn_elevated'))
#sub_know_tfs[sub_know_tfs < -log10(0.05)]   <- NA

dev.off()

heatmap.2(as.matrix(sub_know_tfs)[c(5,1,3,4,2),],density.info = "none", # Remove density legend lines
          trace = "none", # Remove the blue trace lines from heatmap
          Rowv=F,
          dendrogram = "col", # Only plot column dendrogram
          col = c('white',viridis(49)),
          breaks = c(0.3:1.3, 2:50),#breaks = c(0.3:1.3, 1.4:20, seq(21, max(sub_know_tfs) + 10, 10)),
          labCol =NA,cexRow=1,
          na.color=par("bg"),
          margins = c(1, 10),
          #distfun =function(x) dist(x,method =c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1]),
          #hclustfun =function(x) hclust(x,method =c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")[1]),
          key.title = 'log10(P)',
          key.xlab = NA,#'similarity',
          key.ylab = NA, #'Frequency',
          keysize = 1)

dev.print(pdf,file='~/Desktop/fantom/figures/known_tf_heatmap.pdf')

# top:
top_tfs <- know_tfs%>%group_by(cat)%>%top_n(15,(log10p))
top_tfs <- subset(know_tfs, TF%in%top_tfs$TF & as.numeric(V5) < 0.05)

#top_tfs$TF <- factor(top_tfs$TF, levels = c('brain_EleEnhancer_genes','brain_regulated_genes','brain_elevated_genes','adultBrain_EleEnhancer_genes','newbornBrain_EleEnhancer_genes','fetalBrain_EleEnhancer_genes',
 #                                                   'adultBrain_elevated_genes','fetalBrain_elevated_genes','newbornBrain_elevated_genes'))

#top_tfs$cat <- factor(top_tfs$cat, levels = c('enhancer_brain_specific','enhancer_brain','enhancer_brain_adult','enhancer_brain_fetal','enhancer_brain_newborn'))
top_tfs$TF <- factor(top_tfs$TF,levels = unique(top_tfs$TF))

# Expression of TF
promoter_tpm <- read.table('/Users/songlt/Desktop/fantom/promoter_signal.matrix',row.names = 1,header = T,check.names = F)
meta <- read.table('/Users/songlt/Desktop/fantom/metadata.FANTOM/metadata.FANTOM.biosample.txt', header = T,sep='\t',quote = '', row.names = 2)
meta$category <- 'Brain'
meta$category[meta$c!='Brain'] <- 'Non_Brain'
meta$category1 <- meta$category
meta$category1[(grepl('fetal',meta$BiosampleName) | meta$Lifestage%in%c('embryonic') )& meta$category=='Brain'] <- 'fetal'
meta$category1[(grepl('adult',meta$BiosampleName) | meta$Lifestage%in%c('adult') )& meta$category=='Brain'] <- 'adult'
meta$category1[(grepl('newborn',meta$BiosampleName) | meta$Lifestage%in%c('newborn') )& meta$category=='Brain'] <- 'newborn'
meta$category1[meta$category1=='Brain'] <- 'unknown'

gene_nm <- str_split_fixed(rownames(promoter_tpm),'\\|',2)[,2]

tpm_brain <- promoter_tpm[!duplicated(gene_nm),rownames(meta)[meta$category=='Brain']]
tpm_brain <- as.data.frame(rowSums(tpm_brain))
rownames(tpm_brain) <- paste('enhancer_brain',unique(gene_nm))
colnames(tpm_brain) <- 'TPM'

tpm_brain_specific <- promoter_tpm[!duplicated(gene_nm),rownames(meta)[meta$category=='Brain']]
tpm_brain_specific <- as.data.frame(rowSums(tpm_brain))
rownames(tpm_brain_specific) <- paste('enhancer_brain_specific',unique(gene_nm))
colnames(tpm_brain_specific) <- 'TPM'


tpm_adult <- promoter_tpm[!duplicated(gene_nm),rownames(meta)[meta$category1=='adult']]
tpm_adult <- as.data.frame(rowSums(tpm_adult))
rownames(tpm_adult) <- paste('enhancer_brain_adult',unique(gene_nm))
colnames(tpm_adult) <- 'TPM'

tpm_fetal <- promoter_tpm[!duplicated(gene_nm),rownames(meta)[meta$category1=='fetal']]
tpm_fetal <- as.data.frame(rowSums(tpm_fetal))
rownames(tpm_fetal) <- paste('enhancer_brain_fetal',unique(gene_nm))
colnames(tpm_fetal) <- 'TPM'

tpm_newborn <- promoter_tpm[!duplicated(gene_nm),rownames(meta)[meta$category1=='newborn']]
tpm_newborn <- as.data.frame(rowSums(tpm_newborn))
rownames(tpm_newborn) <- paste('enhancer_brain_newborn',unique(gene_nm))
colnames(tpm_newborn) <- 'TPM'

tpm_exp <- rbind(tpm_brain,tpm_brain_specific,tpm_adult,tpm_fetal,tpm_newborn)

top_tfs$TPM <- tpm_exp[top_tfs$cat_tf,'TPM']
top_tfs <- top_tfs[!is.na(top_tfs$TPM),]


top_tfs$cat1 <- car::recode(top_tfs$cat,"'enhancer_brain_specific'='brain_specific'; 'enhancer_brain'='brain_elevated';
                             'enhancer_brain_adult'='adult_elevated';'enhancer_brain_fetal'='fetal_elevated';'enhancer_brain_newborn'='newborn_elevated'")
top_tfs$cat1 <- factor(top_tfs$cat1, levels = rev(c('brain_elevated','brain_specific','fetal_elevated','newborn_elevated','adult_elevated')))


ggplot(top_tfs, aes(x=cat1,y=TF))+
  geom_point(aes(size=log2(TPM),fill=log10p),shape=22)+theme_bw()+
  scale_fill_gradientn(colours=viridis(100), na.value = pal_jco("default")(4)[c(2)],
                       values =rescale(c(-1,-0.5,0,0.5,1)) ,
                       guide="colourbar",name="-log10(P-value)",limits=c(0,50),labels=c(0,10,20,30,40,'>50')) +
  theme(axis.text.x = element_text(angle=90 ,vjust = 0.5,hjust = 1,size=9,face  = 'italic'),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'), 
        legend.key.height = unit (0.3, 'cm'),legend.key.width = unit (0.3, 'cm'),
        legend.title = element_text(size=7), 
        axis.text.y = element_text(size=8),
        text = element_text(size = 9)) +   xlab('')+ylab('')+coord_flip()
  
  
dev.print(pdf,file='~/Desktop/fantom/figures/known_tf.pdf')

### denovo


# TF enrichment
setwd('/Users/songlt/Desktop/fantom/TF_motif/')
denovo_tfs <- c()
for(cat in dir('./')){
  print(cat)
  denovo_tf <- read.table(paste0('./',cat,'/denovo.txt'),sep='\t',quote = '',comment.char = '',header = F)
  denovo_tf$cat <- cat
  denovo_tfs <- rbind(denovo_tfs,denovo_tf)
}

denovo_tfs$TF <- toupper(str_split_fixed(denovo_tfs$V6,'\\(',2)[,1])
denovo_tfs$log10p <- -as.numeric((denovo_tfs$V2))
denovo_tfs$cat_tf <- paste(denovo_tfs$cat, denovo_tfs$TF)
denovo_tfs <- denovo_tfs[!duplicated(denovo_tfs$cat_tf),]

top_tfs <- denovo_tfs#%>%group_by(cat)%>%top_n(15,(log10p))
top_tfs <- subset(know_tfs, TF%in%top_tfs$TF & as.numeric(V5) < 0.05)

#top_tfs$cat <- factor(top_tfs$cat, levels = c('enhancer_brain_specific','enhancer_brain','enhancer_brain_adult','enhancer_brain_fetal','enhancer_brain_newborn'))
top_tfs$TF <- factor(top_tfs$TF,levels = unique(top_tfs$TF))

top_tfs$TPM <- tpm_exp[top_tfs$cat_tf,'TPM']
top_tfs <- top_tfs[!is.na(top_tfs$TPM),]
top_tfs$cat1 <- car::recode(top_tfs$cat,"'enhancer_brain_specific'='brain_specific'; 'enhancer_brain'='brain_elevated';
                             'enhancer_brain_adult'='adult_elevated';'enhancer_brain_fetal'='fetal_elevated';'enhancer_brain_newborn'='newborn_elevated'")
top_tfs$cat1 <- factor(top_tfs$cat1, levels = c('brain_elevated','adult_elevated','fetal_elevated','newborn_elevated','brain_specific'))

ggplot(top_tfs, aes(x=cat1,y=TF))+
  geom_point(aes(size=log2(TPM),fill=log10p),shape=22)+theme_bw()+
  scale_fill_gradientn(colours=viridis(100), na.value = pal_jco("default")(4)[c(2)],
                       values =rescale(c(-1,-0.5,0,0.5,1)) ,
                       guide="colourbar",name="-log10(P-value)",limits=c(0,50),labels=c(0,10,20,30,40,'>50')) +
  theme(axis.text.x = element_text(angle=90 ,vjust = 0.5,hjust = 1,size=9,face  = 'italic'),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'), 
        legend.key.height = unit (0.3, 'cm'),legend.key.width = unit (0.3, 'cm'),
        legend.title = element_text(size=7), 
        axis.text.y = element_text(size=8),
        text = element_text(size = 9)) +   xlab('')+ylab('')+coord_flip()

dev.print(pdf,file='~/Desktop/fantom/figures/denovo_tf.pdf')
