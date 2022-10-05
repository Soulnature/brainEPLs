#umap tsne 

library(stringr)
library(ineq)
library(Rtsne)
library(reshape2)
library(umap)
library(cowplot)
library(ggplot2)

setwd('/Users/songlt/Desktop/fantom')

enhancer_p <- (read.table('~/Desktop/fantom/enhancer_count.txt',  header = T,check.names = F))
promoter_p <- (read.table('~/Desktop/fantom/promoter_count.txt',  header = T,check.names = F))

# normalize to TPM
enhancer_m <- ((t(enhancer_p)/(colSums(enhancer_p)/1000000)))
# normalize to TPM
promoter_m <- ((t(promoter_p)/(colSums(promoter_p)/1000000)))

stat_num <- data.frame(row.names = c('promoter','enhancer'))


# meta data

metadata_dataset <- read.table('metadata.FANTOM.dataset.txt',sep='\t',header = T)
metadata_biosample <- read.table('metadata.FANTOM.biosample.txt',sep='\t',header = T,quote = '')
colorcode <- read.table('colorcode.txt',sep='\t',header = T,comment.char = '',stringsAsFactors = F)
rownames(metadata_biosample) <- metadata_biosample$BiosampleID
metadata_biosample <- metadata_biosample[rownames(promoter_m),]

metadata_biosample <- subset(metadata_biosample,BiosampleGroup!='Cancer' )

enhancer_m <- enhancer_m[rownames(metadata_biosample),]
promoter_m <- promoter_m[rownames(metadata_biosample),]

# function
plot_tsne <- function(input_matrix,perplex){
  
  tsne_out <- Rtsne(input_matrix,pca=FALSE,dims=2,perplexity=perplex) # Run TSNE
  tsne_res <- as.data.frame(tsne_out$Y)
  colnames(tsne_res) <- c("tSNE1","tSNE2")
  tsne_res$BiosampleGroup <- metadata_biosample$BiosampleGroup
  
  tsne_plot <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=BiosampleGroup))+
    geom_point()+
    scale_color_manual(breaks = colorcode$BiosampleGroup,
                       values = colorcode$ColorCode) + theme_bw()
  return(tsne_plot)
}

plot_tsne_nolegend <- function(input_matrix,perplex){
  
  tsne_out <- Rtsne(input_matrix,pca=FALSE,dims=2,perplexity=perplex) # Run TSNE
  tsne_res <- as.data.frame(tsne_out$Y)
  colnames(tsne_res) <- c("tSNE1","tSNE2")
  tsne_res$BiosampleGroup <- metadata_biosample$BiosampleGroup
  
  tsne_plot <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=BiosampleGroup))+
    geom_point()+
    scale_color_manual(breaks = colorcode$BiosampleGroup,
                       values = colorcode$ColorCode) + theme_bw()+
    theme(legend.position = '')
  
  return(tsne_plot)
}

plot_umap <- function(input_matrix){
  umap_out <- umap::umap(input_matrix)
  umap_res <- as.data.frame(umap_out$layout)
  colnames(umap_res) <- c("UMAP_1","UMAP_2")
  umap_res$BiosampleGroup <- metadata_biosample$BiosampleGroup
  
  UMAP_plot <- ggplot(umap_res,aes(UMAP_1,UMAP_2,color=BiosampleGroup))+
    geom_point()+
    scale_color_manual(breaks = colorcode$BiosampleGroup,
                       values = colorcode$ColorCode) + theme_bw()
  
  return(UMAP_plot)
}

plot_umap_nolegend <- function(input_matrix){
  umap_out <- umap::umap(input_matrix)
  umap_res <- as.data.frame(umap_out$layout)
  colnames(umap_res) <- c("UMAP_1","UMAP_2")
  umap_res$BiosampleGroup <- metadata_biosample$BiosampleGroup
  
  UMAP_plot <- ggplot(umap_res,aes(UMAP_1,UMAP_2,color=BiosampleGroup))+
    geom_point()+
    scale_color_manual(breaks = colorcode$BiosampleGroup,
                       values = colorcode$ColorCode) + theme_bw()+theme(legend.position = '')
  
  return(UMAP_plot)
}

# 1 all enhancer promoter 

stat_num['enhancer',1] <- ncol(enhancer_m)
stat_num['promoter',1] <- ncol(promoter_m)

## 1.1 tsne

p1 <- plot_tsne_nolegend(enhancer_m,30)
p2 <- plot_tsne(promoter_m,30)

p <- plot_grid(p1,p2,labels = c('enhancer','promoter'),rel_widths = c(1.25,2))
ggsave(p, file='all.tsne_tpm_nocancer.pdf',width=12,height=5)

## 1.2 umap
p1 <- plot_umap_nolegend(enhancer_m)
p2 <- plot_umap(promoter_m)

p <- plot_grid(p1,p2,labels = c('enhancer','promoter'),rel_widths = c(1.25,2))
ggsave(p, file='all.umap.pdf',width=12,height=5)

# 2 activative enhancer promoter 

promoter_usage <- promoter_m>1


stat_num['enhancer',2] <- table(colSums(enhancer_usage)!=0)[2]
stat_num['promoter',2] <- table(colSums(promoter_usage)!=0)[2]

## 2.1 tsne
p1 <- plot_tsne_nolegend(enhancer_m[,colSums(enhancer_usage)!=0 ],50)
p2 <- plot_tsne(promoter_m[,colSums(promoter_usage)!=0 ],50)
p <- plot_grid(p1,p2,labels = c('enhancer','promoter'),rel_widths = c(1.25,2))

ggsave(p, file='active.tsne_tpm_nocancer_50.pdf',width=12,height=5)

## 2.2 umap
p1 <- plot_umap_nolegend(enhancer_m[,colSums(enhancer_usage)!=0 ])
p2 <- plot_umap(promoter_m[,colSums(promoter_usage)!=0 ])
p <- plot_grid(p1,p2,labels = c('enhancer','promoter'),rel_widths = c(1.25,2))

ggsave(p, file='active.umap.pdf',width=12,height=5)


# 3 select enhancer promoter based on gini index 
promoter_gini <- apply(promoter_m,2,ineq)
enhancer_gini <- apply(enhancer_m,2,ineq)

gini_index <- rbind(as.data.frame(promoter_gini)%>%mutate(group='promoter')%>%mutate(gini=promoter_gini)%>%select(-promoter_gini),
                    as.data.frame(enhancer_gini)%>%mutate(group='enhancer')%>%mutate(gini=enhancer_gini)%>%select(-enhancer_gini))

gini_index_hist <- ggplot(gini_index,aes(x = gini,color=group,fill=group))+
  geom_histogram(alpha=0.2)+
  theme_bw()

ggsave(gini_index_hist, file='gini_index_hist.pdf',width = 5,height = 5)

i <- 3
for(gini in c(0.5,0.6,0.7,0.8,0.9)){
  
  stat_num['enhancer',i] <- length(names(enhancer_gini)[enhancer_gini>gini])
  stat_num['promoter',i] <- length(names(promoter_gini)[promoter_gini>gini])
  
  # 3.1 tsne
  p1 <- plot_tsne_nolegend(enhancer_m[, intersect(colnames(enhancer_m),names(enhancer_gini)[enhancer_gini>gini]) ], 50)
  p2 <- plot_tsne(promoter_m[, intersect(colnames(promoter_m),names(promoter_gini)[promoter_gini>gini]) ], 50)
  
  p <- plot_grid(p1,p2,labels = c('enhancer','promoter'),rel_widths = c(1.25,2))
  
  ggsave(p, file=paste('gini', gini,'.tsne.pdf',sep=''), width=12,height=5)
  
  
  ## 3.2 umap
  p1 <- plot_umap_nolegend(enhancer_m[, intersect(colnames(enhancer_m),names(enhancer_gini)[enhancer_gini>gini]) ])
  p2 <- plot_umap(promoter_m[, intersect(colnames(promoter_m),names(promoter_gini)[promoter_gini>gini]) ])
  p <- plot_grid(p1,p2,labels = c('enhancer','promoter'),rel_widths = c(1.25,2))
  
  ggsave(p, file=paste('gini',gini,'.umap.pdf',sep=''), width=12,height=5)
  
  i <- i+1
  
}

colnames(stat_num) <- c('all','active','gini0.5','gini0.6','gini0.7','gini0.8','gini0.9')
write.table(stat_num, file='stat_num.txt',sep='\t',row.names = T,quote = F)

# 
write.table(t(enhancer_m), file='enhancer_tpm.txt',quote = F,sep='\t')
write.table(t(enhancer_usage), file='enhancer_usage.txt',quote = F,sep='\t')
write.table(t(promoter_m), file='promoter_tpm.txt',quote = F,sep='\t')
write.table(t(promoter_usage), file='promoter_usage.txt',quote = F,sep='\t')










