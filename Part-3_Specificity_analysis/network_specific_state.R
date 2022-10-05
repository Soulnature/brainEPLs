# Subnetworks active in only a single group of FANTOM5 tissues for the enhancers and TSSs on chromosome 1.
# over % state
options(stringsAsFactors = F)
library(reshape2)
library(pheatmap)
library(ClassDiscovery)
library(dendextend)
library(gplots)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(Rtsne)
library(ineq)
library(dplyr)
library(Vennerable)

threshold <- c(0.3,0.5,0.7)[2]

meta<- read.table('/Users/songlt/Desktop/fantom/metadata.FANTOM/metadata.FANTOM.biosample.txt', header = T,sep='\t',quote = '')
meta <- unique(meta[,c(1,4,5,6)])
rownames(meta)<- meta$AggregatedBiosampleID

states <- read.table('/Users/songlt/Desktop/fantom/link.aggregated/stats.txt',sep='\t')
states$edge <- paste(states$V1,states$V2,sep='--')
states$class <- meta[states$V4,'c']
states$BiosampleType <- meta[states$V4,'BiosampleType']

num_tissue <- as.data.frame(table(meta$c))
ggplot(num_tissue,aes(x=reorder(Var1,-Freq),y=Freq))+
  geom_bar(stat = 'identity')+
  theme(axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1), legend.position = '')


states_chr1 <- subset(states, grepl('chr1:',V2) & V3>0 & !class%in%c('Other','Cancer') & class%in%as.character(num_tissue$Var1[num_tissue$Freq>=10]))
states_chr1$source <- states_chr1$class


# specific gene
uniq_gene <- unique(states_chr1[,c('V1','V4','source')]) # 

gene_n <- cast(as.data.frame(xtabs(~V1+source,uniq_gene)),formula = V1~source)
rownames(gene_n)<-gene_n$V1
gene_n <- gene_n[,-1]

gene_prop <- t(t(as.data.frame(gene_n))/as.numeric(table(unique(states_chr1[,c('V4','source')])[,2])))

gene_prop[gene_prop>threshold] <- 1
gene_prop[gene_prop<=threshold] <- 0

table(rowSums(gene_prop))

specif_promoter <- subset(as.data.frame( melt(gene_prop)), X1%in%rownames(gene_prop)[rowSums(gene_prop)==1] & value!=0 )[,1:2]
colnames(specif_promoter) <- c('V1','source')
rownames(specif_promoter) <- specif_promoter$V1
table(specif_promoter$source)

# specific enhancer
uniq_enhancer <- unique(states_chr1[,c('V2','V4','source')]) #
enhancer_n <- cast(as.data.frame(xtabs(~V2+source,uniq_enhancer)),formula = V2~source)
rownames(enhancer_n)<-enhancer_n$V2
enhancer_n <- enhancer_n[,-1]

enhancer_prop <- t(t(as.data.frame(enhancer_n))/as.numeric(table(unique(states_chr1[,c('V4','source')])[,2])))

enhancer_prop[enhancer_prop>threshold] <- 1
enhancer_prop[enhancer_prop<=threshold] <- 0

table(rowSums(enhancer_prop))

specif_enhancer <- subset(as.data.frame( melt(enhancer_prop)), X1%in%rownames(enhancer_prop)[rowSums(enhancer_prop)==1] & value!=0 )[,1:2]
colnames(specif_enhancer) <- c('V2','source')
rownames(specif_enhancer) <- specif_enhancer$V2
table(specif_enhancer$source)

#specif edge
states_chr1$edge <-  paste(states_chr1$V1,states_chr1$V2,sep='--')
uniq_edge <- unique(states_chr1[,c('edge','V4','source')]) # 

edge_n <- cast(as.data.frame(xtabs(~edge+source,uniq_edge)),formula = edge~source)
rownames(edge_n)<-edge_n$edge
edge_n <- edge_n[,-1]

edge_prop <- t(t(as.data.frame(edge_n))/as.numeric(table(unique(states_chr1[,c('V4','source')])[,2])))

edge_prop[edge_prop>0.5] <- 1
edge_prop[edge_prop<=0.5] <- 0

table(rowSums(edge_prop))

specif_edge <- subset(as.data.frame( melt(edge_prop)), X1%in%rownames(edge_prop)[rowSums(edge_prop)==1] & value!=0 )[,1:2]
colnames(specif_edge) <- c('edge','source')
rownames(specif_edge) <- specif_edge$edge
table(specif_edge$source)


#cytoscape
# edge
states_chr1 <- unique(states_chr1[,c('V1','V2','edge')])
rownames(states_chr1) <- states_chr1$edge
states_chr1$source <- 'non-specific'
states_chr1[rownames(specif_edge),'source'] <- as.character(specif_edge$source)
states_chr1 <- states_chr1[,c('V1','V2','source')]

virtual_edge <- cbind(rbind(specif_promoter[,c(1,2)]%>%rename(c(source='V2')),specif_enhancer[,c(1,2)]%>%rename(c(V2='V1',source='V2'))),source='virtual_edge')
write.table(subset(rbind(states_chr1,virtual_edge), V1%in%specif_promoter$V1|V2%in%specif_enhancer$V2|rownames(states_chr1)%in%specif_edge$edge), file=paste('~/Desktop/state_edge',threshold,'.txt',sep=''),sep='\t',quote = F, row.names = F)

# ehhancer promoter
node <- data.frame(row.names = c(unique(states_chr1$V1), unique(states_chr1$V2)) )
node$node <- rownames(node)
node$source <-  'non-specific'
node[rownames(specif_enhancer),'source'] <- as.character(specif_enhancer$source)
node[rownames(specif_promoter),'source'] <- as.character(specif_promoter$source)
node$shape<- 'enhancer'
node$shape[grepl('ENSG',node$node)]<- 'promoter'


write.table(rbind(node,cbind(node=levels(virtual_edge$V2),source='virtual',shape='x')), file=paste('~/Desktop/state_node',threshold,'.txt',sep=''),sep='\t',quote = F, row.names = F)


