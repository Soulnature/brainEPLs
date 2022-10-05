
# prepare files for visualization of specific network
 
NG_t <- read.table('~/Desktop/fantom/NG_tissue_fantom5.txt')
dataset_meta <-  read.table('~/Desktop/fantom/metadata.FANTOM.dataset.txt',header = T, row.names = 6)
NG_t$BiosampleID <- dataset_meta[NG_t$V1,'BiosampleID']
NG_t <- subset(NG_t,!is.na(BiosampleID))
rownames(NG_t) <- NG_t$BiosampleID

states <- read.table('/Users/songlt/Desktop/fantom/link.single_sample/stats.txt',sep='\t')
states$edge <- paste(states$V1,states$V2,sep='--')
states$class <- meta[states$V4,'c']
states$BiosampleType <- meta[states$V4,'BiosampleType']

states_chr1 <- subset(states, grepl('chr1:',V2) & V3>0 & V4%in%NG_t$BiosampleID )
states_chr1$source <- NG_t[states_chr1$V4,'V2']

unique(states_chr1[,c(4:7)])

table(unique(states_chr1[,c('V4','source')])[,2])
states_chr1 <- unique(states_chr1[,c('V1','V2','edge','source')])

# gene
gene_n <- cast(as.data.frame(xtabs(~V1+source,states_chr1)),formula = V1~source)
rownames(gene_n)<-gene_n$V1
gene_n <- gene_n[,-1]
gene_n[gene_n>1] <- 1
table(rowSums(gene_n))
specif_promoter <- subset( as.data.frame(xtabs(~V1+source,states_chr1)), V1%in%rownames(gene_n)[rowSums(gene_n)==1] & Freq!=0 )
rownames(specif_promoter) <- specif_promoter$V1
table(specif_promoter$source)

# enhancer
enh_n <- cast(as.data.frame(xtabs(~V2+source,states_chr1)),formula = V2~source)
rownames(enh_n)<-enh_n$V2
enh_n <- enh_n[,-1]
enh_n[enh_n>1] <- 1
table(rowSums(enh_n))
specif_enhancer <- subset( as.data.frame(xtabs(~V2+source,states_chr1)), V2%in%rownames(enh_n)[rowSums(enh_n)==1] & Freq!=0 )
rownames(specif_enhancer) <- specif_enhancer$V2
table(specif_enhancer$source)

# edge
edge_n <- cast(as.data.frame(xtabs(~edge+source,states_chr1)),formula = edge~source)
rownames(edge_n)<-edge_n$edge
edge_n <- edge_n[,-1]
table(rowSums(edge_n))
specif_edge <- subset( as.data.frame(xtabs(~edge+source,states_chr1)), edge%in%rownames(edge_n)[rowSums(edge_n)==1] & Freq!=0 )
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
write.table(rbind(states_chr1,virtual_edge), file='~/Desktop/cytoscape_edge_ng.txt',sep='\t',quote = F, row.names = F)

# ehhancer promoter
node <- data.frame(row.names = c(unique(states_chr1$V1), unique(states_chr1$V2)) )
node$node <- rownames(node)
node$source <-  'non-specific'
node[rownames(specif_enhancer),'source'] <- as.character(specif_enhancer$source)
node[rownames(specif_promoter),'source'] <- as.character(specif_promoter$source)
node$shape<- 'enhancer'
node$shape[grepl('ENSG',node$node)]<- 'promoter'

write.table(node, file='~/Desktop/cytoscape_node_ng.txt',sep='\t',quote = F, row.names = F)


# pcc background
pcc <- read.table('~/Desktop/fantom/PCC_output.summary.txt')
pcc_chr1 <- subset(pcc, grepl('chr1',V2))
rownames(pcc_chr1) <- paste(pcc_chr1$V1,pcc_chr1$V2, sep='--')
pcc_chr1$source <- 'non-specific'
pcc_chr1[rownames(specif_edge),'source'] <- as.character(specif_edge$source)
pcc_chr1 <- pcc_chr1[,c('V1','V2','source')]

virtual_edge <- cbind(rbind(specif_promoter[,c(1,2)]%>%rename(c(source='V2')),specif_enhancer[,c(1,2)]%>%rename(c(V2='V1',source='V2'))),source='virtual_edge')
write.table(rbind(pcc_chr1,virtual_edge), file='~/Desktop/cytoscape_edge1_pcc.txt',sep='\t',quote = F, row.names = F)


write.table(pcc_chr1, file='~/Desktop/pcc_edge.txt',sep='\t',quote = F, row.names = F)
