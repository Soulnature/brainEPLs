library(reshape2)
library(pheatmap)
library(ClassDiscovery)
library(dendextend)
library(gplots)

meta <- read.table('/Users/songlt/Desktop/fantom/metadata.FANTOM/metadata.FANTOM.biosample.txt', header = T,sep='\t',quote = '')
meta_state <- unique(meta[,c(1,4,5,6)])
meta_state$state_numb <- as.numeric(gsub('State_','',meta_state$AggregatedBiosampleID))
color_code <- read.table('/Users/songlt/Desktop/fantom/metadata.FANTOM/colorcode.txt',sep='\t',header = T,comment.char = '', row.names = 1)
meta_state$color <- color_code[meta_state$c,'ColorCode']

states <- read.table('/Users/songlt/Desktop/fantom/link.aggregated/stats.txt',sep='\t')
states$edge <- paste(states$V1, states$V2, sep='--')
states_m <- dcast(states[,c(4,5,3)],V4~edge)
states_m[is.na(states_m)] <- 0

rownames(states_m) <- states_m$V4
states_m <- states_m[,-1]

# jaccard sim index ab/(|a|^2+|b|^2-ab)
jaccard.sim <- function(ix){
  A <- (X[ix[1],])
  B <- (X[ix[2],])
  return((A %*% B)/(sum(A^2)+sum(B^2)-A %*% B))
}

# cosine sim
cos.sim <- function(ix){
  A = (X[ix[1],])
  B = (X[ix[2],])
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}  

X <- as.matrix(states_m)
n <- nrow(X) 
cmb <- expand.grid(i=1:n, j=1:n) 

cos.sim_matrix <- matrix(apply(cmb,1,cos.sim),n,n)

jaccard.sim_matrix <- matrix(apply(cmb,1,jaccard.sim),n,n)

save(cos.sim_matrix,jaccard.sim_matrix , file='~/Desktop/fantom/sim.RData' )
load('~/Desktop/fantom/sim.RData')

cos.dis_matrix <- 1-cos.sim_matrix
jaccard.dis_matrix  <- 1-jaccard.sim_matrix



# hcluster 
#height <- 0.5
h_method <- c('average','complete')[1]
hd <- hclust(as.dist(1-cos.sim_matrix), members = NULL,method=h_method)
avg_dend_obj <- as.dendrogram(hd)
avg_col_dend <- color_branches(avg_dend_obj, h = height)
plot(avg_col_dend)
#cutree(hd, h = height)

meta_state$cluster  <- cutree(hd, h = height)
meta_state_ht <- meta_state[hd$order,]# 
png(filename="~/Desktop/fantom/figures/cosine.png",width = 4000, height = 4000,res = 300)

heatmap_fig <- heatmap.2(round(cos.sim_matrix,2), trace="none", 
                         scale="none",labRow  =NA,density.info='none',
                         labCol =NA ,revC=T,dendrogram="row",
                         Rowv = avg_dend_obj,
                         Colv = avg_dend_obj,
                         hclust=function(x) hclust(as.dist(1-x),method=h_method),
                         distfun=function(x) (1-x),#
                         ColSideColors=meta_state$color,
                         col=colorRampPalette(c("white","darkblue"))(999),
                         breaks = c(seq(0,0.4,length=150),seq(0.401,0.6,length=350),seq(0.601,1,length=500)),
                         key.title = NA,
                         key.xtickfun=function(){
                           return(list(labels=T,tick=FALSE))},
                         key.ytickfun=function(){
                           return(list(labels=F,tick=FALSE))},
                         key.xlab = NA,#'similarity',
                         key.ylab = NA, #'Frequency',
                         keysize = 0.7,
                         lmat=rbind( c(0,0,5), c(0,1,1 ) ,c(0,4,4),c(3,2,2)), lhei=c(0.1,0.05,0.01,1), lwid=c(1.5,5,1.4 ))#lmat=rbind( c(0,0,4), c(0,3,3 ) ,c(2,1,1)), lhei=c(0.1,0.01,1), lwid=c(1.5,5,1.4 ))
          
dev.off()

dev.print(pdf, file=paste0('~/Desktop/','cosine_',height,h_method,'2.pdf'))




