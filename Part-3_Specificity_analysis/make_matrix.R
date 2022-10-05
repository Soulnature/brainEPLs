

library(stringr)
library(ineq)
library(Rtsne)
library(reshape2)
library(umap)

setwd('/Users/songlt/Desktop/fantom')



#********************* meta data *********************
metadata_dataset <- read.table('metadata.FANTOM.dataset.txt',sep='\t',header = T)
metadata_biosample <- read.table('metadata.FANTOM.biosample.txt',sep='\t',header = T,quote = '')
colorcode <- read.table('colorcode.txt',sep='\t',header = T,comment.char = '',stringsAsFactors = F)
rownames(metadata_biosample) <- metadata_biosample$BiosampleID



#********************* promoter data *********************
metadata_dataset <- read.table('metadata.FANTOM.dataset.txt',sep='\t',header = T)
metadata_biosample <- read.table('metadata.FANTOM.biosample.txt',sep='\t',header = T,quote = '')
colorcode <- read.table('colorcode.txt',sep='\t',header = T,comment.char = '',stringsAsFactors = F)
rownames(metadata_biosample) <- metadata_biosample$BiosampleID

promoter <- read.table('promoter_all.txt')

colnames(promoter) <- c('encode_acc','depth','sample')

promoter_m <- acast(promoter,sample~encode_acc,value.var=c('depth'))

metadata_biosample <- metadata_biosample[rownames(promoter_m),]


#********************* enhancer data *********************
enhancer <- read.table('enhancer_all.txt')

enhancer$chr <- str_split_fixed(enhancer$V1,':',2)[,1]

colnames(enhancer) <- c('encode_acc','depth','sample','chr')

enhancer <- subset(enhancer, !chr%in%c("chrX",  "chrY" , "chrM"))  
enhancer_m <- acast(enhancer,sample~encode_acc,value.var=c('depth'))
#enhancer_m[is.na(enhancer_m)] <- 0

#********************* enhancer usage *********************
enhancer_u <- read.table('enhancer_usage.txt')

enhancer_u$chr <- str_split_fixed(enhancer_u$V1,':',2)[,1]

colnames(enhancer_u) <- c('encode_acc','depth','sample','chr')

enhancer_u <- subset(enhancer_u, !chr%in%c("chrX",  "chrY" , "chrM"))  
enhancer_usage <- acast(enhancer_u,sample~encode_acc,value.var=c('depth'))
enhancer_usage <- enhancer_usage==1

save.image('/Users/songlt/Desktop/fantom/fantom_matrix.RData')


