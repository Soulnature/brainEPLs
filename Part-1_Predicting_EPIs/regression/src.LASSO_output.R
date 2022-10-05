args = commandArgs(trailingOnly=TRUE)

library(glmnet)

#enhancer_file_path = "/Users/yidansun/Dropbox (Personal)/_project_enhancer_promoter/LASSO_010620/LASSO_input.EpiMap_imp.enhancer.1.txt"
#promoter_file_path = "/Users/yidansun/Dropbox (Personal)/_project_enhancer_promoter/LASSO_010620/LASSO_input.EpiMap_imp.promoter.1.txt"

enhancer_file_path = args[1]
promoter_file_path = args[2]

enhancer_matrix <- read.table(enhancer_file_path, header = T, sep="\t", row.names = 1)
promoter_matrix <- read.table(promoter_file_path, header = T, sep="\t", row.names = 1)

row_sub = apply(enhancer_matrix, 1, function(row) !all(row==0))
enhancer_matrix = enhancer_matrix[row_sub, ]

enhancer_matrix = t(enhancer_matrix)
promoter_matrix = t(promoter_matrix)

enhancer_matrix_scale = scale(enhancer_matrix, center=T, scale=T)
promoter_matrix_scale = scale(promoter_matrix, center=T, scale=T)

p = dim(enhancer_matrix_scale)[2]

set.seed(1)
cv_fit = cv.glmnet(enhancer_matrix_scale, promoter_matrix_scale, standardize=F, type.measure="mse", nfolds=100)
lasso_coef = as.matrix(coef(cv_fit, s="lambda.1se")[-1,], p, 1)

#save_file_path = "/Users/yidansun/Dropbox (Personal)/_project_enhancer_promoter/LASSO_010620/LASSO_coef.txt"
save_file_path = args[3]
write.table(lasso_coef,file=save_file_path,sep="\t",row.names=T,col.names=F,quote=FALSE)

