library(bzinb, lib.loc = "/home/users/bmlin/my_r_lib")
library(MASS)
library(matrixStats)
set.seed(1 * 1 * 321) #change this based on replicate number


corr_spearman_vec <- NULL
corr_pearson_vec <- NULL
corr_bnb_vec <- NULL
corr_bzinb_vec <- NULL
corr_completedata_spearman_vec <- NULL
corr_completedata_pearson_vec <- NULL

param_temp <- read.table("param_list.txt",header=FALSE)
par <- as.numeric(param_temp[1,])


for(rep in 1:10){

simul_h0 <- rbzinb(300,param = par)
X <- simul_h0[,1]
Y <- simul_h0[,2]
#X <- round(X/(sd(X)/30))
#Y <- round(Y/(sd(Y)/30))


corr_spearman_vec  <- c(corr_spearman_vec,cor(X,Y,method="spearman"))
corr_pearson_vec  <- c(corr_pearson_vec,cor(X,Y))


#calculate complete data correlation
which.nonzero <- (X != 0 & Y != 0)
corr_completedata_spearman_vec <- c(corr_completedata_spearman_vec,cor(X[which.nonzero],Y[which.nonzero],method="spearman"))
corr_completedata_pearson_vec <- c(corr_completedata_pearson_vec,cor(X[which.nonzero],Y[which.nonzero]))


#calculate bnb correlation
out_list <- NA
tryCatch({out_list <- bnb(X,Y,maxiter=1000)$coefficients[,1]},error=function(e){})
if(is.na(out_list)==FALSE){
corr_tmp <- (out_list[1]/sqrt((out_list[1]+out_list[2])*(out_list[1]+out_list[3])))*
sqrt((out_list[4]*out_list[5])/((out_list[4]+1)*(out_list[5]+1)))
}else{
corr_tmp <- NA
}
corr_bnb_vec <- c(corr_bnb_vec,corr_tmp)


#calculate bzinb correlation
out_list <- NA
tryCatch({out_list <- bzinb(X,Y,maxiter=1000)$coefficients[,1]},error=function(e){})
if(is.na(out_list)==FALSE){
corr_tmp <- (out_list[1]/sqrt((out_list[1]+out_list[2])*(out_list[1]+out_list[3])))*
sqrt((out_list[4]*out_list[5])/((out_list[4]+1)*(out_list[5]+1)))
}else{
corr_tmp <- NA
}
corr_bzinb_vec <- c(corr_bzinb_vec,corr_tmp)

}

vec_list <- list(corr_spearman_vec,corr_pearson_vec,corr_bnb_vec,corr_bzinb_vec,corr_completedata_spearman_vec,corr_completedata_pearson_vec)
save(vec_list,file="res/vec_list_1_1.RData")
