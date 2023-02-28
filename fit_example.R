########################Code to run 10 pairs of count vectors for ECC=1 group########################
library(bzinb, lib.loc = "/home/users/bmlin/my_r_lib")
library(MASS)
library(matrixStats)
vec_list <- read.csv("vec_list.csv") #table of metabolite and species indices
set_in <- vec_list[1:10,]
set_in
#   metab microb
#1      1      1
#2      1      2
#3      1      3
#4      1      4
#5      1      5
#6      1      6
#7      1      7
#8      1      8
#9      1      9
#10     1     10

load("counts_ecc1.RData")
#counts_ecc1 is a list of (1) metabolome counts matrix, (2) microbiome counts matrix for ECC=1 group.
## replace the metabolome counts matrix with microbiome counts matrix for correlation between species.
metab_counts <- counts_ecc1[[1]][,-1]
metab_counts[is.na(metab_counts)] <- 0
microb_counts <- counts_ecc1[[2]][,-1]
corr_out <- NULL

for(line in 1:10){
metab <- set_in[line,1]
microb <- set_in[line,2]
X <- metab_counts[,metab]
Y <- microb_counts[,microb]
X <- round(X/(sd(X[X>0])/30))
Y <- round(Y/(sd(Y[Y>0])/30))


#BZINB correlation
out_list <- bzinb(X,Y,maxiter=10000)$coefficients[,1]
corr_tmp <- (out_list[1]/sqrt((out_list[1]+out_list[2])*(out_list[1]+out_list[3])))*
sqrt((out_list[4]*out_list[5])/((out_list[4]+1)*(out_list[5]+1)))



out_line <- c(metab,microb,corr_tmp)
corr_out <- rbind(corr_out,out_line)

}
save(corr_out,file="res_ecc1/out_ecc1_1.RData")
#example results
#                        a0
#out_line 1  1 0.0778856972
#out_line 1  2 0.0001503718
#out_line 1  3 0.1031816508
#out_line 1  4 0.0008504096
#out_line 1  5 0.0011139067
#out_line 1  6 0.0008793986
#out_line 1  7 0.0313911807
#out_line 1  8 0.0833830678
#out_line 1  9 0.0009769941
#out_line 1 10 0.0024401345

