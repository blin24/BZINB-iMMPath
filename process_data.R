### Filter microbiome data, merge microbiome and metabolome data
filter_file <- read.csv("filter_file.csv")
data_dna <- data$otu[,,1]
data_dna_filtered <- data_dna[rownames(data_dna) %in% filter_file$bacteria[filter_file$filter.D==TRUE],]
data_dna_filtered_t <- t(data_dna_filtered)
data_metabolome_filtered <- data_metabolome[data_metabolome$ID %in% data_out$ID,]

# Separate datasets by ECC disease group. Use microb_t3_b as the binary ECC trait.
temp1 <- data_out[,c("ID","micro_t3_b")]
temp2 <- merge(temp1,data_metabolome_filtered,by.x="ID",by.y="ID")
temp3 <- merge(temp1,data_dna_filtered_t,by.x="ID",by.y="row.names")
sum(temp3$ID==temp2$ID) #check

temp4 <- temp2[temp2$micro_t3_b==0,-2]
temp5 <- temp2[temp2$micro_t3_b==1,-2]
temp6 <- temp3[temp3$micro_t3_b==0,-2]
temp7 <- temp3[temp3$micro_t3_b==1,-2]

sum(temp4$ID==temp6$ID) #check
sum(temp5$ID==temp7$ID) #check
temp8 <- temp2[,-2]
temp9 <- temp3[,-2]

counts_all <- list(metabolome=temp8,microbiome=temp9)
counts_ecc0 <- list(metabolome=temp4,microbiome=temp6)
counts_ecc1 <- list(metabolome=temp5,microbiome=temp7)
counts_list <- list(counts_all=counts_all,counts_ecc0=counts_ecc0,counts_ecc1=counts_ecc1)
saveRDS(counts_list,file="counts_list.RDS")
save(counts_list,file="counts_list.RData",version=2)
