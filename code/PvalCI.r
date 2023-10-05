library("readxl")
library("SNPassoc")
library(dplyr)
library("data.table")

data_path <- "./data/data_Ali.xlsx"

# for saving...
p_vals <- "./data/R_stat_p.csv"
conf_file_path <- './data/confidence.txt'

setwd("~/projects/manhattan_pca/code")
# xls files
my_data <- read_excel(data_path)

# renv.packages('IRkernel')

# renv::install("IRkernel")

# renv::activate()

# my_data

data(my_data, package = "SNPassoc")
# str(my_data, list.len=9)

idx <- grep("^rs", colnames(my_data))
my_data.s <- setupSNP(data=my_data, colSNPs=idx, sep="")

# d <- summary(my_data.s, print=FALSE)

# write.csv(d, "./data/R_stat.csv", row.names=TRUE)

# # plotMissing(my_data.s, print.labels.SNPs = FALSE)

# hwe <- tableHWE(my_data.s)
# # hwe

# hwe2 <- tableHWE(my_data.s, Status)
# # hwe2
# write.csv(hwe2, "./data/hwe.csv", row.names=TRUE)

# snpNHWE <- hwe2[,1]>0.05 & hwe2[,2]<0.05
# rownames(hwe2)[snpNHWE]
# hwe2[snpNHWE,]

# snps.ok <- rownames(hwe2)[hwe2[,1]>=0.01]


# pos <- which(colnames(my_data)%in%snps.ok, useNames = FALSE)
# my_data.s <- setupSNP(my_data, pos, sep="")
# association(Status ~ rs225131, data = my_data.s)

ans1 <- WGassociation(Status, data=my_data.s[, c(1:500)])
ans2 <- WGassociation(Status, data=my_data.s[, c(2, 501:ncol(my_data))])
ans <- rbind(as.data.frame(ans1), as.data.frame(ans2))
ans <- ans[!is.na(ans$codominant),]
ans <- ans[order(ans$codominant), ]
# print(dim(ans))

# plot(ans)

write.csv(ans, p_vals, row.names=TRUE)

ans <- ans[ans$codominant < 0.05, ]
rows = rownames(ans)
snp_list <- rows
results <- list()  # Initialize a list to store results for each SNP

for (snp in snp_list) {
  formula <- as.formula(paste("Status ~", snp))
  result <- association(formula, data = my_data.s)
  results[[snp]] <- result
}

# Capture the printed output
output <- capture.output(print(results))

# Write the captured output to a text file
writeLines(output, con = conf_file_path)


# pos <- which(colnames(my_data.s), useNames = FALSE)
# pos
# pos <- colnames(my_data.s)
# pos
# pos <- which(colnames(my_data.s), useNames = FALSE)
# my_data.s <- setupSNP(my_data, pos, sep="")
# association(Status ~ rs225131, data = my_data.s)
# 
# 
# # Splitting the SNP subset dataframe into smaller subsets
# num_snps <- ncol(my_data.s)
# subset_size <- 2000
# num_subsets <- ceiling(num_snps / subset_size)
# 
# subset_list <- vector("list", num_subsets)
# for (i in 1:num_subsets) {
#   start <- (i - 1) * subset_size + 1
#   end <- min(i * subset_size, num_snps)
#   DATA <- cbind(my_data[["Status"]], my_data.s[, start:end])
#   idx <- grep("^rs", colnames(DATA))
#   DATA <- setupSNP(data=DATA, colSNPs=idx, sep="")
#   subset_list[[i]] <- DATA
# }
# 
# # Running WGassociation on each subset
# result_list <- vector("list", num_subsets)
# for (i in 1:num_subsets) {
#   result_list[[i]] <- WGassociation(subset_list[[i]][, -1], data = subset_list[[i]][, 1])
# }
# 
# # Combining the results
# combined_result <- do.call("rbind", result_list)
# 
# print(combined_result)
# resuols <- WGassociation(Status, my_data.s[, 1:1000])
# 
# # resuols[order(resuols$codominant), ]
# 
# 
# resuols[order(resuols[, 2]), ]
# 
# ans <- WGassociation(Status, data=my_data.s)
# head(ans)

# plot(ans)



# ans

# write.csv(ans, "./data/R_stat_p.csv", row.names=TRUE)

