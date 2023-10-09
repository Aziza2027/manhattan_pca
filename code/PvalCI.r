library("readxl")
library("SNPassoc")
library(dplyr)
library("data.table")
# library(HardyWeinberg)

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

# Comment next 5 lines if don't wanna use MAF
maf <- summary(my_data.s, print=FALSE)
snps.ok <- rownames(maf)[maf[,2]<=95]
my_data <- cbind(my_data[,1], my_data[,'Status'], my_data[,snps.ok])
idx <- grep("^rs", colnames(my_data))
my_data.s <- setupSNP(data=my_data, colSNPs=idx, sep="")

# my_data.s <- setupSNP(my_data, pos, sep="")

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




