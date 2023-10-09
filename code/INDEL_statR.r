library("readxl")
library("SNPassoc")
library(dplyr)
library("data.table")

conf_file_path <- './data/confidence_INDEL_r.txt'
path <- './data/significant_GENOTYPES_encoded.xlsx'

setwd("~/projects/manhattan_pca/code")

my_data <- read_excel(path)

idx <- grep("^rs", colnames(my_data))

my_data.s <- setupSNP(data=my_data, colSNPs=idx, sep="")

for (i in 2:length(colnames(my_data.s))) {
  snp <- colnames(my_data.s)[i]
  formula <- as.formula(paste("Status ~", snp))
  result <- association(formula, data = my_data.s)
  results[[snp]] <- result
}

# Capture the printed output
output <- capture.output(print(results))

# Write the captured output to a text file
writeLines(output, con = conf_file_path)


