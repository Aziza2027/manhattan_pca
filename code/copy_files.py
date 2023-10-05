import shutil


src = "/home/bio/projects/GWAS_analysis-20230810T033758Z-001/GWAS_analysis/variants/all_original.csv"
src2 = "/home/bio/projects/GWAS_analysis-20230810T033758Z-001/GWAS_analysis/variants/removed.csv"

dest = "./code/data/original.csv"
dest2 = "./code/data/filled_dropped.csv"


shutil.copy(src, dest)
shutil.copy(src2, dest2)



