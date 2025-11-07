dat = mtcars
head(dat)

dat |>
  filter(am == 1) |>
  mutate(ratio = hp/wt) |>
  arrange()

install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

gse <- getGEO("GSE33000", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])       # expression matrix
pheno <- pData(gse[[1]])      # sample metadata
head(pheno$`disease status:`) # disease column (if exists)
head(expr)

colnames(pheno) 
table(pheno$patient_id)
table(pheno$`disease status:ch2`)

