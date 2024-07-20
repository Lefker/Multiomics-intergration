library(TCGAbiolinks)
library(SummarizedExperiment)


# Loads query for the project we need
query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450"
)

# Downloads the files
GDCdownload(
  query = query,
  files.per.chunk = 100
)

# Reads the downloaded data and makes an R object
PRADDnaMethexp <- GDCprepare(
  query = query, 
  save = TRUE,
  save.filename = "PRADDnaMethexp3.rda"
)