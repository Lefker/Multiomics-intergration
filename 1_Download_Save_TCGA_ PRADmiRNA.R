                                            ###  miRNA-Seq data  ###
library(TCGAbiolinks)
library(SummarizedExperiment)


# Loads query for the project we need
query <- GDCquery(
  project = "TCGA-PRAD", 
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  experimental.strategy = "miRNA-Seq"
)

# Downloads the files
GDCdownload(
  query = query,
  files.per.chunk = 100
)

# Reads the downloaded data and makes an SummarizedExperiment object
PRADmiRNA_Seq <- GDCprepare(
  query = query, 
  save = TRUE, 
  save.filename = "PRADmiRNA_Seq.rda"
)
