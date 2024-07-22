                                                    ###  Gene expression data  ###
library(TCGAbiolinks)
library(SummarizedExperiment)


# Loads query for the project we need
query.PRAD.TranscrProf.genexp <- GDCquery(
  project = "TCGA-PRAD", 
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# Downloads the files
GDCdownload(
  query = query.PRAD.TranscrProf.genexp,
  files.per.chunk = 100
)

# Reads the downloaded data and makes an SummarizedExperiment object
PRADgenexp <- GDCprepare(
  query = query.PRAD.TranscrProf.genexp, 
  save = TRUE, 
  save.filename = "PRADgenexp.rda"
)
