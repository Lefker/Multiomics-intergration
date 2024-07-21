library(TCGAbiolinks)
library(SummarizedExperiment)


# Loads the METH data from the downloaded and saved file.
PRADDnaMethexp <- load("file_path\\PRADDnaMethexp3.rda")
PRADDnaMethexp <- data 
remove(data)

# Delete all rows where the beta values are equal to 0
PRADMeth <- PRADDnaMethexp[rowSums(is.na(assay(PRADDnaMethexp))) == 0,]
PRADMethProbes <- assay(PRADMeth)

# Loads the Gene expression data from the downloaded and saved file.
load(file="file_path\\PRADgenexp.rda")
PRADgenexp<- data 
remove(data)

PRADgene <- as.data.frame(assay(PRADgenexp))

# Loads the mi-RNA data from the downloaded and saved file.
load(file="file_path\\PRADmiRNA_Seq.rda")
PRADmiRNA_Seq <- data
remove(data)
rownames(PRADmiRNA_Seq) <- PRADmiRNA_Seq$miRNA_ID

# using read_count's data 
read_countData <-  colnames(PRADmiRNA_Seq)[grep("count", colnames(PRADmiRNA_Seq))]
PRADmiRNA_Seq <- PRADmiRNA_Seq[,read_countData]
colnames(PRADmiRNA_Seq) <- gsub("read_count_","", colnames(PRADmiRNA_Seq))
rm(read_countData)

# get subtype information
information.subtype <- TCGAquery_subtype(tumor = "PRAD")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-PRAD",type = "clinical") 

# gets the barcodes according to the condition of the samples and orders by row names
PRADgeneids <- get_IDs(PRADgenexp)
PRADgeneids <- PRADgeneids[c(which(grepl("normal",PRADgeneids$condition)),which(grepl("cancer",PRADgeneids$condition))),]
PRADgeneids <- PRADgeneids[order(PRADgeneids$barcode),]

PRADMethids <- get_IDs(PRADDnaMethexp)
PRADMethids <- PRADMethids[c(which(grepl("normal",PRADMethids$condition)),which(grepl("cancer",PRADMethids$condition))),]
PRADMethids  <- PRADMethids [order(PRADMethids$barcode),]

PRADmiRNAids <- get_IDs(PRADmiRNA_Seq)
PRADmiRNAids <- PRADmiRNAids[c(which(grepl("normal",PRADmiRNAids$condition)),which(grepl("cancer",PRADmiRNAids$condition))),]
PRADmiRNAids <- PRADmiRNAids[order(PRADmiRNAids$barcode),]

# Αdds N in Normal samples and C in tumor samples
PRADgeneidsnormal <- PRADgeneids[c(which(grepl("normal",PRADgeneids$condition))),]
PRADgeneidsnormal$patient <- paste("N",sep="_",PRADgeneidsnormal$patient)
PRADgeneidsnormal$participant <- paste("N",sep="_",PRADgeneidsnormal$participant)
PRADgeneidscancer <- PRADgeneids[c(which(grepl("cancer",PRADgeneids$condition))),]
PRADgeneidscancer$patient <- paste("C",sep="_",PRADgeneidscancer$patient)
PRADgeneidscancer$participant <- paste("C",sep="_",PRADgeneidscancer$participant)

PRADMethidsnormal <- PRADMethids[c(which(grepl("normal",PRADMethids$condition))),]
PRADMethidsnormal$patient <- paste("N",sep="_",PRADMethidsnormal$patient)
PRADMethidsnormal$participant <- paste("N",sep="_",PRADMethidsnormal$participant)
PRADMethidscancer <- PRADMethids[c(which(grepl("cancer",PRADMethids$condition))),]
PRADMethidscancer$patient <- paste("C",sep="_",PRADMethidscancer$patient)
PRADMethidscancer$participant <- paste("C",sep="_",PRADMethidscancer$participant)

PRADmiRNAidsnormal <- PRADmiRNAids[c(which(grepl("normal",PRADmiRNAids$condition))),]
PRADmiRNAidsnormal$patient <- paste("N",sep="_",PRADmiRNAidsnormal$patient)
PRADmiRNAidsnormal$participant <- paste("N",sep="_",PRADmiRNAidsnormal$participant)
PRADmiRNAidscancer <- PRADmiRNAids[c(which(grepl("cancer",PRADmiRNAids$condition))),]
PRADmiRNAidscancer$patient <- paste("C",sep="_",PRADmiRNAidscancer$patient)
PRADmiRNAidscancer$participant <- paste("C",sep="_",PRADmiRNAidscancer$participant)

# Finds duplicate samples and removes them
PRADgeneidscancer <- PRADgeneidscancer[!duplicated(PRADgeneidscancer$patient),]
PRADMethidscancer <- PRADMethidscancer[!duplicated(PRADMethidscancer$patient),]
PRADmiRNAidscancer <- PRADmiRNAidscancer[!duplicated(PRADmiRNAidscancer$patient),]

# Find common patient ids
common_patients_normal <- Reduce(intersect, list(PRADgeneidsnormal$patient,PRADMethidsnormal$patient,PRADmiRNAidsnormal$patient))
common_patients_cancer <- Reduce(intersect, list(PRADgeneidscancer$patient,PRADMethidscancer$patient,PRADmiRNAidscancer$patient))

# Keeps only the common patient's sample
PRADgeneidsnormal <- PRADgeneidsnormal[ which(PRADgeneidsnormal$patient %in% common_patients_normal),]

PRADMethidsnormal <- PRADMethidsnormal[ which(PRADMethidsnormal$patient %in% common_patients_normal),]

PRADmiRNAidsnormal <- PRADmiRNAidsnormal[ which(PRADmiRNAidsnormal$patient %in% common_patients_normal),]

PRADgeneidscancer <- PRADgeneidscancer[c(which(PRADgeneidscancer$patient %in% common_patients_cancer)),]

PRADMethidscancer <- PRADMethidscancer[c(which(PRADMethidscancer$patient %in% common_patients_cancer)),]

PRADmiRNAidscancer <- PRADmiRNAidscancer[c(which(PRADmiRNAidscancer$patient %in% common_patients_cancer)),]

# Combine the common patient's ids and orders them by barcode
PRADgeneids <- rbind(PRADgeneidscancer,PRADgeneidsnormal)
PRADgeneids <- PRADgeneids[order(PRADgeneids$barcode),]
PRADMethids <- rbind(PRADMethidscancer,PRADMethidsnormal)
PRADMethids <- PRADMethids[order(PRADMethids$barcode),]
PRADmiRNAids <- rbind(PRADmiRNAidscancer,PRADmiRNAidsnormal)
PRADmiRNAids <- PRADmiRNAids[order(PRADmiRNAids$barcode),]

# From the datasets brings only the omics data of the common patients
PRADgene <- PRADgene[,which(colnames(PRADgene) %in% PRADgeneids$barcode)]
PRADgene <- t(PRADgene)
PRADgene <- PRADgene[order(rownames(PRADgene)),]
rownames(PRADgene) <- PRADgeneids$patient
PRADgene <- t(PRADgene)

PRADMethProbes <- PRADMethProbes[,which(colnames(PRADMethProbes) %in% PRADMethids$barcode)]
PRADMethProbes <- t(PRADMethProbes)
PRADMethProbes <- PRADMethProbes[order(rownames(PRADMethProbes)),]
PRADMethids <- PRADMethids[order(PRADMethids$barcode),]
rownames(PRADMethProbes) <- PRADMethids$patient
PRADMethProbes <- t(PRADMethProbes)

PRADmiRNA_Seq <- PRADmiRNA_Seq[,which(colnames(PRADmiRNA_Seq) %in% PRADmiRNAids$barcode)]
PRADmiRNA_Seq <- t(PRADmiRNA_Seq)
PRADmiRNA_Seq <- PRADmiRNA_Seq[order(rownames(PRADmiRNA_Seq)),]
PRADmiRNAids <- PRADmiRNAids[order(PRADmiRNAids$barcode),]
rownames(PRADmiRNA_Seq) <- PRADmiRNAids$patient
PRADmiRNA_Seq <- t(PRADmiRNA_Seq)

rm(PRADgeneidscancer,PRADgeneidsnormal,PRADMethidscancer,PRADMethidsnormal,
   PRADmiRNAidscancer,PRADmiRNAidsnormal,common_patients_cancer,
   common_patients_normal)

#Saves the ids in files

#setwd("file_path\\")
#write.table(PRADgeneids, "PRADgeneids.tab", sep = "\t", quote = FALSE)

#setwd("file_path\\")
#write.table(PRADMethids, "PRADMethids.tab", sep = "\t", quote = FALSE)

#setwd("file_path\\")
#write.table(PRADmiRNAids, "PRADmiRNAids.tab", sep = "\t", quote = FALSE)

# Saves the data in files without normalization and standarization

#setwd("file_path\\")
#write.table(PRADgene, "PRADgene.tab", sep = "\t", quote = FALSE)

#setwd("file_path\\")
#write.table(PRADMethProbes, "PRADMethProbes.tab", sep = "\t", quote = FALSE)

#setwd("file_path\\")
#write.table(PRADmiRNA_Seq, "PRADmiRNA_Seq.tab", sep = "\t", quote = FALSE)

#export<-PRADall
#setwd("file_path\\")
#write.table(export, "PRADall.tab", sep = "\t",col.names = NA, quote = FALSE)

gc()
