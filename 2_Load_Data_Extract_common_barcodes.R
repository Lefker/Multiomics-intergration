library(SummarizedExperiment)
library(TCGAbiolinks)

# Loads the data if has been saved first in a file.
LUADDnaMethexp <- load(file="file_path_of_meth_data.rda")
LUADDnaMethexp <- data 
remove(data)

# na.omit 
LUADMeth <- LUADDnaMethexp[rowSums(is.na(assay(LUADDnaMethexp))) == 0,]
LUADMethProbes <- assay(LUADMeth)

# Loads the data if has been saved first in a file.
load(file="file_path_of_gene_expr_data.rda")
LUADgenexp<- data 
remove(data)

LUADgene <- as.data.frame(assay(LUADgenexp))

# Loads the data if has been saved first in a file.
load(file="file_path_of_miRNA_seq_data.rda")
LUADmiRNA_Seq <- data
remove(data)
rownames(LUADmiRNA_Seq) <- LUADmiRNA_Seq$miRNA_ID

# using read_count's data 
read_countData <-  colnames(LUADmiRNA_Seq)[grep("count", colnames(LUADmiRNA_Seq))]
LUADmiRNA_Seq <- LUADmiRNA_Seq[,read_countData]
colnames(LUADmiRNA_Seq) <- gsub("read_count_","", colnames(LUADmiRNA_Seq))

# get subtype information
information.subtype <- TCGAquery_subtype(tumor = "LUAD")

# get clinical data
information.clinical <- GDCquery_clinic(project = "TCGA-LUAD",type = "clinical") 

# gets the barcodes and the condition of the samples
LUADgeneids <- get_IDs(LUADgenexp)
LUADgeneids <- LUADgeneids[c(which(grepl("normal",LUADgeneids$condition)),which(grepl("cancer",LUADgeneids$condition))),]
LUADgeneids <- LUADgeneids[order(LUADgeneids$patient),]
LUADMethids <- get_IDs(LUADDnaMethexp)
LUADMethids <- LUADMethids[c(which(grepl("normal",LUADMethids$condition)),which(grepl("cancer",LUADMethids$condition))),]
LUADMethids  <- LUADMethids [order(LUADMethids $patient),]
LUADmiRNAids <- get_IDs(LUADmiRNA_Seq)
LUADmiRNAids <- LUADmiRNAids[c(which(grepl("normal",LUADmiRNAids$condition)),which(grepl("cancer",LUADmiRNAids$condition))),]
LUADmiRNAids <- LUADmiRNAids[order(LUADmiRNAids$patient),]

# Extract the common patients IDs
LUADgeneidsnormal <- LUADgeneids[c(which(grepl("normal",LUADgeneids$condition))),]
LUADgeneidscancer <- LUADgeneids[c(which(grepl("cancer",LUADgeneids$condition))),]
LUADMethidsnormal <- LUADMethids[c(which(grepl("normal",LUADMethids$condition))),]
LUADMethidscancer <- LUADMethids[c(which(grepl("cancer",LUADMethids$condition))),]
LUADmiRNAidsnormal <- LUADmiRNAids[c(which(grepl("normal",LUADmiRNAids$condition))),]
LUADmiRNAidscancer <- LUADmiRNAids[c(which(grepl("cancer",LUADmiRNAids$condition))),]

common_patients_normal <- Reduce(intersect, list(LUADgeneidsnormal$patient,LUADMethidsnormal$patient,LUADmiRNAidsnormal$patient))
common_patients_cancer <- Reduce(intersect, list(LUADgeneidscancer$patient,LUADMethidscancer$patient,LUADmiRNAidscancer$patient))

LUADgeneidsnormal <- LUADgeneidsnormal[ which(LUADgeneidsnormal$patient %in% common_patients_normal),]

LUADMethidsnormal <- LUADMethidsnormal[ which(LUADMethidsnormal$patient %in% common_patients_normal),]

LUADmiRNAidsnormal <- LUADmiRNAidsnormal[ which(LUADmiRNAidsnormal$patient %in% common_patients_normal),]

LUADgeneidscancer <- LUADgeneidscancer[ which(unique(LUADgeneidscancer$patient) %in% common_patients_cancer),]
LUADgeneidscancer <- LUADgeneidscancer[ which(common_patients_cancer %in% LUADgeneidscancer$patient, arr.ind=TRUE),]

LUADMethidscancer <- LUADMethidscancer[ which(unique(LUADMethidscancer$patient) %in% unique(common_patients_cancer)),]
LUADMethidscancer <- LUADMethidscancer[ which(common_patients_cancer %in% LUADMethidscancer$patient, arr.ind=TRUE),]

LUADmiRNAidscancer <- LUADmiRNAidscancer[ which(unique(LUADmiRNAidscancer$patient) %in% common_patients_cancer),]
LUADmiRNAidscancer <- LUADmiRNAidscancer[ which(common_patients_cancer %in% LUADmiRNAidscancer$patient, arr.ind=TRUE),]

LUADgeneids <- rbind(LUADgeneidscancer,LUADgeneidsnormal)
LUADMethids <- rbind(LUADMethidscancer,LUADMethidsnormal)
LUADmiRNAids <- rbind(LUADmiRNAidscancer,LUADmiRNAidsnormal)

a <- vector()
b <- vector()
c <- vector()
for (i in 1:length(common_patients)) {
  a[i] <- which(grepl(common_patients[i], LUADgeneids$patient))
}

for (i in 1:length(common_patients)) {
  b[i] <- which(grepl(common_patients[i], LUADMethids$patient))
}

for (i in 1:length(common_patients)) {
  c[i] <- which(grepl(common_patients[i], LUADmiRNAids$patient))
}

LUADgeneids <- LUADgeneids[ c(a),]
LUADMethids <- LUADMethids[ c(b),]
LUADmiRNAids <- LUADmiRNAids[ c(c),]

a<- setdiff(LUADgeneids$patient,LUADMethids$patient)
d <- setdiff(LUADMethids$patient,LUADgeneids$patient)
b<-setdiff(LUADgeneids$patient,LUADmiRNAids$patient)
c<-setdiff(LUADMethids$patient,LUADmiRNAids$patient)

subasd <- LUADgeneids[(-c(setdiff(LUADgeneids$patient,LUADMethids$patient),setdiff(LUADgeneids$patient,LUADmiRNAids$patient))),]

common_patients <- Reduce(intersect, list(LUADgeneids$patient,LUADMethids$patient,LUADmiRNAids$patient))

LUADgeneids <- LUADgeneids[ which(common_patients %in% LUADgeneids$patient),]

LUADMethids <- LUADMethids[ which(common_patients %in% LUADMethids$patient),]

LUADmiRNAids <- LUADmiRNAids[ which(common_patients %in% LUADmiRNAids$patient),]


#export<-LUADgeneids
# Creates a tab file
#setwd("H:/Ergasia")
#write.table(export, "LUADgeneids.tab", sep = "\t",col.names = NA, quote = FALSE)

#export<-LUADMethids
# Creates a tab file
#setwd("H:/Ergasia")
#write.table(export, "LUADMethids.tab", sep = "\t",col.names = NA, quote = FALSE)

#export<-LUADmiRNAids
# Creates a tab file
#setwd("H:/Ergasia")
#write.table(export, "LUADmiRNAids.tab", sep = "\t",col.names = NA, quote = FALSE)

# Keep the common patients samples
LUADgene <- LUADgene[,which(colnames(LUADgene) %in% LUADgeneids$barcode)]
LUADgene <- LUADgene[,order(colnames(LUADgene))]
colnames(LUADgene) <- LUADgeneids$patient
LUADMethProbes <- LUADMethProbes[,which(colnames(LUADMethProbes) %in% LUADMethids$barcode)]
LUADMethProbes <- LUADMethProbes[,order(colnames(LUADMethProbes))]
colnames(LUADMethProbes) <- LUADMethids$patient
LUADmiRNA_Seq <- LUADmiRNA_Seq[,which(colnames(LUADmiRNA_Seq) %in% LUADmiRNAids$barcode)]
LUADmiRNA_Seq <- LUADmiRNA_Seq[,order(colnames(LUADmiRNA_Seq))]
colnames(LUADmiRNA_Seq) <- LUADmiRNAids$patient

LUADall <- rbind(LUADgene,LUADMethProbes,LUADmiRNA_Seq)

#export<-LUADgene
# Creates a tab file
#setwd("file_path")
#write.table(export, "LUADgene.tab", sep = "\t",col.names = NA, quote = FALSE)

#export<-LUADMethProbes
# Creates a tab file
#setwd("file_path")
#write.table(export, "LUADMethProbes.tab", sep = "\t",col.names = NA, quote = FALSE)

#export<-LUADmiRNA_Seq
# Creates a tab file
#setwd("file_path")
#write.table(export, "LUADmiRNA_Seq.tab", sep = "\t",col.names = NA, quote = FALSE)

# Keep the information from clinical and subtypes of the common patients
inf <- information.subtype[which(information.subtype$patient %in% common_patients),]
infc <- information.clinical[which(information.clinical$submitter_id %in% common_patients),]
normMeth <- LUADMethids[c(which(grepl("normal",LUADMethids$condition))),]
cancMeth <- LUADMethids[c(which(grepl("cancer",LUADMethids$condition))),]
length(which(normMeth$patient %in% inf$patient))/length(which(cancMeth$patient %in% inf$patient))
length(which(grepl("normal",LUADMethids$condition)))/length(which(grepl("cancer",LUADMethids$condition)))

for (i in 1:length(a)) {
  t[i] <- which(LUADgeneids$patient == a[i], arr.ind=TRUE)
}

