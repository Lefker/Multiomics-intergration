# Loads the tabs of the expression values and the patient's ids
PRADgene <- read.table("file_path_PRADgene.tab", row.names=1)
PRADgene <- t(PRADgene)
PRADgene <- PRADgene[!apply(PRADgene, 1, function(x) all(x == 0)), ]
PRADgeneids <- read.delim("file_path_PRADgeneids.tab", row.names=1)

# Loads the tabs of the beta values and the patient's ids
PRADMethProbes <- read.table("file_path_PRADMethProbes.tab", row.names=1)
names(PRADMethProbes) <- gsub(x = names(PRADMethProbes), pattern = "\\.", replacement = "-")
PRADMethids <- read.delim("file_path_PRADMethids.tab", row.names=1)

# Loads the tabs of the miRNA values and the patient's ids
PRADmiRNA_Seq <- read.delim("file_path_PRADmiRNA_Seq.tab", row.names=1)
names(PRADmiRNA_Seq) <- gsub(x = names(PRADmiRNA_Seq), pattern = "\\.", replacement = "-")
PRADmiRNA_Seq <- PRADmiRNA_Seq[!apply(PRADmiRNA_Seq, 1, function(x) all(x == 0)), ]
PRADrow <- rownames(PRADmiRNA_Seq)
PRADcol <- colnames(PRADmiRNA_Seq)
PRADmiRNAids <- read.delim("file_path_PRADmiRNAids.tab", row.names=1)
PRADmiRNA_Seq <- matrix(unlist(PRADmiRNA_Seq), ncol = nrow(PRADmiRNAids))
rownames(PRADmiRNA_Seq) <- PRADrow
colnames(PRADmiRNA_Seq) <- PRADcol

# normalization between and within columns
Norm_data_Rseq <- norm_fun(PRADgene)
Norm_data_Meth <- norm_fun(PRADMethProbes)
#PRADMethProbes_norm <- as.data.frame(PRADMethProbes_norm)
Norm_data_miRna <- norm_fun(PRADmiRNA_Seq)
#PRADmiRNA_Seq_norm <- as.data.frame(PRADmiRNA_Seq_norm)

# Combines the three normalized datasets
PRADall <- rbind(Norm_data_Rseq,Norm_data_Meth,Norm_data_miRna)

# Saves the integrated dataset
#setwd("file_path")
#write.table(PRADall, "PRADall.tab", sep = "\t", quote = FALSE)
#PRADall <- read.table("file_path_PRADall.tab", row.names=1)

#rm(PRADcol,PRADrow,PRADgene,PRADMethProbes,PRADmiRNA_Seq,PRADgeneids,
#   PRADMethids,PRADall, PRADmiRNA_Seq_norm, PRADMethProbes_norm,PRADgene_norm)
gc()

# splits the ids to healthy and disease
healthyids <- PRADgeneids$patient[which(grepl("normal",PRADgeneids$condition))]
cancerids <- PRADgeneids$patient[which(grepl("cancer",PRADgeneids$condition))]

# Splits the dataset to normal and disease
normal <- PRADall[,which(colnames(PRADall) %in% healthyids)]
cancer <- PRADall[,which(colnames(PRADall) %in% cancerids)]

# Combines with the condition
PRADall <- t(PRADall)
PRADall <- PRADall[order(row.names(PRADall)), ]
PRADgeneids <- PRADgeneids[order(PRADgeneids$patient),]
PRADall_Cond <- cbind("Condition"=PRADgeneids$condition, as.data.frame(PRADall))
# PRADall<-as.data.frame(PRADall)
PRADall_Cond <- as.data.frame(PRADall_Cond)

# Gets the p-value from student's t-test
pvalue <- vector()
for( i in 1:nrow(normal)){
  stat.test <- t.test(normal[i,], cancer[i,],var.equal = TRUE)
  pvalue[i] <- stat.test$p.value
}
PRADall <- t(PRADall)
pvalue <- as.matrix(pvalue)
rownames(pvalue)<-rownames(PRADall)
colnames(pvalue)<-c("pvalue")

# Adjusts the p-value with fdr method
set.seed(42)
padj <- p.adjust(pvalue, method = "fdr")
padj <- as.matrix(padj)
rownames(padj) <- rownames(PRADall)
colnames(padj) <- c("p_adjusted")

# Combines the p-adjusted value and the normalized expression data
# Function of the fold change.
# In brackets is needed in the order: normal table, cancer table, combined table
foldchange <- FoldChange_fun(normal,cancer,PRADall)

# Combines the p-adjusted value and the normalized expression data
Comb_data_all <- cbind(padj,foldchange,as.data.frame(PRADall))
Comb_data_all <- Comb_data_all[order(Comb_data_all$FoldChange,decreasing = TRUE),]

# Orders of the table based on the adjusted p-value 
# and keeps the features with adjusted p_value < 0.05
Comb_data_all <- subset(Comb_data_all, p_adjusted < 0.05 & abs(FoldChange) >=1.5, select = c(-1,-2))
Comb_data_all <- as.data.frame(Comb_data_all)

Comb_data_all <- t(Comb_data_all)
Comb_data_all <- Comb_data_all[order(row.names(Comb_data_all),decreasing = FALSE),]
Comb_data_all_Cond <- cbind("Condition"=PRADgeneids$condition,as.data.frame(Comb_data_all))

# Creates a tab file
setwd("file_path")
write.table(Comb_data_all_Cond, "significantPRADALL.tab", sep = "\t", quote = FALSE)
