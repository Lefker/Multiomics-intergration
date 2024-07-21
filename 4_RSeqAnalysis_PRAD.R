# Loads the tabs of the expression values and the patient's ids
PRADgene <- read.table("file_path_PRADgene.tab", row.names=1)
PRADgene <- t(PRADgene)
PRADgene <- PRADgene[!apply(PRADgene, 1, function(x) all(x == 0)), ]
PRADgeneids <- read.delim("file_path_PRADgeneids.tab", row.names=1)

# Normalize the expression matrix

Norm_data_Rseq <- norm_fun(PRADgene)

# splits the ids to healthy and disease
healthyids <- PRADgeneids$patient[which(grepl("normal",PRADgeneids$condition))]
cancerids <- PRADgeneids$patient[which(grepl("cancer",PRADgeneids$condition))]

# Splits the dataset to normal and disease
normal <- Norm_data_Rseq[,which(colnames(Norm_data_Rseq) %in% healthyids)]
cancer <- Norm_data_Rseq[,which(colnames(Norm_data_Rseq) %in% cancerids)]

# Combines with the condition
Norm_data_Rseq <- t(Norm_data_Rseq)
Norm_data_Rseq <- Norm_data_Rseq[order(row.names(Norm_data_Rseq)), ]
PRADgeneids <- PRADgeneids[order(PRADgeneids$patient),]
Norm_data_Rseq_Cond <- cbind("Condition"=PRADgeneids$condition, as.data.frame(Norm_data_Rseq))
# Norm_data_Rseq<-as.data.frame(Norm_data_Rseq)
Norm_data_Rseq_Cond <- as.data.frame(Norm_data_Rseq_Cond)

# Gets the p-value from student's t-test
pvalue <- vector()
for( i in 1:nrow(normal)){
  stat.test <- t.test(normal[i,], cancer[i,],var.equal = TRUE)
  pvalue[i] <- stat.test$p.value
}
Norm_data_Rseq <- t(Norm_data_Rseq)
pvalue <- as.matrix(pvalue)
rownames(pvalue)<-rownames(Norm_data_Rseq)
colnames(pvalue)<-c("pvalue")

# Adjusts the p-value with fdr method
set.seed(42)
padj <- p.adjust(pvalue, method = "fdr")
padj <- as.matrix(padj)
rownames(padj) <- rownames(Norm_data_Rseq)
colnames(padj) <- c("p_adjusted")

# Function of the fold change.
# In brackets is needed in the order: normal table, cancer table, combined table
foldchange <- FoldChange_fun(normal,cancer,Norm_data_Rseq)

# Combines the p-adjusted value and the normalized expression data
Comb_data_Rseq <- cbind(padj,foldchange,as.data.frame(Norm_data_Rseq))
Comb_data_Rseq <- Comb_data_Rseq[order(Comb_data_Rseq$FoldChange,decreasing = TRUE),]

# Orders of the table based on the adjusted p-value 
# and keeps the features with adjusted p_value < 0.05
Comb_data_Rseq <- subset(Comb_data_Rseq, p_adjusted < 0.05 & abs(FoldChange) >=1.5, select = c(-1,-2))
Comb_data_Rseq <- as.data.frame(Comb_data_Rseq)

Comb_data_Rseq <- t(Comb_data_Rseq)
Comb_data_Rseq <- Comb_data_Rseq[order(row.names(Comb_data_Rseq),decreasing = FALSE),]
FinalDataCond_Rseq <- cbind("Condition"=PRADgeneids$condition,as.data.frame(Comb_data_Rseq))

# Creates a tab file
setwd("file_path")
write.table(FinalDataCond_Rseq, "significantRSeq.tab", sep = "\t",col.names = NA, quote = FALSE)
