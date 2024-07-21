# Loads the tabs of the expression values and the patient's ids
PRADMethProbes <- read.table("file_path_PRADMethProbes.tab", row.names=1)
names(PRADMethProbes) <- gsub(x = names(PRADMethProbes), pattern = "\\.", replacement = "-")
PRADMethProbes <- PRADMethProbes[!apply(PRADMethProbes, 1, function(x) all(x == 0)), ]

PRADMethids <- read.delim("file_path_PRADMethids.tab", row.names=1)

# Normalize the expression matrix
Norm_data_Meth <- norm_fun(PRADMethProbes)

# splits the ids to healthy and disease
healthyids <- PRADMethids$patient[which(grepl("normal",PRADMethids$condition))]
cancerids <- PRADMethids$patient[which(grepl("cancer",PRADMethids$condition))]

# Splits the dataset to normal and disease
normal <- Norm_data_Meth[,which(colnames(Norm_data_Meth) %in% healthyids)]
cancer <- Norm_data_Meth[,which(colnames(Norm_data_Meth) %in% cancerids)]

# Combines with the condition
Norm_data_Meth <- t(Norm_data_Meth)
Norm_data_Meth <- Norm_data_Meth[order(row.names(Norm_data_Meth)), ]
PRADMethids <- PRADMethids[order(PRADMethids$patient),]
Norm_data_Meth_Cond <- cbind("Condition"=PRADMethids$condition, as.data.frame(Norm_data_Meth))
# Norm_data_Meth<-as.data.frame(Norm_data_Meth)
Norm_data_Meth_Cond <- as.data.frame(Norm_data_Meth_Cond)

# Gets the p-value from student's t-test
pvalue <- vector()
for( i in 1:nrow(normal)){
  stat.test <- t.test(normal[i,], cancer[i,],var.equal = TRUE)
  pvalue[i] <- stat.test$p.value
}
Norm_data_Meth <- t(Norm_data_Meth)
pvalue <- as.matrix(pvalue)
rownames(pvalue)<-rownames(Norm_data_Meth)
colnames(pvalue)<-c("pvalue")

# Adjusts the p-value with fdr method
set.seed(42)
padj <- p.adjust(pvalue, method = "fdr")
padj <- as.matrix(padj)
rownames(padj) <- rownames(Norm_data_Meth)
colnames(padj) <- c("p_adjusted")

# Function of the fold change.
# In brackets is needed in the order: normal table, cancer table, combined table
foldchange <- FoldChange_fun(normal,cancer,Norm_data_Meth)

# Combines the p-adjusted value and the normalized expression data
Comb_data_Meth <- cbind(padj, foldchange, as.data.frame(Norm_data_Meth))
Comb_data_Meth <- Comb_data_Meth[order(Comb_data_Meth$FoldChange,decreasing = TRUE),]

# Orders the table based on the adjusted p-value 
# and keeps the features with adjusted p_value < 0.05
Comb_data_Meth <- subset(Comb_data_Meth, p_adjusted < 0.05 & abs(FoldChange) >=1.5, select = c(-1,-2))

Comb_data_Meth <- t(Comb_data_Meth)
Comb_data_Meth <- Comb_data_Meth[order(row.names(Comb_data_Meth),decreasing = FALSE),]
Comb_data_Meth_Cond <- cbind("Condition"=PRADMethids$condition,as.data.frame(Comb_data_Meth))


# Creates a tab file
setwd("file_path")
write.table(Comb_data_Meth_Cond, "significantMethPRAD.tab", sep = "\t",col.names = NA, quote = FALSE)
