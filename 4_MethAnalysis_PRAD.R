# Loads the tabs of the expression values and the patient's ids
PRADMethProbes <- read.table("K:\\METAPTYXIAKO\\Diplomatikh\\Ergasia\\PRADMethProbes.tab", row.names=1)
names(PRADMethProbes) <- gsub(x = names(PRADMethProbes), pattern = "\\.", replacement = "-")
PRADMethProbes <- PRADMethProbes[!apply(PRADMethProbes, 1, function(x) all(x == 0)), ]

PRADMethids <- read.delim("K:\\METAPTYXIAKO\\Diplomatikh\\Ergasia\\PRADMethids.tab", row.names=1)

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

# Function of the fold change.
# In brackets is needed in the order: normal table, cancer table, combined table
foldchange <- FoldChange_fun(normal,cancer,Norm_data_Meth)

# Combines the p-adjusted value and the normalized expression data
Comb_data_Meth <- cbind(padj, foldchange, as.data.frame(Norm_data_Meth))
Comb_data_Meth <- Comb_data_Meth[order(Comb_data_Meth$FoldChange,decreasing = TRUE),]

# Orders of the table based on the adjusted p-value 
# and keeps the features with adjusted p_value < 0.05
Comb_data_Meth <- subset(Comb_data_Meth, p_adjusted < 0.05 & abs(FoldChange) >=1.5, select = c(-1,-2))

#Comb_data_Meth <- Comb_data_Meth[,c(-1,-2)]
#Comb_data_Meth <- as.matrix(Comb_data_Meth)
Comb_data_Meth <- t(Comb_data_Meth)
Comb_data_Meth <- Comb_data_Meth[order(row.names(Comb_data_Meth),decreasing = FALSE),]
FinalDataCond_Meth <- cbind("Condition"=PRADMethids$condition,as.data.frame(Comb_data_Meth))

rm(Comb_data_Meth,Norm_data_Meth,Norm_data_Meth_Cond,normal,cancer)
gc()

# Creates a tab file
setwd("file_path")
write.table(FinalDataCond_Meth, "significantMethPRAD.tab", sep = "\t",col.names = NA, quote = FALSE)

