# Loads the tabs of the expression values and the patient's ids
PRADmiRNA_Seq <- read.delim("H:/Ergasia/PRADmiRNA_Seq.tab", row.names=1)
names(PRADmiRNA_Seq) <- gsub(x = names(PRADmiRNA_Seq), pattern = "\\.", replacement = "-")
PRADmiRNA_Seq <- PRADmiRNA_Seq[!apply(PRADmiRNA_Seq, 1, function(x) all(x == 0)), ]
PRADrow <- rownames(PRADmiRNA_Seq)
PRADcol <- colnames(PRADmiRNA_Seq)
PRADmiRNAids <- read.delim("H:/Ergasia/PRADmiRNAids.tab", row.names=1)
PRADmiRNA_Seq <- matrix(unlist(PRADmiRNA_Seq), ncol = nrow(PRADmiRNAids))
rownames(PRADmiRNA_Seq) <- PRADrow
colnames(PRADmiRNA_Seq) <- PRADcol


# Standardization by row and by column
# Normalize the expression matrix

Norm_data_miRna <- norm_fun(PRADmiRNA_Seq)

# splits the ids to healthy and disease
healthyids <- PRADmiRNAids$patient[which(grepl("normal",PRADmiRNAids$condition))]
cancerids <- PRADmiRNAids$patient[which(grepl("cancer",PRADmiRNAids$condition))]

# Splits the dataset to normal and disease
normal <- Norm_data_miRna[,which(colnames(Norm_data_miRna) %in% healthyids)]
cancer <- Norm_data_miRna[,which(colnames(Norm_data_miRna) %in% cancerids)]

# Combines with the condition
Norm_data_miRna <- t(Norm_data_miRna)
Norm_data_miRna <- Norm_data_miRna[order(row.names(Norm_data_miRna)), ]
PRADmiRNAids <- PRADmiRNAids[order(PRADmiRNAids$patient),]
Norm_data_miRna_Cond <- cbind("Condition"=PRADmiRNAids$condition, as.data.frame(Norm_data_miRna))
# Norm_data_miRna<-as.data.frame(Norm_data_miRna)
Norm_data_miRna_Cond <- as.data.frame(Norm_data_miRna_Cond)

# Gets the p-value from student's t-test
pvalue <- vector()
for( i in 1:nrow(normal)){
  stat.test <- t.test(normal[i,], cancer[i,],var.equal = TRUE)
  pvalue[i] <- stat.test$p.value
}
Norm_data_miRna <- t(Norm_data_miRna)
pvalue <- as.matrix(pvalue)
rownames(pvalue)<-rownames(Norm_data_miRna)
colnames(pvalue)<-c("pvalue")

# Adjusts the p-value with fdr method
set.seed(42)
padj <- p.adjust(pvalue, method = "fdr")
padj <- as.matrix(padj)
rownames(padj) <- rownames(Norm_data_miRna)
colnames(padj) <- c("p_adjusted")

# Function of the fold change.
# In brackets is needed in the order: normal table, cancer table, combined table
foldchange <- FoldChange_fun(normal,cancer,Norm_data_miRna)

# Combines the p-adjusted value and the normalized expression data
Comb_data_miRna <- cbind(padj,foldchange,as.data.frame(Norm_data_miRna))
Comb_data_miRna <- Comb_data_miRna[order(Comb_data_miRna$FoldChange,decreasing = TRUE),]

# Orders of the table based on the adjusted p-value 
# and keeps the features with adjusted p_value < 0.05
Comb_data_miRna <- subset(Comb_data_miRna, p_adjusted < 0.05 & abs(FoldChange) >=1.5, select = c(-1,-2))
Comb_data_miRna <- as.data.frame(Comb_data_miRna)

Comb_data_miRna <- t(Comb_data_miRna)
Comb_data_miRna <- Comb_data_miRna[order(row.names(Comb_data_miRna),decreasing = FALSE),]
Comb_data_miRna_Cond <- cbind("Condition"=PRADmiRNAids$condition,as.data.frame(Comb_data_miRna))



# Creates a tab file
setwd("H:/Ergasia")
write.table(Comb_data_miRna_Cond, "significantmiRNA.tab", sep = "\t", quote = FALSE)
