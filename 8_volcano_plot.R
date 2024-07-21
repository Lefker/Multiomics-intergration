library(ggrepel)
#123123123
dataVolc <- cbind(padj,foldchange,as.data.frame(Norm_data_Rseq))
dataVolc <- cbind(padj,foldchange,as.data.frame(PRADall))
dataVolc <- cbind(padj,foldchange,as.data.frame(Norm_data_miRna))
dataVolc <- cbind(padj,foldchange,as.data.frame(Norm_data_Meth))

ggplot(data=as.data.frame(dataVolc), aes(x=FoldChange, y=-log10(p_adjusted))) + 
  geom_point() + theme_minimal() + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# add a column of NAs
dataVolc$diffexpressed <- "NO"
# if log2Foldchange >= 1.5 and pvalue < 0.05, set as "UP" 
dataVolc$diffexpressed[dataVolc$FoldChange >= 1.5 & dataVolc$p_adjusted < 0.05] <- "UP"
# if log2Foldchange <= -0.6 and pvalue < 0.05, set as "DOWN"
dataVolc$diffexpressed[dataVolc$FoldChange <= -1.5 & dataVolc$p_adjusted < 0.05] <- "DOWN"

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("red", "darkgreen", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")

# Create a new column "delabel" to dataVolc, that will contain the name of genes differentially expressed (NA in case they are not)
dataVolc$delabel <- NA
dataVolc$delabel[dataVolc$diffexpressed != "NO"] <- rownames(dataVolc[dataVolc$diffexpressed != "NO",])

# Volcano plot with dots text
ggplot(data=as.data.frame(dataVolc), aes(x=FoldChange, y=-log10(p_adjusted), col=diffexpressed, label = delabel)) + 
  geom_point() + theme_minimal() + scale_color_manual(values = mycolors) + geom_text_repel(max.overlaps = Inf,hjust=1,direction = "y") +
  geom_vline(xintercept=c(-1.5, 1.5), col="blue") +
  geom_hline(yintercept=-log10(0.05), col="blue")

# Volcano plot without dots text
ggplot(data=as.data.frame(dataVolc), aes(x=FoldChange, y=-log10(p_adjusted), col=diffexpressed, label = delabel)) + 
  geom_point(size = 1.5) + theme_minimal() + scale_color_manual(values = mycolors) + 
  geom_vline(xintercept=c(-1.5, 1.5), col="blue") +
  geom_hline(yintercept=-log10(0.05), col="blue") + 
  labs(title = "Volcano Plot for Multiomics", color= "Diffential Expression", x = expression("log"[2]*"FC"), y = expression("-log"[10]*"FDR")) +
  theme(plot.title = element_text(size = 25, face = "bold",hjust = 0.5)
  ) + xlim(-3,3)
