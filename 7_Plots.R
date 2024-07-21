library(ggplot2)
library(Rtsne)
library(gplots)

# Boxplot of normalized expression matrix

## Needs the Comb_data_**** table from 4_files.R ##
mp<-colorRampPalette(c("yellow" ,"red", "black"))
boxplot(Comb_data_miRna, outline = FALSE, las = 2,col = mp(n=10))

# Boxplot and histogram of normal and cancer datasets

## Needs the normal and cancer tables from 4_ files.R )
boxplot(normal,las=2,outline=FALSE,col=mp(n=10))
boxplot(cancer,las=2,outline=FALSE,col=mp(n=10))
hist(normal)
hist(cancer)



# Dimensionality reduction with t-SNE

## Needs the Norm_data_****  from 4_files.R ##
Norm_data_miRna <- t(Norm_data_miRna)
PRADmiRNA_Seq <- t(PRADmiRNA_Seq)
tSNEall<-Rtsne(PRADmiRNA_Seq,dims = 2, normalize = FALSE, pca_center = FALSE, max_iter = 5000, perplexity = 50,theta = 0.0)

# Dimensionality reduction with t-SNE on the selected data

#Comb_data_Rseq <- t(Comb_data_Rseq)
## Needs the Comb_data_***** from 4_files.R ##
tSNEfinal <- Rtsne(Comb_data_miRna, dims = 2, normalize = FALSE, pca_center = FALSE, max_iter = 5000, perplexity=50, theta = 0.0)

# Visualization of t-SNE results

## Needs the PRADALL_cond or Norm_data_****_Cond from 4_files.R
tsne_plot1 <- data.frame(x = tSNEall$Y[,1], y = tSNEall$Y[,2], col = Norm_data_miRna_Cond$Condition)
ggplot(tsne_plot1,aes(x=x, y=y, color=col)) + geom_point() + labs(title = "Initial Tsne Plot") + theme(axis.title = element_blank(), axis.title.y = element_blank())

# Visualization of t-SNE results to the selected data
tsne_plot2 <- data.frame(x = tSNEfinal$Y[,1], y = tSNEfinal$Y[,2], col = FinalDataCond_miRNA$Condition)
ggplot(tsne_plot2,aes(x=x, y=y, color=col)) + geom_point() + labs(title = "Final Tsne Plot") + theme(axis.title = element_blank(), axis.title.y = element_blank())



# Density and Histograms of pvalues and adjusted pvalues

## Needs the pvalue and padh from 4_files.R
ggplot(as.data.frame(pvalue), aes(x=pvalue)) + 
  geom_histogram(aes(y=after_stat(density)), colour="black", fill="red",bins = 15)+
  geom_density(alpha=.5, fill="yellow", linewidth= .2, linetype = 2) +
  labs(
    title = "Density plot and Histogram of pvalues",
    x = "pvalues",
    y = "Counts"
  )

ggplot(as.data.frame(padj), aes(x=padj)) + 
  geom_histogram(aes(y=after_stat(density)), colour="black", fill="white",bins = 15)+
  geom_density(alpha=.3, fill="blue") +
  labs(
    title = "Density plot and Histogram of Adjusted pvalues",
    x = "Adjusted pvalues",
    y = "Counts"
  )



# Histogram of Fold change

## Needs the foldchange from 4_files.R
hist(foldchange)
hist(abs(foldchange))
ggplot(as.data.frame(foldchange), aes(x=abs(FoldChange))) + 
  geom_histogram(aes(y=after_stat(density)),color="black", fill="white", binwidth = 0.2,position = "dodge") + 
  geom_density(alpha=.2, fill="yellow", linewidth= .2, linetype = 2) +
  labs(
    title = "Density plot and Histogram of Multiomics Absolute of FoldChange",
    x = "FoldChange",
    y = "Density"
  ) + theme_minimal()

# Gets the heatmap data
#jpeg(file="saving_plot4.jpg",width= 3000, height = 1200)
mp1 <-colorRampPalette(c("red","white", "green"))
heat_data <- t(as.data.frame(lapply(FinalDataCond_miRNA[,-1], as.numeric)))
heatmap.2(as.matrix(heat_data), col=mp1(n=10), key.xlab = paste("log2data"), 
          xlab = "Samples", dendrogram = "none", Colv = "NA" , ylab = "Features", 
          trace=c("none"), main = "Multiomics Expression Data", 
          labCol = FALSE, labRow = FALSE,
          keysize = 0.85, ColSideColors = rep(c("purple","black"),c(493,35))) 
par(lend = 2)           # square line ends for the color legend
legend(y=1.073, x=.888, xpd=TRUE,      # location of the legend on the heatmap plot
       legend = c("Cancer", "Normal"), # category labels
       col = c("purple", "black"),  # color key
       lty= 1,             # line style
       lwd = 15            # line width
)


