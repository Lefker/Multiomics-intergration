# Gets the results of the training_Meth models
results_knn <- resamples(list(Rseq_knn=fit.knn_Rseq, miRNA_knn=fit.knn_miRNA, Meth_knn=fit.knn_Meth, Multi_knn=fit.knn))
results_svm <- resamples(list(Rseq_svm=fit.svm_Rseq, miRNA_svm=fit.svm_miRNA, Meth_svm=fit.svm_Meth, Multi_svm=fit.svm))
results_rpart <- resamples(list(Rseq_rpart=fit.rpart_Rseq, miRNA_rpart=fit.rpart_miRNA, Meth_rpart=fit.rpart_Meth, Multi_rpart=fit.rpart))
results_final <- resamples(list(Rseq_knn=fit.knn_Rseq, miRNA_knn=fit.knn_miRNA, Meth_knn=fit.knn_Meth, Multi_knn=fit.knn,
                          Rseq_svm=fit.svm_Rseq, miRNA_svm=fit.svm_miRNA, Meth_svm=fit.svm_Meth, Multi_svm=fit.svm,
                          Rseq_rpart=fit.rpart_Rseq, miRNA_rpart=fit.rpart_miRNA, Meth_rpart=fit.rpart_Meth, Multi_rpart=fit.rpart))

# Comparison of models performances
summary(results_knn)
summary(results_svm)
summary(results_rpart)
summary(results_final)

# Dotplot of the models performances comparison
dotplot(results_knn, main="Knn Performance")
dotplot(results_svm, main="SVM Performance")
dotplot(results_rpart, main="Rpart Performance")
dotplot(results_final, main="Performance",scales=list(cex=1.4))

# Parallel plot of accuracy
parallelplot(results_final)
parallelplot(results_final, metric = "Kappa")

# Shows the differences of the models performances
diffs_knn<-diff(results_knn)
diffs_svm<-diff(results_svm)
diffs_rpart<-diff(results_rpart)

# Shows the differences of the models
summary(diffs_knn)
summary(diffs_svm)
summary(diffs_rpart)

library(xlsx)
write.xlsx(svm_resample, "svm_resample.xlsx", sheetName = "svm_resample", col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(svm_resample_Rseq, "svm_resample_Rseq.xlsx", sheetName = "svm_resample_Rseq", col.names = TRUE, row.names = TRUE, append = FALSE)
