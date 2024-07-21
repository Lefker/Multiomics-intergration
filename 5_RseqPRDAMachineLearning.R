# Loads the library
library(caret)
library(rpart)
library(rpart.plot)
library(smotefamily)
# Acquires the data
FinalDataCond_Rseq <- read.delim("H:/Ergasia/ssignificantRSeq.tab", row.names=1)

# Prepares the model's parametres
set.seed(50)
Smote_data_Rseq <- SMOTE(FinalDataCond_Rseq[,-1],FinalDataCond_Rseq$Condition,K = 5)
Smote_data_Rseq <- Smote_data_Rseq$data # extract only the balanced dataset
Smote_data_Rseq$class <- as.factor(Smote_data_Rseq$class)
inTrain <- createDataPartition(Smote_data_Rseq$class, p = 0.75, list = FALSE, times = 1)
training_Rseq<-Smote_data_Rseq[inTrain,]
testing_Rseq <-Smote_data_Rseq[-inTrain,]

control <-trainControl(method='cv', number = 25, classProbs = TRUE)
metric <- "Accuracy"
# Trains the three models and gets the confusion matrices with the classification types
set.seed(50)
#modelLookup("knn")
fit.knn_Rseq<-train(class~., data=training_Rseq, method="knn", metric=metric, trControl=control)
predknn_Rseq<-predict(fit.knn_Rseq,testing_Rseq)
confusionMatrix(predknn_Rseq,as.factor(testing_Rseq$class))
knnRaw_Rseq<-predict(fit.knn_Rseq, newdata = testing_Rseq, type = "raw")
knnRaw_Rseq<-as.data.frame(knnRaw_Rseq)
colnames(knnRaw_Rseq)<- c("Classification class")
knnR_Rseq<-cbind("Original class"=testing_Rseq$class, knnRaw_Rseq)
row.names(knnR_Rseq)<-rownames(testing_Rseq)

set.seed(50)
# modelLookup("svmR_Rseqadial")
fit.svm_Rseq<-train(class ~., data=training_Rseq, method="svmRadial", metric=metric, trControl=control, tuneLength = 10)
predsvm_Rseq<-predict(fit.svm_Rseq,testing_Rseq)
confusionMatrix(predsvm_Rseq,as.factor(testing_Rseq$class))
svmRaw_Rseq<-predict(fit.svm_Rseq, newdata = testing_Rseq, type = "raw")
svmRaw_Rseq<-as.data.frame(svmRaw_Rseq)
colnames(svmRaw_Rseq)<- c("Classification Condition")
svmR_Rseq<-cbind("Original Condition"=testing_Rseq$class, svmRaw_Rseq)
row.names(svmR_Rseq)<-rownames(testing_Rseq)

set.seed(50)
# modelLookup("rpart2")
fit.rpart_Rseq<-train(class~., data=training_Rseq, method="rpart2", metric=metric, trControl=control)
predrpart_Rseq<-predict(fit.rpart_Rseq,testing_Rseq)
confusionMatrix(predrpart_Rseq,as.factor(testing_Rseq$class))
rpartRaw_Rseq<-predict(fit.rpart_Rseq, newdata = testing_Rseq, type = "raw")
rpartRaw_Rseq<-as.data.frame(rpartRaw_Rseq)
colnames(rpartRaw_Rseq)<- c("Classification Condition")
rpartR_Rseq<-cbind("Original Condition"=testing_Rseq$class, rpartRaw_Rseq)
row.names(rpartR_Rseq)<-rownames(testing_Rseq)
rpart.plot(fit.rpart_Rseq$finalModel,type = 5, extra = 101, branch.lty = 3, shadow.col = "grey",cex=1)

# Gets the results of the training_Rseq models
results_Rseq<- resamples(list(rpart=fit.rpart_Rseq, knn=fit.knn_Rseq, svm=fit.svm_Rseq))
summary(results_Rseq)
# Estimates the performance of the 3 models
dotplot(results_Rseq, main="Model Performance")
diffs_Rseq<-diff(results_Rseq)
summary(diffs_Rseq)
