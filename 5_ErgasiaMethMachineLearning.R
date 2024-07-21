# Loads the library
library(caret)
library(rpart)
library(rpart.plot)
library(smotefamily)
# Acquires the data
FinalDataCond_Meth <- read.delim("file_path_significantMethPRAD.tab", row.names=1)

# Prepares the model's parametres
set.seed(50)
Smote_data_Meth <- SMOTE(FinalDataCond_Meth[,-1],FinalDataCond_Meth$Condition,K = 5)
Smote_data_Meth <- Smote_data_Meth$data # extract only the balanced dataset
Smote_data_Meth$class <- as.factor(Smote_data_Meth$class)
inTrain <- createDataPartition(Smote_data_Meth$class, p = 0.75, list = FALSE, times = 1)
training_Meth<-Smote_data_Meth[inTrain,]
testing_Meth <-Smote_data_Meth[-inTrain,]

control <-trainControl(method='cv', number = 25, classProbs = TRUE)
metric <- "Accuracy"
# Trains the three models and gets the confusion matrices with the classification types
set.seed(50)
#modelLookup("knn")
fit.knn_Meth<-train(class~., data=training_Meth, method="knn", metric=metric, trControl=control)
predknn_Meth<-predict(fit.knn_Meth,testing_Meth)
confusionMatrix(predknn_Meth,as.factor(testing_Meth$class))
knnRaw_Meth<-predict(fit.knn_Meth, newdata = testing_Meth, type = "raw")
knnRaw_Meth<-as.data.frame(knnRaw_Meth)
colnames(knnRaw_Meth)<- c("Classification class")
knnR_Meth<-cbind("Original class"=testing_Meth$class, knnRaw_Meth)
row.names(knnR_Meth)<-rownames(testing_Meth)

set.seed(50)
# modelLookup("svmR_Methadial")
fit.svm_Meth<-train(class ~., data=training_Meth, method="svmRadial", metric=metric, trControl=control, tuneLength = 10)
predsvm_Meth<-predict(fit.svm_Meth,testing_Meth)
confusionMatrix(predsvm_Meth,as.factor(testing_Meth$class))
svmRaw_Meth<-predict(fit.svm_Meth, newdata = testing_Meth, type = "raw")
svmRaw_Meth<-as.data.frame(svmRaw_Meth)
colnames(svmRaw_Meth)<- c("Classification Condition")
svmR_Meth<-cbind("Original Condition"=testing_Meth$class, svmRaw_Meth)
row.names(svmR_Meth)<-rownames(testing_Meth)

set.seed(50)
# modelLookup("rpart2")
fit.rpart_Meth<-train(class~., data=training_Meth, method="rpart2", metric=metric, trControl=control)
predrpart_Meth<-predict(fit.rpart_Meth,testing_Meth)
confusionMatrix(predrpart_Meth,as.factor(testing_Meth$class))
rpartRaw_Meth<-predict(fit.rpart_Meth, newdata = testing_Meth, type = "raw")
rpartRaw_Meth<-as.data.frame(rpartRaw_Meth)
colnames(rpartRaw_Meth)<- c("Classification Condition")
rpartR_Meth<-cbind("Original Condition"=testing_Meth$class, rpartRaw_Meth)
row.names(rpartR_Meth)<-rownames(testing_Meth)
rpart.plot(fit.rpart_Meth$finalModel,type = 5, extra = 101, branch.lty = 3, shadow.col = "grey",cex=1)

# Gets the results of the training_Meth models
results_Meth<- resamples(list(rpart=fit.rpart_Meth, knn=fit.knn_Meth, svm=fit.svm_Meth))
summary(results_Meth)
# Estimates the performance of the 3 models
dotplot(results_Meth, main="Model Performance")
diffs_Meth<-diff(results_Meth)
summary(diffs_Meth)
