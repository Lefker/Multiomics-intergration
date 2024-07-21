# Loads the library
library(caret)
library(rpart)
library(rpart.plot)
library(smotefamily)
# Acquires the data
FinalDataCond_miRNA <- read.delim("K:\\METAPTYXIAKO\\Diplomatikh\\Ergasia\\significantmiRNA.tab", row.names=1)


# Prepares the model's parametres
set.seed(50)
Smote_data_miRNA <- SMOTE(FinalDataCond_miRNA[,-1],FinalDataCond_miRNA$Condition,K = 5)
Smote_data_miRNA <- Smote_data_miRNA$data # extract only the balanced dataset
Smote_data_miRNA$class <- as.factor(Smote_data_miRNA$class)
inTrain <- createDataPartition(Smote_data_miRNA$class, p = 0.75, list = FALSE, times = 1)
training_miRNA<-Smote_data_miRNA[inTrain,]
testing_miRNA <-Smote_data_miRNA[-inTrain,]

control <-trainControl(method='cv', number = 10)
metric <- "Accuracy"
# Trains the three models and gets the confusion matrices with the classification types
set.seed(50)
fit.knn<-train(class~., data=training_miRNA, method="knn", metric=metric, trControl=control)
predknn<-predict(fit.knn,testing_miRNA)
confusionMatrix(predknn,as.factor(testing_miRNA$class))
knnRaw<-predict(fit.knn, newdata = testing_miRNA, type = "raw")
knnRaw<-as.data.frame(knnRaw)
colnames(knnRaw)<- c("Classification Condition")
knnR<-cbind("Original Condition"=testing_miRNA$class, knnRaw)
row.names(knnR)<-rownames(testing_miRNA)

set.seed(50)
fit.svm<-train(class ~., data=training_miRNA, method="svmRadial", metric=metric, trControl=control, tuneLength = 10)
predsvm<-predict(fit.svm,testing_miRNA)
confusionMatrix(predsvm,as.factor(testing_miRNA$class))
svmRaw<-predict(fit.svm, newdata = testing_miRNA, type = "raw")
svmRaw<-as.data.frame(svmRaw)
colnames(svmRaw)<- c("Classification Condition")
svmR<-cbind("Original Condition"=testing_miRNA$class, svmRaw)
row.names(svmR)<-rownames(testing_miRNA)

set.seed(50)
fit.rpart<-train(class~., data=training_miRNA, method="rpart", metric=metric, trControl=control)
predrpart<-predict(fit.rpart,testing_miRNA)
confusionMatrix(predrpart,as.factor(testing_miRNA$class))
rpartRaw<-predict(fit.rpart, newdata = testing_miRNA, type = "raw")
rpartRaw<-as.data.frame(rpartRaw)
colnames(rpartRaw)<- c("Classification Condition")
rpartR<-cbind("Original Condition"=testing_miRNA$class, rpartRaw)
row.names(rpartR)<-rownames(testing_miRNA)
rpart.plot(fit.rpart$finalModel)

# Gets the results of the training models
results<- resamples(list(rpart=fit.rpart, knn=fit.knn, svm=fit.svm))
summary(results)
# Estimates the performance of the 3 models
dotplot(results, main="Model Performance")
diffs<-diff(results)
summary(diffs)
