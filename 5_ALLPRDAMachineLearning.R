# Loads the library
library(caret)
library(rpart)
library(rpart.plot)
library(smotefamily)
# Acquires the data
FinalDataCond_ALL <- read.delim("file_path//significantPRADALL.tab", row.names=1)

# Prepares the model's parametres
set.seed(50)
Smote_data <- SMOTE(FinalDataCond_ALL[,-1],FinalDataCond_ALL$Condition,K = 5)
Smote_data <- Smote_data$data # extract only the balanced dataset
Smote_data$class <- as.factor(Smote_data$class)
inTrain <- createDataPartition(Smote_data$class, p = 0.75, list = FALSE, times = 1)
training<-Smote_data[inTrain,]
testing <-Smote_data[-inTrain,]

control <-trainControl(method='cv', number = 25, classProbs = TRUE)
metric <- "Accuracy"
# Trains the three models and gets the confusion matrices with the classification types
set.seed(50)
#modelLookup("knn")
fit.knn<-train(class~., data=training, method="knn", metric=metric, trControl=control)
predknn<-predict(fit.knn,testing)
confusionMatrix(predknn,as.factor(testing$class))
knnRaw<-predict(fit.knn, newdata = testing, type = "raw")
knnRaw<-as.data.frame(knnRaw)
colnames(knnRaw)<- c("Classification class")
knnR<-cbind("Original class"=testing$class, knnRaw)
row.names(knnR)<-rownames(testing)

set.seed(50)
# modelLookup("svmRadial")
fit.svm<-train(class ~., data=training, method="svmRadial", metric=metric, trControl=control, tuneLength = 10)
predsvm<-predict(fit.svm,testing)
confusionMatrix(predsvm,as.factor(testing$class))
svmRaw<-predict(fit.svm, newdata = testing, type = "raw")
svmRaw<-as.data.frame(svmRaw)
colnames(svmRaw)<- c("Classification Condition")
svmR<-cbind("Original Condition"=testing$class, svmRaw)
row.names(svmR)<-rownames(testing)

set.seed(50)
# modelLookup("rpart2")
fit.rpart<-train(class~., data=training, method="rpart2", metric=metric, trControl=control)
predrpart<-predict(fit.rpart,testing)
confusionMatrix(predrpart,as.factor(testing$class))
rpartRaw<-predict(fit.rpart, newdata = testing, type = "raw")
rpartRaw<-as.data.frame(rpartRaw)
colnames(rpartRaw)<- c("Classification Condition")
rpartR<-cbind("Original Condition"=testing$class, rpartRaw)
row.names(rpartR)<-rownames(testing)
rpart.plot(fit.rpart$finalModel,type = 5, extra = 101, branch.lty = 3, shadow.col = "grey",cex=1)

# Gets the results of the training models
results<- resamples(list(rpart=fit.rpart, knn=fit.knn, svm=fit.svm))
summary(results)
# Estimates the performance of the 3 models
dotplot(results, main="Model Performance")
diffs<-diff(results)
summary(diffs)
