# Loads the library
library(caret)
library(rpart)
library(rpart.plot)
library(smotefamily)
# Acquires the data
FinalDataCond_miRNA <- read.delim("file_path_significantmiRNA.tab", row.names=1)


# Creates the training and testing set
inTrain <- createDataPartition(FinalDataCond_miRNA$Condition, p = 0.75, list = FALSE, times = 1)
training<-FinalDataCond_miRNA[inTrain,]
testing <-FinalDataCond_miRNA[-inTrain,]
# Prepares the model's parametres
set.seed(50)

train.smote <- SMOTE(training[,-1],training$Condition,K = 5)
train.smote <- train.smote$data # extract only the balanced dataset
train.smote$class <- as.factor(train.smote$class)
prop.table(table(train.smote$class))

control <-trainControl(method='cv', number = 10)
metric <- "Accuracy"
# Trains the three models and gets the confusion matrices with the classification types
set.seed(50)
fit.knn<-train(class~., data=train.smote, method="knn", metric=metric, trControl=control)
predknn<-predict(fit.knn,testing)
confusionMatrix(predknn,as.factor(testing$Condition))
knnRaw<-predict(fit.knn, newdata = testing, type = "raw")
knnRaw<-as.data.frame(knnRaw)
colnames(knnRaw)<- c("Classification Condition")
knnR<-cbind("Original Condition"=testing$Condition, knnRaw)
row.names(knnR)<-rownames(testing)

set.seed(50)
fit.svm<-train(class ~., data=train.smote, method="svmRadial", metric=metric, trControl=control, tuneLength = 10)
predsvm<-predict(fit.svm,testing)
confusionMatrix(predsvm,as.factor(testing$Condition))
svmRaw<-predict(fit.svm, newdata = testing, type = "raw")
svmRaw<-as.data.frame(svmRaw)
colnames(svmRaw)<- c("Classification Condition")
svmR<-cbind("Original Condition"=testing$Condition, svmRaw)
row.names(svmR)<-rownames(testing)

set.seed(50)
fit.rpart<-train(class~., data=train.smote, method="rpart", metric=metric, trControl=control)
predrpart<-predict(fit.rpart,testing)
confusionMatrix(predrpart,as.factor(testing$Condition))
rpartRaw<-predict(fit.rpart, newdata = testing, type = "raw")
rpartRaw<-as.data.frame(rpartRaw)
colnames(rpartRaw)<- c("Classification Condition")
rpartR<-cbind("Original Condition"=testing$Condition, rpartRaw)
row.names(rpartR)<-rownames(testing)
rpart.plot(fit.rpart$finalModel)

# Gets the results of the training models
results<- resamples(list(rpart=fit.rpart, knn=fit.knn, svm=fit.svm))
summary(results)
# Estimates the performance of the 3 models
dotplot(results, main="Model Performance")
diffs<-diff(results)
summary(diffs)
