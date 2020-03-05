
# Using the same structure than in gmb script approach 2, I tried different algorithms:

# load libraries
library(caret)
library(pROC)
library(plyr)
library(ROCR)

# load data
# If necessary, convert the output variable into factor.
#pgpAgeGenderSNPs$Status <- as.factor(pgpAgeGenderSNPs$Status)

pgpDF <- pgpAgeGenderSNPs
#pgpDF <- PGPFeatA2SNPs_v2


# Here we select which variables we want to consider as features (feature selection) for
# the classification.

# pgpDF <- pgpDF[,c('Age','Sex','rs...',..., 'Status')]# Status will be our output
# Set 1: SNPs proposed in paper 2
#pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134","rs10195252")]

# Set 2:SNPs paper2 + 2 new: "rs10090154" and "rs1447295" (In svm, AUC=0.8683)
# pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134" , "rs10195252", "rs10090154", "rs1447295")]
# Set 3:(In svm, AUC=0.8645)
# pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134" , "rs10195252", "rs10090154", "rs1447295","rs10045413","rs2007084")]
# Set 4: Set proposed in paper2 + "rs17104630". This SNP was suggested by Boruta after correcting values from PGPFeaturesV2.csv(new PGPFeaturesV2_v2.csv)(In svm, AUC=0.8875)
# pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134","rs17104630" )]# "rs10195252" removed for testing purposes, it seems to be not relevant.
# Set 5: same as set 2 + "rs17104630" (In svm, AUC=0.8862) 
# pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134" , "rs10195252", "rs10090154", "rs1447295","rs17104630")]

# Set 6: Same as set4 + two new SNPs suggested by Boruta: "rs1447295" and "rs4242382".(In svm, AUC=0.8887)
# pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134","rs17104630","rs1447295", "rs4242382" )]

# Set 7: Set 6 + "rs10090154"(In svm, AUC=0.8913)
# pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134","rs10195252","rs17104630","rs1447295", "rs4242382","rs10090154" )]

# Set 8 (gets the best AUC results so far): Set 7 + "Age" (AUC=0.9028)
#  pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134","rs10195252","rs17104630","rs1447295", "rs4242382","rs10090154","Age")]

# Set 9: Set 7 + "rs12978500" (In svm, AUC=0.9054)
 pgpDF <- pgpDF[,c("Status","rs12567355", "rs3001167", "rs586688", "rs11241130", "rs7525133", "rs2076529", "rs12970134","rs10195252","rs17104630","rs1447295", "rs4242382","rs10090154","rs12978500")]


# test results from univariate filter: pgpDF <- pgpDF[,c("Status","rs10090154", "rs10506525" , "rs11074889", "rs12101726", "rs12567355")]

# generalize outcome and predictor variables for easy recycling of the code for subsequents needs.
outcomeName <- 'Status'
predictorsNames <- names(pgpDF)[names(pgpDF) != outcomeName]

#################################################
# model it
#################################################
# get names of all caret supported models. Useful if we need to know other available models. 
# names(getModelInfo())

pgpDF$Status <- ifelse(pgpDF$Status==1,'Normal','Risk') # Assign 'Normal' to Status=1 and 'Risk' to Status=2
pgpDF$Status <- as.factor(pgpDF$Status)

# we need to split our data into two portions: a training and a testing portion.
# The function createDataPartition can be used to create a stratified random sample of the data into training and test sets
# Setting the seed is paramount for reproducibility as createDataPartition will shuffle the data randomly before splitting it.

set.seed(1234)
splitIndex <- createDataPartition(pgpDF[,outcomeName], p = .75, list = FALSE, times = 1)
trainDF <- pgpDF[ splitIndex,]
testDF  <- pgpDF[-splitIndex,]

# MODEL TRAINING 

# Here we are declaring the number of folds for the training function
FoldNumber <- 10

# NOTE: I removed returnResamp='none' from all trainControl() in order to use resamples() with all the results at the end of this script.

# 1- gbm: Gradient Boost Model (Gradient Boosting Classifier Algorithm)

gbmctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)

# run model
set.seed(1234)
gbm.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]), 
                  method='gbm', 
                  trControl=gbmctrl,  
                  metric = "ROC",
                  preProc = c("center", "scale"))

# Model details
gbm.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(gbm.mod, results[order(results$ROC, decreasing=T),]))

# Variable Importance
gbmImp <- varImp(gbm.mod, scale = FALSE)
gbmImp

# Variable importance plot
plot(varImp(gbm.mod),main="gbm Variable Importance")

### gbm Model Predictions and Performance
# Make predictions using the test data set
gbm.pred <- predict(gbm.mod,testDF[,predictorsNames])

#Look at the confusion matrix  
confusionMatrix(gbm.pred,testDF$Status, positive = "Risk")   

#Draw the ROC curve 
gbm.probs <- predict(gbm.mod,testDF[,predictorsNames],type="prob")
head(gbm.probs)

gbm.ROC <- roc(predictor=gbm.probs$Normal,
                 response=testDF$Status,
                 levels=levels(testDF$Status))
gbm.ROC$auc

plot(gbm.ROC,main="gbm ROC")

# It is also possible to add a threshold value in the ROC plot. In this case, 50%.
#plot(gbm.ROC, type="S", print.thres = .5) 

# 2- glmnet: Generalized Linear Model Lasso and Elastic-Net Regularized
glmnetctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(1234)
glmnet.mod <- train(trainDF[,predictorsNames], trainDF[,outcomeName], method='glmnet',  metric = "ROC", trControl=glmnetctrl)

# Model details
glmnet.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(glmnet.mod, results[order(results$ROC, decreasing=T),]))

# Variable Importance
glmnetImp <- varImp(glmnet.mod, scale = FALSE)
glmnetImp
# Variable importance plot
plot(varImp(glmnet.mod),main="glmnet Variable Importance")

### glmnet Model Predictions and Performance
# Make predictions using the test data set
glmnet.pred <- predict(glmnet.mod,testDF[,predictorsNames])

#Look at the confusion matrix  
confusionMatrix(glmnet.pred,testDF$Status,positive = "Risk")   

#Draw the ROC curve 
glmnet.probs <- predict(glmnet.mod,testDF[,predictorsNames],type="prob")
head(glmnet.probs)

glmnet.ROC <- roc(predictor=glmnet.probs$Normal,
               response=testDF$Status,
               levels=levels(testDF$Status))
glmnet.ROC$auc

plot(glmnet.ROC,main="glmnet ROC")

# 3- CLASSIFICATION AND REGRESION TREES (CART)
# create caret trainControl object to control the number of cross-validations performed
rpartctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)

# With this tunning parameters I obtained a very good AUC in SVM
#objControl <- trainControl(method = "repeatedcv", repeats = 5, summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(1234)
rpart.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]),
                  method="rpart", metric="ROC", trControl=rpartctrl,
                  preProc = c("center", "scale"))

# Model details
rpart.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(rpart.mod, results[order(results$ROC, decreasing=T),]))

# Variable Importance
rpartImp <- varImp(rpart.mod, scale = FALSE)
rpartImp
# Variable importance plot
plot(varImp(rpart.mod),main="CART Variable Importance")

### rpart Model Predictions and Performance
# Make predictions using the test data set
rpart.pred <- predict(rpart.mod,testDF[,predictorsNames])

#Look at the confusion matrix  
confusionMatrix(rpart.pred,testDF$Status,positive = "Risk")   

#Draw the ROC curve 
rpart.probs <- predict(rpart.mod,testDF[,predictorsNames],type="prob")
head(rpart.probs)

rpart.ROC <- roc(predictor=rpart.probs$Normal,
               response=testDF$Status,
               levels=levels(testDF$Status))
rpart.ROC$auc

plot(rpart.ROC,main="CART ROC")

# 4- k-Nearest Neighbors (kNN).
knnctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(1234)
knn.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]),
                  method="knn", metric="ROC", trControl=knnctrl,
                  preProc = c("center", "scale"))

# Model details
knn.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(knn.mod, results[order(results$ROC, decreasing=T),]))

# Variable Importance
knnImp <- varImp(knn.mod, scale = FALSE)
knnImp
# Variable importance plot
plot(varImp(knn.mod),main="knn Variable Importance")

# Make predictions using the test data set
knn.pred <- predict(knn.mod,testDF[,predictorsNames])
confusionMatrix(knn.pred,testDF$Status,positive = "Risk")   

#Draw the ROC curve 
knn.probs <- predict(knn.mod,testDF[,predictorsNames],type="prob")
head(knn.probs)

knn.ROC <- roc(predictor=knn.probs$Normal,
                 response=testDF$Status,
                 levels=levels(testDF$Status))
knn.ROC$auc

plot(knn.ROC,main="knn ROC")

# 5- Support Vector Machines (SVM) with a linear kernel. [svmRAdial --> Support Vector Machines with Radial Basis Function Kernel]
SVMctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(1234)
svm.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]),
                  method="svmRadial", metric="ROC", trControl=SVMctrl,
                  preProc = c("center", "scale"))

# Model details
svm.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(svm.mod, results[order(results$ROC, decreasing=T),]))

# Variable Importance
svmImp <- varImp(svm.mod, scale = FALSE)
svmImp
# Variable importance plot
plot(varImp(svm.mod),main="svm Variable Importance")

# Make predictions using the test data set
svm.pred <- predict(svm.mod,testDF[,predictorsNames])
confusionMatrix(svm.pred,testDF$Status,positive = "Risk")   

#Draw the ROC curve 
svm.probs <- predict(svm.mod,testDF[,predictorsNames],type="prob")
head(svm.probs)

svm.ROC <- roc(predictor=svm.probs$Normal,
               response=testDF$Status,
               levels=levels(testDF$Status))
svm.ROC$auc

plot(svm.ROC,main="svm ROC")

# 6- Random Forest (RF)
RFctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)

set.seed(1234)
rf.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]),
                  method="rf", metric="ROC", trControl=RFctrl,
                  preProc = c("center", "scale"))

# Model details
rf.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(rf.mod, results[order(results$ROC, decreasing=T),]))

# Variable Importance
rfImp <- varImp(rf.mod, scale = FALSE)
rfImp
# Variable importance plot
plot(varImp(rf.mod),main="rf Variable Importance")

# Make predictions using the test data set
rf.pred <- predict(rf.mod,testDF[,predictorsNames])
confusionMatrix(rf.pred,testDF$Status,positive = "Risk")   

#Draw the ROC curve 
rf.probs <- predict(rf.mod,testDF[,predictorsNames],type="prob")
head(rf.probs)

rf.ROC <- roc(predictor=rf.probs$Normal,
               response=testDF$Status,
               levels=levels(testDF$Status))
rf.ROC$auc

plot(rf.ROC,main="rf ROC")

# 7- Neural Network

# If we want to include nnet in the resamples list created when comparing multiple models, we need to use the same trainControl parameters than in the other models (metric must be ROC).
nnetctrl <- trainControl(method='cv', number=FoldNumber, summaryFunction = twoClassSummary, classProbs = TRUE)
set.seed(1234)
nnet.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]), method = "nnet", metric = "ROC", trControl = nnetctrl, algorithm = "backprop", linear.output = T) # Algorithm --> 'backprop'  Backpropagation neural network.

# nnetctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 30)
# set.seed(1234)
# nnet.mod <- train(trainDF[,predictorsNames], as.factor(trainDF[,outcomeName]), method = "nnet", metric = "Accuracy", trControl = nnetctrl, algorithm = "rprop+", linear.output = T)

# Model details
nnet.mod

# Results orederd from top to bottom (of the ROC results). The top 6 results are printed using head().
head(with(nnet.mod, results[order(results$ROC, decreasing=T),]))
# For Accuracy will be this instead...
#head(with(nnet.mod, results[order(results$Accuracy, decreasing=T),]))

# Variable Importance
nnetImp <- varImp(nnet.mod, scale = FALSE)
nnetImp
# Variable importance plot
plot(varImp(nnet.mod),main="nnet Variable Importance")

# Make predictions using the test data set
nnet.pred <- predict(nnet.mod,testDF[,predictorsNames])
confusionMatrix(nnet.pred,testDF$Status,positive = "Risk")   

#Draw the ROC curve 
nnet.probs <- predict(nnet.mod,testDF[,predictorsNames],type="prob")
head(nnet.probs)

nnet.ROC <- roc(predictor=nnet.probs$Normal,
               response=testDF$Status,
               levels=levels(testDF$Status))
nnet.ROC$auc

plot(nnet.ROC,main="nnet ROC")

# Comparing Multiple Models
# Having set the same seed before running the models
# we have generated paired samples and are in a position to compare models 
# using a resampling technique.
# (See Hothorn et al, "The design and analysis of benchmark experiments
# -Journal of Computational and Graphical Statistics (2005) vol 14 (3) 
# pp 675-699) 

rValues <- resamples(list(gbm=gbm.mod,glmnet=glmnet.mod ,rpart=rpart.mod,knn=knn.mod, svm=svm.mod, rf=rf.mod, nnet=nnet.mod))
rValues$values
summary(rValues)

# ROC summary (predictions):

ROCsummary <- data.frame(Model=c("gbm", "glmnet", "rpart", "knn", "svm", "rf", "nnet"),ROC=c(gbm.ROC$auc, glmnet.ROC$auc, rpart.ROC$auc,knn.ROC$auc, svm.ROC$auc, rf.ROC$auc, nnet.ROC$auc))

ROCsummary

# Plots
bwplot(rValues,metric="ROC",main="model comparison")	# boxplot
dotplot(rValues,metric="ROC",main="model comparison")	# dotplot



#Better plots for comparison
scales <- list(x=list(relation="free"), y=list(relation="free"))

bwplot(rValues, scales=scales,main="Boxplot: model comparison using ROC as metric")

dotplot(rValues,scales=scales,main="Dotplot: model comparison using ROC as metric")

#This create a scatterplot matrix of all fold-trial results for an algorithm compared to the same fold-trial results for all other algorithms. All pairs are compared.
splom(rValues,scales=scales)

#1# All ROC curves in one plot:

plot(gbm.ROC, col= 'red',lty=1)
par(new=T)
plot(glmnet.ROC, col='green',lty=2)
par(new=T)
plot(rpart.ROC, col='blue',lty=3)
par(new=T)
plot(knn.ROC,col='darkorange1',lty=4)
par(new=T)
plot(svm.ROC,col= 'navyblue',lty=5)
par(new=T)
plot(rf.ROC,col='gray46',lty=6)
par(new=T)
plot(nnet.ROC,main="ROC Comparison",col='purple',lty=7)
# Custom legend location in the plot
legend(0.2,0.6,c("Gbm","Glmnet","CART","Knn","SVM","RF","Nnet"), lty=1:7,lwd=2.5,col=c("red","green","blue","darkorange1","navyblue","gray46","purple"), cex=.75, bty='n')

#legend('right',c("Gbm","Glmnet","CART","Knn","SVM","RF","Nnet"), lty=1,lwd=2.5,col=c("red","green","blue","orange","yellow","grey","purple"), bty='n')
#legend('bottomright', c("Gbm","Glmnet","Rpart","Knn","SVM","RF","Nnet"), lty=1, col=c('red', 'blue', 'green',' brown'), bty='n', cex=.75)

#2# All ROC curves in the same plot but individually represented:

# We prepare the plot for the 7 ROC graphs. We create 3 rows and 3 columns
par(mfrow=c(3,3))

plot(gbm.ROC,main="gbm ROC", col= 'red')
plot(glmnet.ROC,main="glmnet ROC", col='green')
plot(rpart.ROC,main="rpart ROC", col='blue')
plot(knn.ROC,main="knn ROC",col='orange')
plot(svm.ROC,main="svm ROC",col= 'yellow')
plot(rf.ROC,main="rf ROC",col='grey')
plot(nnet.ROC,main="nnet ROC",col='purple')

# Variable Importance plots:

par(mfrow=c(3,3))

 plot(gbmImp,main="gbm Variable Importance", col= 'red')
 plot(glmnetImp,main="glmnet Variable Importance", col='green')
 plot(rpartImp,main="CART Variable Importance", col='blue')
 plot(knnImp,main="knn Variable Importance",col='orange')
 plot(svmImp,main="svm Variable Importance",col= 'yellow')
 plot(rfImp,main="rf Variable Importance",col='grey')
 plot(nnetImp,main="nnet Variable Importance",col='purple')



# MODEL PREDICTION PLOTS (confusion matrix representation) - It looks nice to see participants predictions based on the two classes. 
# web link example: http://stats.stackexchange.com/questions/109079/confusion-between-caret-randomforest-predict-results-and-reported-model-perfor
# example with gbm:

results <- data.frame(pred = predict(gbm.mod,testDF[,predictorsNames]), obs = testDF[,outcomeName])

p <- ggplot(results, aes(x = pred, y = obs))
p <- p + geom_jitter(position = position_jitter(width = 0.25, height = 0.25))
p
