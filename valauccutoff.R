valauccutoff <- function(data, seed = 100, output = "AUC", alg = "RF", hid = 200, mtr = round(sqrt(ncol(data)-1)), ntree = 1000, nodesize = 2){
  
  #library(mxnet)
  
  library(caret)
  
  library(pROC)
  
  library(ROCR)
  
  library(RODBC)
  
  library(e1071)
  
  library(randomForest)
  
  k = 5
  
  r = 1
  
  data[,1] <- as.factor(data[,1])
  
  data[,1] <- as.character(data[,1])
  
  data[,1] <- as.factor(data[,1])
  
  ultauc <- matrix(0, ncol = k, nrow = r)
  
  ultfm <- matrix(0, ncol = k, nrow = r)
  
  ultbr <- matrix(0, ncol = k, nrow = r)
  
  datalist <- list()
  
  crosspred <- as.data.frame(matrix(0, nrow = 0, ncol = length(levels(data$Type))+1))
  
  for (j in 1:r){
    
    if (featsel == FALSE){
      
      print(paste("Commencing repeat...       [", j, "/", r, "]", sep = ""))
      
    }
    
    set.seed(seed[j])
    
    obs <- length(data[,1])
    
    for (z in 1:length(levels(data$Type))){
      
      datalist[[z]] <- data[which(data[,1] == levels(data$Type)[z]),]
      
      datalist[[z]] <- datalist[[z]][sample(1:nrow(datalist[[z]]), nrow(datalist[[z]])),]
      
      rownames(datalist[[z]]) <- 1:nrow(datalist[[z]])
      
    }
    
    totauc <- rep(0, k)
    
    totfm <- rep(0, k)
    
    totbr <- rep(0, k)
    
    samplist <- list()
    
    trainlist <- list()
    
    testlist <- list()
    
    size <- c()
    
    if (train == TRUE){
      
      for (z in 1:length(levels(data$Type))){
        
        trainlist[[z]] <- datalist[[z]]
        
        size <- c(size, length(trainlist[[z]][,1]))
        
      }
      
      minsize <- min(size)
      
      maxsize <- max(size)
      
      ratsize <- ceiling(maxsize / minsize)
      
      for (z in 1:length(levels(data$Type))){
        
        trainlist[[z]] <- trainlist[[z]][rep(seq.int(1, nrow(trainlist[[z]])), ratsize),][1:maxsize,]
        
        rownames(trainlist[[z]]) <- 1:maxsize
        
      }
      
      traindata <- trainlist[[1]]
      
      for (z in 2:length(levels(data$Type))){
        
        traindata <- rbind(traindata, trainlist[[z]])
        
      }
      
      if (alg == "SVM"){
        
        if (quiet == FALSE){
          
          print("Pretuning SVM...")
          
        }
        
        obj <- tune(svm, train.x = traindata[,-1], train.y = traindata[,1], ranges = list(cost = 2^(0:5)), tune.control(sampling = "fix"), kernel = "linear")
        
        bcost <- as.numeric(obj$best.parameters[1])
        
        if (quiet == FALSE){
          
          print(paste("Tuning Complete. Best Cost = ", bcost, ".", sep = ""))
          
        }
        
      }
      
      if (alg == "RF"){
        
        model <- randomForest(Type ~ ., data = traindata, mtry = mtr, ntree = ntree, nodesize = nodesize, importance = T, proximity = T)
        
        write.table(model$importance, file = "F:/Documents/PhD/ImpTest/importance_train.csv", sep=",")
        
      }
      else if (alg == "SVM"){
        
        model <- svm(Type ~ ., data = traindata, kernel = "linear", cost = bcost, gamma = 1, probability = TRUE)
        
      }
      else if (alg == "NN"){
        
        model <- mlptrain(traindata, testdata, hidden = hid)
        
      }
      else{
        
        stop(paste("'", alg, "' is an unacceptable algorithm. Exiting...", sep = ""))
        
      }
      
    }
    else{
      
      for (i in 1:k){
        
        if (featsel == FALSE){
          
          print(paste("Executing Cross-Validation ", i, "...", sep=""))
          
        }
        
        for (z in 1:length(levels(data$Type))){
          
          samplist[[z]] <- (1:nrow(datalist[[z]]) %% k != (i-1))
          
        }
        
        for (z in 1:length(levels(data$Type))){
          
          trainlist[[z]] <- datalist[[z]][samplist[[z]],]
          
          testlist[[z]] <- datalist[[z]][samplist[[z]] == F,]
          
          size <- c(size, length(trainlist[[z]][,1]))
          
        }
        
        minsize <- min(size)
        
        maxsize <- max(size)
        
        ratsize <- ceiling(maxsize / minsize)
        
        for (z in 1:length(levels(data$Type))){
          
          trainlist[[z]] <- trainlist[[z]][rep(seq.int(1, nrow(trainlist[[z]])), ratsize),][1:maxsize,]
          
          rownames(trainlist[[z]]) <- 1:maxsize
          
        }
        
        traindata <- trainlist[[1]]
        
        testdata <- testlist[[1]]
        
        for (z in 2:length(levels(data$Type))){
          
          traindata <- rbind(traindata, trainlist[[z]])
          
          testdata <- rbind(testdata, testlist[[z]])
          
        }
        
        ground <- matrix(0, nrow = nrow(testdata), ncol = length(levels(data$Type)))
        
        for (z in 1:nrow(testdata)){
          
          ground[z,which(levels(data$Type) == testdata[z,1])] <- 1
          
        }
        
        if (alg == "SVM" & i == 1 & j == 1){
          
          if (quiet == FALSE){
            
            print("Pretuning SVM and using values for all models (This seems reasonable in testing)...")
            
          }
          
          obj <- tune(svm, train.x = traindata[,-1], train.y = traindata[,1], ranges = list(cost = 2^(0:5)), tune.control(sampling = "fix"), kernel = "linear")
          
          bcost <- as.numeric(obj$best.parameters[1])
          
          if (quiet == FALSE){
            
            print(paste("Tuning Complete. Best Cost = ", bcost, ".", sep = ""))
            
          }
          
        }
        
        if (alg == "RF"){
          
          model <- randomForest(Type ~ ., data = traindata, mtry = mtr, ntree = ntree, nodesize = nodesize, importance = T, proximity = T)
          
          pred <- predict(model, testdata[,-1], type = "prob")
          
          totbr[i] <- (1/nrow(testdata)) * sum(sum((ground - pred)^2.0))
          
          fm <- try(as.numeric(as.vector(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1]))$byClass[,7])), TRUE)
          
          if (class(fm) == "try-error"){
            
            fm <- as.numeric(as.vector(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1]))$byClass[7]))
            
          }
          
          fm[which(is.na(fm) | is.nan(fm))] <- 0
          
          totfm[i] <- mean(fm)
          
          crosspredi <- cbind(testdata[,1], pred)
          
          crosspred <- rbind(crosspred, crosspredi)
          
          if (quiet == FALSE){
            
            print(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1])))
            
            print(paste("F1-Score for this Crossfold: ", totfm[i], sep=""))
            
            print(paste("Brier-Score for this Crossfold: ", totbr[i], sep=""))
            
          }
          
          roctest <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet)
          
        }
        else if (alg == "SVM"){
          
          model <- svm(Type ~ ., data = traindata, kernel = "linear", cost = bcost, gamma = 1, probability = TRUE)
          
          pred <- predict(model, testdata[,-1], probability = TRUE)
          
          pred <- attr(pred, "probabilities")
          
          roctest <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet)
          
        }
        else if (alg == "NN"){
          
          model <- mlptrain(traindata, testdata, hidden = hid)
          
          pred <- t(predict(model, data.matrix(testdata[,-1])))
          
          totbr[i] <- (1/nrow(testdata)) * sum(sum((ground - pred)^2.0))
          
          predmax <- max.col(pred) - 1
          
          print(confusionMatrix(predmax, as.numeric(testdata[,1])))
          
          fm <- as.numeric(as.vector(confusionMatrix(predmax, as.numeric(testdata[,1]))$byClass[,7]))
          
          fm[which(is.na(fm) | is.nan(fm))] <- 0
          
          totfm[i] <- mean(fm)
          
          roctest <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet)
          
        }
        else{
          
          stop(paste("'", alg, "' is an unacceptable algorithm. Exiting...", sep = ""))
          
        }
        
        totauc[i] <- as.numeric(roctest$auc)/100
        
      }
      
      finauc <- mean(totauc)
      
      finfm <- mean(totfm)
      
      finbr <- mean(totbr)
      
      if (quiet == FALSE){
        
        print(paste("The mean CV AUC value for ", k, "-crossfolds on repeat ", j, " is: ", finauc, ".", sep = ""))
        
        print(paste("The mean CV F1 value for ", k, "-crossfolds on repeat ", j, " is: ", finfm, ".", sep = ""))
        
        print(paste("The mean CV Brier-Score for ", k, "-crossfolds on repeat ", j, " is: ", finbr, ".", sep = ""))
        
      }
      
      ultauc[j,] <- totauc
      
      ultfm[j,] <- totfm
      
      ultbr[j,] <- totbr
      
    }
    
    mnauc <- mean(ultauc)
    
    sdauc <- sd(ultauc)
    
    mnfm <- mean(ultfm)
    
    sdfm <- sd(ultfm)
    
    mnbr <- mean(ultbr)
    
    sdbr <- sd(ultbr)
    
    if (quiet == FALSE){
      
      print(paste("The mean value of the AUC of the ", r, " repeats is ", mnauc, ".", sep = ""))
      
      print(paste("The standard deviation of the AUC of the ", r, " repeats is ", sdauc, ".", sep = ""))
      
      print(paste("The mean value of the F1 of the ", r, " repeats is ", mnfm, ".", sep = ""))
      
      print(paste("The standard deviation of the F1 of the ", r, " repeats is ", sdfm, ".", sep = ""))
      
      print(paste("The mean value of the Brier Score of the ", r, " repeats is ", mnbr, ".", sep = ""))
      
      print(paste("The standard deviation of the Brier Score of the ", r, " repeats is ", sdbr, ".", sep = ""))
      
    }
    
  }
  
  if (cvpred == TRUE){
    
    ultauc <- crosspred
    
  }
  
  if (train == FALSE){
    
    if (output == "AUC"){
      
      ultauc
      
    }
    else if (output == "F1"){
      
      ultfm
      
    }
    else if (output == "BS"){
      
      ultbr
      
    }
    else{
      
      print(paste("'", output, "' is an invalid output. Defaulting to AUC...", sep=""))
      
      ultauc
      
    }
    
  }
  else{
    
    model
    
  }
  
}