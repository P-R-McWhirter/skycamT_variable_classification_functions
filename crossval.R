crossval <- function(data, k = 10, r = 1, seed = 100, output = "AUC", alg = "RF", scaleit = FALSE, hid = 200, mtr = round(sqrt(ncol(data)-1)), ntree = 1000, nodesize = 2, oset = 0, quiet = FALSE, featsel = FALSE, train = FALSE, cvpred = FALSE){
  
  #library(mxnet)
  
  library(caret)
  
  library(pROC)
  
  library(ROCR)
  
  library(RODBC)
  
  library(e1071)
  
  library(randomForest)
  
  library(gbm)
  
  data[,1] <- as.factor(data[,1])
  
  data[,1] <- as.character(data[,1])
  
  data[,1] <- as.factor(data[,1])
  
  if (scaleit == TRUE){
    
    data[,-1] <- scale(data[,-1])
    
  }
  
  if (length(seed) != r){
    
    if (quiet == FALSE){
      
      print("Seed vector must be the same length as the number of repeats.")
      
      print(paste("Setting number of repeats to ", length(seed), ".", sep=""))
      
    }
    
    r <- length(seed)
    
  }
  else{
    
    if (quiet == FALSE){
      
      print(paste("Number of repeats is: ", r, ".", sep = ""))
      
    }
    
  }
  
  if (train == TRUE){
    
    print("Over-riding options as we are training a single model...")
    
    seed <- 100
    
    r <- 1
    
  }
  
  if (cvpred == TRUE){
    
    print("Over-riding options as we are calibrating probabilities using Brier Score...")
    
    seed <- 100
    
    r <- 1
    
  }
  
  ultauc <- matrix(0, ncol = k, nrow = r)
  
  ultfm <- matrix(0, ncol = k, nrow = r)
  
  ultbr <- matrix(0, ncol = k, nrow = r)
  
  fincuts <- list()
  
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
    
    cutoffs <- matrix(0, nrow = k, ncol = length(levels(data$Type)))
    
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
        
        model <- svm(Type ~ ., data = traindata, kernel = "polynomial", cost = bcost, gamma = bgamma, degree = bdegree, probability = TRUE)
        
      }
      else if (alg == "NN"){
        
        model <- neuralnet(traindata, testdata, hidden = hid)
        
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
          
          obj <- tune(svm, train.x = traindata[,-1], train.y = traindata[,1], ranges = list(cost = 2^(-2:5), gamma = 2^(-8:0), degree = 1:3), tune.control(sampling = "fix"), kernel = "polynomial")
          
          bcost <- as.numeric(obj$best.parameters[1])
          
          bgamma <- as.numeric(obj$best.parameters[2])
          
          bdegree <- as.numeric(obj$best.parameters[3])
          
          if (quiet == FALSE){
            
            print(paste("Tuning Complete. Best Cost = ", bcost, ". Best Gamma = ", bgamma, ". Best Degree = ", bdegree, ".", sep = ""))
            
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
          
          rocset <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet, oset = oset)
          
          roctest <- rocset$roctest
          
          cutoffs[i,] <- rocset$cutoff
          
        }
        else if (alg == "SVM"){
          
          model <- svm(Type ~ ., data = traindata, kernel = "polynomial", cost = bcost, gamma = bgamma, degree = bdegree, probability = TRUE)
          
          pred <- predict(model, testdata[,-1], probability = TRUE)
          
          pred <- attr(pred, "probabilities")
          
          totbr[i] <- (1/nrow(testdata)) * sum(sum((ground - pred)^2.0))
          
          fm <- try(as.numeric(as.vector(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1]))$byClass[,7])), TRUE)
          
          if (class(fm) == "try-error"){
            
            fm <- as.numeric(as.vector(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1]))$byClass[7]))
            
          }
          
          fm[which(is.na(fm) | is.nan(fm))] <- 0
          
          totfm[i] <- mean(fm)
          
          rocset <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet, oset = oset)
          
          roctest <- rocset$roctest
          
          print(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1])))
          
          cutoffs[i,] <- rocset$cutoff
          
        }
        else if (alg == "NN"){
          
          neuralout <- neuralnet(traindata, testdata, hid = hid, epochs = 1000)
          
          model <- neuralout$model
          
          totfm[i] <- mean(neuralout$fm)
          
          pred <- predict_proba(model, neuralout$x_test)
          
          totbr[i] <- (1/nrow(testdata)) * sum(sum((ground - pred)^2.0))
          
          rocset <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet, oset = oset)
          
          roctest <- rocset$roctest
          
          cutoffs[i,] <- rocset$cutoff
          
        }
        else if (alg == "DNN"){
          
          neuralout <- deepneuralnet(traindata, testdata, hid = hid, epochs = 500)
          
          model <- neuralout$model
          
          totfm[i] <- mean(neuralout$fm)
          
          pred <- predict_proba(model, neuralout$x_test)
          
          rocset <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet, oset = oset)
          
          roctest <- rocset$roctest
          
          cutoffs[i,] <- rocset$cutoff
          
        }
        else if (alg == "NB"){
          
          model <- naiveBayes(Type ~ ., data = traindata)
          
          pred <- predict(model, testdata[,-1], "raw")
          
          totbr[i] <- (1/nrow(testdata)) * sum(sum((ground - pred)^2.0))
          
          fm <- try(as.numeric(as.vector(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1]))$byClass[,7])), TRUE)
          
          if (class(fm) == "try-error"){
            
            fm <- as.numeric(as.vector(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1]))$byClass[7]))
            
          }
          
          fm[which(is.na(fm) | is.nan(fm))] <- 0
          
          totfm[i] <- mean(fm)
          
          rocset <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet, oset = oset)
          
          roctest <- rocset$roctest
          
          print(confusionMatrix(predict(model, testdata[,-1]), (testdata[,1])))
          
          cutoffs[i,] <- rocset$cutoff
          
        }
        else if (alg == "kNN"){
          
          pred <- knn3Train(traindata[,-1], testdata[,-1], traindata[,1], k = 10, prob = TRUE)
          
          #pred <- predict(model, testdata[,-1], "raw")
          
          pred <- attr(pred, "prob")
          
          print(dim(pred))
          
          totbr[i] <- (1/nrow(testdata)) * sum(sum((ground - pred)^2.0))
          
          fm <- try(as.numeric(as.vector(confusionMatrix(pred, (testdata[,1]))$byClass[,7])), TRUE)
          
          if (class(fm) == "try-error"){
            
            fm <- as.numeric(as.vector(confusionMatrix(pred, (testdata[,1]))$byClass[7]))
            
          }
          
          fm[which(is.na(fm) | is.nan(fm))] <- 0
          
          totfm[i] <- mean(fm)
          
          rocset <- roccurves(pred, as.numeric(testdata[,1]), NULL, testdata[,1], quiet = quiet, oset = oset)
          
          roctest <- rocset$roctest
          
          print(confusionMatrix(pred, (testdata[,1])))
          
          cutoffs[i,] <- rocset$cutoff
          
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
      
      fincuts[[j]] <- cutoffs
      
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
  
  if (cvpred == TRUE | train == TRUE){
    
    fincuts <- NULL
    
  }
  else{
    
    fincuts <- do.call(rbind, fincuts)
    
    fincuts <- colMeans(fincuts, na.rm = T)
    
  }
  
  if (train == FALSE){
    
    if (output == "AUC"){
      
      result <- list(result = ultauc, cutoffs = fincuts)
      
      result
      
    }
    else if (output == "F1"){
      
      result <- list(result = ultfm, cutoffs = fincuts)
      
      result
      
    }
    else if (output == "BS"){
      
      result <- list(result = ultbr, cutoffs = fincuts)
      
      result
      
    }
    else{
      
      print(paste("'", output, "' is an invalid output. Defaulting to AUC...", sep=""))
      
      result <- list(result = ultauc, cutoffs = fincuts)
      
      result
      
    }
    
  }
  else{
    
    model
    
  }
  
}