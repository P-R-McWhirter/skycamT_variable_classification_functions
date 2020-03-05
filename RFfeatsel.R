RFfeatsel <- function(data, k = 5, r = 1, seedno = 100, mtry = 4, coff = 5){
  
  library(caret)
  
  library(pROC)
  
  library(ROCR)
  
  library(RODBC)
  
  library(e1071)
  
  library(randomForest)
  
  start <- Sys.time()
  
  print(paste("Commencing Random Forest ", k, "-fold with ", r, " repeat(s) cross-validated feature selection...", sep = ""))
  
  classes <- subset(data, select = Type)
  
  features <- subset(data, select = -Type)
  
  data <- cbind(classes, features)
  
  finauc <- rep(0, (ncol(data)-coff))
  
  print("Evaluating full feature set model (iteration 1)...")
  
  runauc <- crossval(data, k = k, r = r, seed = seedno, alg = "RF", mtr = min(coff, mtry), quiet = T, featsel = T)
  
  finauc[1] <- mean(runauc)
  
  crossvaldata <- data
  
  rmfeat <- c("None")
  
  for (i in 2:(ncol(data)-coff)){
    
    print(paste("Commencing evaluation iteration ", i, "...", sep = ""))
    
    loopauc <- rep(0, (ncol(crossvaldata)))
    
    for (j in 2:(ncol(crossvaldata))){
      
      runauc <- crossval(crossvaldata[,-j], k = k, r = r, seed = seedno, alg = "RF", mtr = min(coff, mtry), quiet = T, featsel = T)
      
      loopauc[j-1] <- mean(runauc)
      
    }
    
    finauc[i] <- max(loopauc)
    
    dropfeat <- which.max(loopauc) + 1
    
    print(paste("Poorest feature determined to be: '", colnames(crossvaldata)[dropfeat], "'. Feature dropped for next loop...", sep = ""))
    
    rmfeat <- c(rmfeat, colnames(crossvaldata)[dropfeat])
    
    crossvaldata <- crossvaldata[,-dropfeat]
    
  }
  
  out <- list(finauc = finauc, rmfeat = rmfeat)
  
  print(Sys.time() - start)
  
  out
  
}