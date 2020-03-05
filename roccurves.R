roccurves <- function(x, a, samp, type, oset = 0, quiet = FALSE) {
  
  library(ROCR)
  
  library(pROC)
  
  x <- as.matrix(x)
  
  if (is.null(samp)){
    
    lev <- as.character(as.vector(levels(type)))
    
  }
  else{
    
    lev <- as.character(as.vector(levels(type[-samp])))
    
  }
  
  nocl <- as.numeric(dim(x)[2])
  
  maxlim <- as.numeric(nocl + 10)
  
  auc <- rep(0, nocl)
  
  cutoff <- rep(0, nocl)
  
  for (i in 1:nocl){
    
    actual <- a
    
    actual[which(actual == i)] <- maxlim
    
    actual[which(actual != maxlim)] <- 0
    
    actual[which(actual == maxlim)] <- 1
    
    if (is.null(samp)){
      
      pred <- prediction(x[,i], actual)
      
    }
    else{
      
      pred <- prediction(x[,i], actual[-samp])
      
    }
    
    perf <- performance(pred, measure = "tpr", x.measure = "fpr")
    
    auc.perf <- performance(pred, measure = "auc")
    
    auc[i] <- (auc.perf@y.values[[1]])
    
    ss <- performance(pred, measure = "sens", x.measure = "spec")
    
    iuc <- abs(ss@y.values[[1]] + oset - auc[i]) + abs(ss@x.values[[1]] - oset - auc[i])
    
    cutoff[i] <- ss@alpha.values[[1]][which.min(iuc)]
    
    if (quiet == FALSE){
      
      if (i == 1){
        
        plot.new()
        
        frame()
        
        plot(perf, col = rainbow(nocl)[i], main = "One-vs-many classification accuracy ROC curve")
        
        legend("bottomright", inset = .05, cex = 1, title = "Classes", lev, lty = c(1,1), lwd = c(2,2), col = rainbow(nocl), bg = "grey96")
        
        par(new=TRUE)
        
      }
      else{
        
        plot(perf, col = rainbow(nocl)[i], main = "", sub = "", xlab = "", ylab = "", xaxt = "n", yaxt = "n", axes = FALSE, ann = FALSE, col.lab = "black")
        
        par(new=TRUE)
        
      }
      
    }
    
  }
  
  if (quiet == FALSE){
    
    abline(a=0, b=1, lty=2, lwd=3)
    
    print("Printing results...")
    
    for (i in 1:nocl){
      
      print(paste("The AUC for Class ", i, ", named '", lev[i], "' is: ", auc[i], ".", sep=""))
      
    }
    
  }
  
  prediction_class <- apply(x, 1, which.max)
  
  roctest <- multiclass.roc(predictor=prediction_class, response=as.numeric(a),percent=T)
  
  if (quiet == FALSE){
    
    print("Calculating Multi-class area under the curve...")
    
    print(roctest$auc)
    
  }
  
  cutoff[which(cutoff < 0)] <- 0.000
  
  cutoff[which(cutoff > 1)] <- 1.000
  
  result <- list(roctest = roctest, auc = auc, cutoff = cutoff)
  
  result
  
}