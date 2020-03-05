bigparamsel <- function(input, optfn = "AUC", mtrs = seq(5, 50, by = 5), ntrees = c(500, 1000, 1500), nodesizes = c(25, 50, 75), oset = 0, quiet = T){
  
  result <- list()
  
  cutoffs <- list()
  
  #mtrs <- 10
  
  #ntrees <- 250
  
  #nodesizes <- 18
  
  config <- data.frame("", "", "", stringsAsFactors = F)
  
  colnames(config) <- c("mtrs", "ntrees", "nodesizes")
  
  z <- 1
  
  for (i in 1:length(mtrs)){
    
    for (j in 1:length(ntrees)){
      
      for (k in 1:length(nodesizes)){
        
        print(paste("Running set ", z, "... mtr = ", mtrs[i], " ntree = ", ntrees[j], " nodesize = ", nodesizes[k], ".", sep = ""))
        
        cross <- crossval(input, k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = mtrs[i], ntree = ntrees[j], nodesize = nodesizes[k], quiet = quiet, featsel = T, output = optfn, oset = oset)
        
        result[[z]] <- cross$result
        
        print(mean(cross$result))
        
        cutoffs[[z]] <- cross$cutoffs
        
        config[z,] <- c(mtrs[i], ntrees[j], nodesizes[k])
        
        z <- z + 1
        
      }
      
    }
    
  }
  
  out <- list(result = result, config = config, cutoffs = cutoffs)
  
}


bigparamselmean <- function(input){
  
  n <- length(input)
  
  result <- rep(0, n)
  
  result2 <- rep(0, n)
  
  for (i in 1:n){
    
    result[i] <- mean(input[[i]])
    
    result2[i] <- sd(input[[i]])
    
  }
  
  output <- list(meanvals = result, sdvals = result2)
  
  output
  
}




plotbigparamsel <- function(input){
  
  plot(2:7, input[c(2, 20, 38, 56, 74, 92)], pch = 19, type = "b", lwd = 4, lty = 1, ylim = c(min(input), max(input)), main = "Random Forest tuning parameters", ylab = "Mean AUC", xlab = "mtry")
  
  lines(2:7, input[c(4, 22, 40, 58, 76, 94)], pch = 18, type = "b", lwd = 4, col = "red", lty = 2)
  
  lines(2:7, input[c(6, 24, 42, 60, 78, 96)], pch = 17, type = "b", lwd = 4, col = "blue", lty = 3)
  
  lines(2:7, input[c(8, 26, 44, 62, 80, 98)], pch = 15, type = "b", lwd = 4, col = "darkgreen", lty = 4)
  
  lines(2:7, input[c(10, 28, 46, 64, 82, 98)], pch = 20, type = "b", lwd = 4, col = "darkred", lty = 5)
  
  lines(2:7, input[c(12, 30, 48, 66, 84, 100)], pch = 3, type = "b", lwd = 4, col = "darkblue", lty = 6)
  
  lines(2:7, input[c(14, 32, 50, 68, 86, 102)], pch = 0, type = "b", lwd = 4, col = "cyan4", lty = 7)
  
  lines(2:7, input[c(16, 34, 52, 70, 88, 106)], pch = 2, type = "b", lwd = 4, col = "purple", lty = 8)
  
  lines(2:7, input[c(18, 36, 54, 72, 90, 108)], pch = 5, type = "b", lwd = 4, col = "goldenrod4", lty = 9)
  
  leg <- c("ntree = 500, nodesize = 3", "ntree = 500, nodesize = 5", "ntree = 500, nodesize = 7",
           "ntree = 1000, nodesize = 3", "ntree = 1000, nodesize = 5", "ntree = 1000, nodesize = 7",
           "ntree = 2000, nodesize = 3", "ntree = 2000, nodesize = 5", "ntree = 2000, nodesize = 7")
  
  colset <- c("black", "red", "blue", "darkgreen", "darkred", "darkblue", "cyan4", "purple", "goldenrod4")
  
  pchset <- c(19, 18, 17, 15, 20, 3, 0, 2, 5)
  
  legend("topright", inset = .05, cex = 1, title = "RF Parameters", leg, lty = c(1,2,3,4,5,6,7,8,9), lwd = c(4,4,4,4,4,4,4,4,4), col = colset, pch = pchset, bg = "grey96")
  
}