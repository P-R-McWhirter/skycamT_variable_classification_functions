testpolyfit <- function(data, p, seed = 1:10, gen = T){
  
  if (gen == T){
    
    out <- matrix(0, nrow = length(seed), ncol = 128)
    
    out <- as.data.frame(out)
  
    for (i in 1:length(seed)){
    
      polyfit <- try(lcgbinpolyfit(data, spur = T, pero = p, delta = 1, seedno = seed[i]), TRUE)
    
      if (class(polyfit) == "try-error"){
    
        out[i,] <- rep(NA, 128)
      
      }
      else{
      
        out[i,] <- polyfit[9:136]
      
        colnames(out) <- colnames(polyfit)[9:136]
      
      }
    
    }
    
  }
  else{
    
    out <- matrix(0, nrow = length(seed), ncol = 20)
    
    out <- as.data.frame(out)
    
    for (i in 1:length(seed)){
      
      polyfit <- try(lcbinpolyfit(data, spur = T, pero = p, delta = 1, seedno = seed[i]), TRUE)
      
      if (class(polyfit) == "try-error"){
        
        out[i,] <- rep(NA, 20)
        
      }
      else{
        
        out[i,] <- polyfit[4:23]
        
        colnames(out) <- colnames(polyfit)[4:23]
        
      }
      
    }
    
  }
  
  out
  
}



bigtestpolyfit <- function(data, p, seed = 1:10, gen = T){
  
  start <- Sys.time()
  
  p <- as.numeric(as.vector(p))
  
  if (gen == T){
    
    testsetmean <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 128))
    
    testsetsd <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 128))
    
    results <- list()
  
    for (i in 1:nrow(data)){
  
      res <- testpolyfit(data[i,], p[i], seed = seed, gen = gen)
    
      results[[i]] <- res
    
      testsetmean[i,] <- colMeans(res, na.rm = T)
    
      for (j in 1:128){
      
        testsetsd[i,j] <- sd(res[,j], na.rm = T)
      
      }
  
    }
    
  }
  else{
    
    testsetmean <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 20))
    
    testsetsd <- as.data.frame(matrix(0, nrow = nrow(data), ncol = 20))
    
    results <- list()
    
    for (i in 1:nrow(data)){
      
      res <- testpolyfit(data[i,], p[i], seed = seed, gen = gen)
      
      results[[i]] <- res
      
      testsetmean[i,] <- try(colMeans(res, na.rm = T), TRUE)
      
      for (j in 1:20){
        
        testsetsd[i,j] <- try(sd(res[,j], na.rm = T), TRUE)
        
      }
      
    }
    
  }
  
  colnames(testsetmean) <- colnames(res)
  
  colnames(testsetsd) <- colnames(res)
  
  out <- list(results = results, testsetmean = testsetmean, testsetsd = testsetsd)
  
  print(Sys.time() - start)
  
  out
  
}