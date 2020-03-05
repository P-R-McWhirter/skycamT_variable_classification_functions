bootstrapping_error <- function(data, runs = 1000, ep = 0.01, seedno = 1){
  
  set.seed(seedno)
  
  resulterr <- list()
  
  out <- matrix(0, nrow = runs, ncol = 3)
    
  bb <- synthpermodes(data, ep = ep, quiet = T)
    
  hitsampmean <- sum(bb$Hit)/nrow(bb)
  
  multsampmean <- (sum(bb$Multiple) + sum(bb$Sub.Multiple))/nrow(bb)
  
  alsampmean <- (sum(bb$Alias.1) + sum(bb$Alias.0.5))/nrow(bb)
  
  unsampmean <- sum(bb$Unknown)/nrow(bb)
  
  print(c(hitsampmean, multsampmean, alsampmean, unsampmean))
  
  resultmean <- c(hitsampmean, multsampmean, alsampmean, unsampmean)
  
  n <- nrow(data)
  
  setdata <- matrix(0, nrow = nrow(data), ncol = 4)
  
  for (i in 1:nrow(bb)){
  
    setdata[i,1] <- bb$Hit[i]
  
    setdata[i,2] <- max(bb$Multiple[i], bb$Sub.Multiple[i])
  
    setdata[i,3] <- max(bb$Alias.1[i], bb$Alias.0.5[i])
    
    setdata[i,4] <- bb$Unknown[i]
  
  }
  
  hittrialset <- matrix(0, nrow = runs, ncol = nrow(setdata))
  
  multtrialset <- matrix(0, nrow = runs, ncol = nrow(setdata))
  
  aliastrialset <- matrix(0, nrow = runs, ncol = nrow(setdata))
  
  untrialset <- matrix(0, nrow = runs, ncol = nrow(setdata))
  
  for (i in 1:runs){
    
    trialset <- setdata[sample(1:nrow(setdata), nrow(setdata), replace = T),]
    
    hittrialset[i,] <- as.vector(trialset[,1])
    
    multtrialset[i,] <- as.vector(trialset[,2])
    
    aliastrialset[i,] <- as.vector(trialset[,3])
    
    untrialset[i,] <- as.vector(trialset[,4])
    
  }
  
  hittrialset <- rowMeans(hittrialset)
  
  multtrialset <- rowMeans(multtrialset)
  
  aliastrialset <- rowMeans(aliastrialset)
  
  untrialset <- rowMeans(untrialset)
  
  hittrialset <- sort(hittrialset - hitsampmean)
  
  multtrialset <- sort(multtrialset - multsampmean)
  
  aliastrialset <- sort(aliastrialset - alsampmean)
  
  untrialset <- sort(untrialset - unsampmean)
  
  hits <- quantile(hittrialset, c(0.05, 0.95))
  
  mults <- quantile(multtrialset, c(0.05, 0.95))
  
  aliass <- quantile(aliastrialset, c(0.05, 0.95))
  
  unknowns <- quantile(untrialset, c(0.05, 0.95))
  
  print(as.numeric(hits))
  
  print(as.numeric(mults))
  
  print(as.numeric(aliass))
  
  print(as.numeric(unknowns))
  
  resulterr[[1]] <- as.numeric(hits)
  
  resulterr[[2]] <- as.numeric(mults)
  
  resulterr[[3]] <- as.numeric(aliass)
  
  resulterr[[4]] <- as.numeric(unknowns)
  
  fin <- list(resultmean = resultmean, resulterr = resulterr)
  
  fin
  
}