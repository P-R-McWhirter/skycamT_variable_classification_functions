lcebchk <- function(ts, per = 1.0, ep = 0.01, bins = 100, quiet = FALSE) {
  
  library(e1071)
  
  start <- Sys.time()
  
  lowper <- per * (1 - ep)
  
  highper <- per * (1 + ep)
  
  pers <- seq(lowper, highper, length.out = 1000)
  
  perf <- rep(0, length(pers))
  
  res <- matrix(0, nrow = length(pers), ncol = 3)
  
  ts[,2] <- ts[,2] - median(ts[,2])
  
  m <- median(ts)
  
  ts <- ts[which(ts[,2] > 0),]
  
  for (i in 1:length(pers)){
    
    fts <- ftsb
    
    fts1 <- ftsb
    
    p <- (t / (1 / freq[k])) %% 1
    
    x[, 1] <- p
    
    x <- x[order(p),]
    
    #x[,2] <- (x[,2] - min(x[,2])) / (2.0 * (max(x[,2]) - min(x[,2]))) - 0.25
    
    x1 <- x
    
    for (i in 2:length(x[,1])) {
      
      x1[i,] <- x[i-1,]
      
    }
    
    x1[1,] <- x[length(x[,1]),]
    
    #LK[j] <- sum(((x1[,2] - x[,2])^2.0 + (x1[,1] - x[,1])^2.0)^0.5) / sum(((x[,2] - me)^2.0))
    
    pwr[k] <- sum((x1[,2] - x[,2])^2) / sum((x[,2] - me)^2)
    
  }
  
  plot(pers, perf, type = "l")
  
  print(res[which.max(perf),])
  
  plot((((((ts[,1] - ts[which.max(ts[,2]),1])/ pers[which.max(perf)]) + 0.5) %% 1) - 0.5), ts[,2], pch=19)
  
  #lines(seq(-0.5, 0.5, length.out = 10000), res[which.max(perf),1] * exp(-0.5 * ((seq(-0.5, 0.5, length.out = 10000) - res[which.min(perf),2])/res[which.min(perf),3])^2), col = "red")
  
  print(Sys.time() - start)
  
  pers[which.min(perf)]
  
}






lcebchk_old <- function(ts, per = 1.0, ep = 0.01, quiet = FALSE) {
  
  start <- Sys.time()
  
  lowper <- per * (1 - ep)
  
  highper <- per * (1 + ep)
  
  pers <- seq(lowper, highper, length.out = 1000)
  
  perf <- rep(0, length(pers))
  
  res <- matrix(0, nrow = length(pers), ncol = 3)
  
  ts[,2] <- ts[,2] - median(ts[,2])
  
  ts <- ts[which(ts[,2] > 0),]
  
  for (i in 1:length(pers)){
    
    fts <- ts
    
    fts[,1] <- ((((fts[,1] - fts[which.min(fts[,2]),1])/ pers[i]) + 0.5) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    fts <- fts[which(fts[,1] > -0.1 & fts[,1] < 0.1),]
    
    r <- fts[,2]
    
    f <- function(par)
    {
      m <- par[1]
      sd <- par[2]
      k <- par[3]
      rhat <- k * exp(-0.5 * ((fts[,1] - m)/sd)^2)
      sum((r - rhat)^2)
    }
    
    
    
    oldcurr <- 100000000
    
    curr <- 1000000
    
    #plot(fts, pch=19)
    
    for (k in 1:100){
      
      parm <- rep(0, 3)
      
      parm[1] <- runif(1, 0.0, 2)
      
      parm[2] <- runif(1, -0.5, 0.5)
      
      parm[3] <- runif(1, 0.001, 0.05)
      
      SSTlm <- lm(fts[,2] ~ -1 + exp(-0.5 * ((fts[,1])/0.05)^2), data = as.data.frame(fts), weights = rep(1, nrow(fts)))
      
      #print(summary(SSTlm))
      
      co <- SSTlm$coefficients
      
      co[is.na(co)] <- 0
      
      curr <- sum((fts[,2] - (co[1] + exp(-0.5 * ((fts[,1] - 0.017)/0.02)^2)))^2/rep(1, nrow(fts))^2)
      
      #print(curr)
      
      if (curr < oldcurr){
        
        oldcurr <- curr
        
        res[i,] <- c(co[1], 0.017, 0.02)
        
      }
      
    }
    
    perf[i] <- curr
    
    #print(curr)
    
    #lines(seq(-0.5, 0.5, length.out = 10000), parm[1] * exp(-0.5 * ((seq(-0.5, 0.5, length.out = 10000) - parm[2])/parm[3])^2), col = "red")
    
  }
  
  plot(pers, perf, type = "l")
  
  print(res[which.min(perf),])
  
  plot((((((ts[,1] - ts[which.max(ts[,2]),1])/ pers[which.min(perf)]) + 0.5) %% 1) - 0.5), ts[,2], pch=19)
  
  lines(seq(-0.5, 0.5, length.out = 10000), res[which.min(perf),1] * exp(-0.5 * ((seq(-0.5, 0.5, length.out = 10000) - res[which.min(perf),2])/res[which.min(perf),3])^2), col = "red")
  
  pers[which.min(perf)]
  
}