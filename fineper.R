fineper <- function(ts, f, type) {
  
  print("Fine tuning top 10 periods...")
  
  f <- f[1:10]
  
  if (any(is.na(f))){
    
    f <- f[!is.na(f)]
    
    if (length(f) < 1){
      
      f <- 1 / (max(ts[,1]) - min(ts[,1]))
      
    }
    
  }
  
  trialpwr <- rep(0, length(f))
  
  maxtime <- max(ts[,1]) - min(ts[,1])
  
  if (type == "LSP"){
  
    for (k in 1:length(f)){
    
      fineans <- lsp(ts[,1:2], from = (0.9975 * f[k]), to = (1.0025 * f[k]), type = "frequency", ofac = 100, plot = FALSE)
    
      ftrial <- fineans$scanned[order(fineans$power, decreasing = T)]
    
      f[k] <- ftrial[1]
    
      trialpwr[k] <- max(fineans$power)
    
    }
  
  }
  else{
    
    f <- f[1]
    
    trialpwr <- trialpwr[1]
    
    fineans <- lsp(ts[,1:2], from = (0.9975 * f[1]), to = (1.0025 * f[1]), type = "frequency", ofac = 100, plot = FALSE)
    
    ftrial <- fineans$scanned[order(fineans$power, decreasing = T)]
    
    f[1] <- ftrial[1]
    
    trialpwr[1] <- max(fineans$power)
    
  }
  
  f <- f[1]
  
  pwr <- max(trialpwr)
  
  df <- 0.01 / maxtime
  
  if (!is.na(f)){
    
    scanlen <- length(seq((0.9975 * f), (1.0025 * f), by = df))
    
    pwr <- max(trialpwr)
    
  }
  else{
    
    scanlen <- 1
    
    pwr <- 0
    
    f <- 1 / maxtime
    
  }
  
  sp.out <- list(freq = f, power = pwr, len = scanlen)
  
}