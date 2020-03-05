SLLK2 <- function(x, from = 0.001, to = 20, ofac = 1, plot = TRUE) {
  
  type = "frequency"
  alpha = 0.01
  level = 0
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  
  n <- length(y)
  tspan <- t[n] - t[1]
  step <- 1/(tspan * ofac)
  
  freq <- seq(from, to, by = step)
  
  n = length(freq)
  LK <- rep.int(0, n)
  
  x[,2] <- x[,2] - mean(x[,2])
  
  me = mean(x[,2])
  
  print(n)
  
  for (j in 1:n) {
    
    if (j %% 1000 == 0){
      
      print(j)
      
    }
    
    x <- xb
    
    x1 <- xb
    
    p <- (t / (1 / freq[j])) %% 1
    
    x[, 1] <- p
    
    x <- x[order(p),]
    
    #x[,2] <- (x[,2] - min(x[,2])) / (2.0 * (max(x[,2]) - min(x[,2]))) - 0.25
    
    x1 <- x
    
    for (i in 2:length(x[,1])) {
      
      x1[i,] <- x[i-1,]
      
    }
    
    x1[1,] <- x[length(x[,1]),]
    
    #LK[j] = sum(((x1[,2] - x[,2])^2.0 + (x1[,1] - x[,1])^2.0)^0.5) / sum(((x[,2] - me)^2.0))
    
    LK[j] = sum((x1[,2] - x[,2])^2) / sum((x[,2] - me)^2)
    
    LK[j] = LK[j] * ((length(x[,2]) - 1) / (2 * length(x[,2])))
    
  }
  
  #print(freq[which.max(LK)])
  
  #print(max(LK))
  
  peak.at = c(freq[which.min(LK)], 1/freq[which.min(LK)])
  
  p = 0
  
  if (plot == TRUE){
    
    plot(freq, LK, pch=19, type="n", main="SLLK", xlab = "Frequency (cycles/days)", ylab = "Normalized String Length", ylim = c(min(LK), max(LK)))
    
    lines(freq, LK, type="l")
    
  }
  
  sp.out <- list(scanned = freq, power = (max(LK) - LK), data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(max(LK) - LK), peak.at = peak.at, p.value = p)
}