stringlength <- function(x, from = 0.01, to = 1000, val = 0.01) {
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  
  freq <- seq(from, to, by = val)
  n <- length(freq)
  
  LK <- rep.int(0, n)
  print(n)
  
  for (j in 1:n) {
    
    x <- xb
    
    p <- (((t - t[which.max(t)]) / freq[j]) + 0.25) %% 1
    x[, 1] <- p
    
    x <- x[order(x[,1]),]
    
    x1 <- x
    
    for (i in 2:length(x[,1])) {
      
      x1[i,] <- x[i-1,]
      
    }
    
    x1[1,] <- x[length(x[,1]),]
    
    NORM_CONSTANT <- 2
    
    COLUMN_VALUE <- 1
    
    normalization_constant <- ((length(x[,COLUMN_VALUE]) - 1) / (NORM_CONSTANT * length(x[,COLUMN_VALUE])))
    
    dispersion <- (x1[,COLUMN_VALUE] - x[,COLUMN_VALUE])^2.0
    
    variance <- (x[,COLUMN_VALUE] - mean(x[,COLUMN_VALUE]))^2.0
    
    dispersion2 <- (x1[,2] - x[,2])^2.0
    
    variance2 <- (x[,2] - mean(x[,2]))^2.0
    
    LK[j] = ((sum(dispersion + dispersion2)^0.5 / sum(variance + variance2)^0.5 * normalization_constant))
    
    x1 <- rep.int(0, length(x))
    
    x <- rep.int(0, length(x))
    
  }
  
  print(freq[which.min(LK)])
  
  plot(freq, LK, pch=19, type="n", main="SLLK", xlab = "Period (days)", ylab = "Log Lafler-Kinman value", ylim = c(max(LK), min(LK)))
  lines(freq, LK, type="l")
}