spurper <- function(x, from = 0.001, to = 20, ofac = 1, name = "", plot = TRUE) {
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  
  fmin <- from
  fmax <- to
  
  n <- length(t)
  
  ttot <- max(t) - min(t)
  
  df <- 1 / (ttot * ofac)
  
  f <- seq(fmin, fmax, by = df)
  
  m <- length(f)
  
  p <- rep(0, m)
  
  print("Calculating Spurious Frequencies...")
  
  for (j in 1:m){
    
    p[j] <- abs(sum(exp(1i * 2.0 * pi * f[j] * t)))^2 / n
    
  }
  
  peak.at = c(f[which.max(p)], 1/f[which.max(p)])
  
  if (plot == TRUE){
    
    plot(f, p, pch=19, type="n", main=paste("Spurious Frequencies for ", name, "", sep=""), xlab = "Frequency (Cycles/days)", ylab = "Spectral Power", ylim = c(min(p), max(p)))
    
    lines(f, p, type="l")
    
  }
  
  sp.out <- list(scanned = f, power = p, data = datanames, n = n, ofac = ofac, peak = max(p), peak.at = peak.at)
  
}








spurperlsp <- function(x, from = 0.001, to = 20, ofac = 1, name = "", plot = TRUE) {
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  
  fmin <- from
  fmax <- to
  
  n <- length(t)
  
  ttot <- max(t) - min(t)
  
  df <- 1 / (ttot * ofac)
  
  f <- seq(fmin, fmax, by = df)
  
  m <- length(f)
  
  p <- rep(0, m)
  
  print("Calculating Spurious Frequencies...")
  
  norm <- 1/(2 * var(x[,2]))
  
  for (j in 1:m){
    
    wi <- 2 * pi * f[j]
    
    tau <- 0.5 * atan2(sum(sin(wi * x[,1])), sum(cos(wi * x[,1])))/wi
    
    arg <- wi * (x[,1] - tau)
    
    cs <- cos(arg)
    
    sn <- sin(arg)
    
    A <- (sum(cs))^2
    
    B <- sum(cs * cs)
    
    C <- (sum(sn))^2
    
    D <- sum(sn * sn)
    
    PN <- A/B + C/D
    
    p[j] <- PN
    
  }
  
  p <- norm * p
  
  peak.at = c(f[which.max(p)], 1/f[which.max(p)])
  
  if (plot == TRUE){
    
    plot(f, p, pch=19, type="n", main=paste("Spurious Frequencies for ", name, "", sep=""), xlab = "Frequency (Cycles/days)", ylab = "Spectral Power", ylim = c(min(p), max(p)))
    
    lines(f, p, type="l")
    
  }
  
  sp.out <- list(scanned = f, power = p, data = datanames, n = n, ofac = ofac, peak = max(p), peak.at = peak.at)
  
}