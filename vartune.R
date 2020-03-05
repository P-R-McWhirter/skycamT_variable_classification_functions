vartune <- function(x, from = 0.001, to = 20, ofac = 1, plot = TRUE) {
  
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
  
  #freq <- seq(from, to, by = step)
  
  #n = length(freq)
  
  #print(n)
    
  freq <- seq(from, to, length.out = 100)
    
  n <- 100
  
  #print(n)
  
  pwr <- rep.int(0, n)
  
  x[,2] <- x[,2] - mean(x[,2])
  
  me = mean(x[,2])
  
  ts <- x
  
  tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(mean(ts[,2]), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, ts)
  
  tsin <- tsin[order(tsin[,1]),]
  
  for (k in 1:n) {
    
    model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*freq[k]*tsin[,1]) + cos(2*pi*freq[k]*tsin[,1]) +
                  sin(4*pi*freq[k]*tsin[,1]) + cos(4*pi*freq[k]*tsin[,1]) + sin(6*pi*freq[k]*tsin[,1]) +
                  cos(6*pi*freq[k]*tsin[,1]) + sin(8*pi*freq[k]*tsin[,1]) + cos(8*pi*freq[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    coeff <- model$coefficients
    
    coeff[is.na(coeff)] <- 0
    
    yp <- y - (coeff[2]*t + coeff[3]*sin(2*pi*freq[k]*t) + coeff[4]*cos(2*pi*freq[k]*t) +
                 coeff[5]*sin(4*pi*freq[k]*t) + coeff[6]*cos(4*pi*freq[k]*t) + coeff[7]*sin(6*pi*freq[k]*t) +
                 coeff[8]*cos(6*pi*freq[k]*t) + coeff[9]*sin(8*pi*freq[k]*t) + coeff[10]*cos(8*pi*freq[k]*t))
    
    pwr[k] <- var(yp) / var(y)
    
  }
  
  pwr <- 1 - pwr
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  if (plot == TRUE){
    
    plot(freq, pwr, pch=19, type="n", main="VRP", xlab = "Frequency (cycles/days)", ylab = "Variance Ratio", ylim = c(min(pwr), max(pwr)))
    
    lines(freq, pwr, type="l")
    
  }
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}