VRP <- function(x, filt, ofac = 1, from = 0.001, to = 20, nt = 150, b = 3, plot = TRUE) {
  
  type = "frequency"
  alpha = 0.01
  level = 0
  
  eps = 1e-3
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  yb <- y
  
  n <- length(y)
  tspan <- t[n] - t[1]
  
  df <- 1 / (tspan * ofac)
  
  print(paste("Starting trial period estimation. # of trials per bin: ", nt, " and # of bins: ", b, ".", sep=""))
  
  d <- rep(0, n-1)
  
  for (k in 1:(n-1)){
    
    d[k] <- (y[k+1] - y[k]) / (t[k+1] - t[k])
    
  }
  
  sortedmags <- sort(y)
  
  E1 <- sortedmags[ceiling(0.1 * n)]
  
  E2 <- sortedmags[ceiling(0.2 * n)]
  
  E3 <- sortedmags[ceiling(0.3 * n)]
  
  E4 <- sortedmags[ceiling(0.4 * n)]
  
  E5 <- sortedmags[ceiling(0.5 * n)]
  
  E6 <- sortedmags[ceiling(0.6 * n)]
  
  E7 <- sortedmags[ceiling(0.7 * n)]
  
  E8 <- sortedmags[ceiling(0.8 * n)]
  
  E9 <- sortedmags[ceiling(0.9 * n)]
  
  Dj <- rep(0, 10)
  
  band <- rep(0, (n-1))
  
  for (i in 1:(n-1)){
    
    if ((y[i] >= min(y)) & (y[i] < E1)){
      
      Dj[1] <- Dj[1] + abs(d[k])
      
      band[i] <- 1
      
    }
    else if ((y[i] >= E1) & (y[i] < E2)){
      
      Dj[2] <- Dj[2] + abs(d[k])
      
      band[i] <- 2
      
    }
    else if ((y[i] >= E2) & (y[i] < E3)){
      
      Dj[3] <- Dj[3] + abs(d[k])
      
      band[i] <- 3
      
    }
    else if ((y[i] >= E3) & (y[i] < E4)){
      
      Dj[4] <- Dj[4] + abs(d[k])
      
      band[i] <- 4
      
    }
    else if ((y[i] >= E4) & (y[i] < E5)){
      
      Dj[5] <- Dj[5] + abs(d[k])
      
      band[i] <- 5
      
    }
    else if ((y[i] >= E5) & (y[i] < E6)){
      
      Dj[6] <- Dj[6] + abs(d[k])
      
      band[i] <- 6
      
    }
    else if ((y[i] >= E6) & (y[i] < E7)){
      
      Dj[7] <- Dj[7] + abs(d[k])
      
      band[i] <- 7
      
    }
    else if ((y[i] >= E7) & (y[i] < E8)){
      
      Dj[8] <- Dj[8] + abs(d[k])
      
      band[i] <- 8
      
    }
    else if ((y[i] >= E8) & (y[i] < E9)){
      
      Dj[9] <- Dj[9] + abs(d[k])
      
      band[i] <- 9
      
    }
    else if ((y[i] >= E9) & (y[i] <= max(y))){
      
      Dj[10] <- Dj[10] + abs(d[k])
      
      band[i] <- 10
      
    }
    
  }
  
  Dj <- as.matrix(cbind(Dj, c(1:10)))
  
  Dj <- Dj[order(Dj[,1], decreasing = T),]
  
  Dj <- Dj[(1:b),]
  
  freq <- seq(from, to, by = df)
  
  m = length(freq)
  
  Sj <- matrix(0, nrow = m, ncol = b)
  
  for (i in 1:m){
    
    #gaup <- (1 / (sqrt(2 * pi) * sigt)) * exp(-1.0 * (2.0 * (sin(pi * freq[i] * (dift)))^2 / (sigt^2)))
    
    #gaups <- sum(sum(gaup))
    
    for (j in 1:b){
      
      Sj[i,j] <- abs(sum(exp(1i * 2.0 * pi * freq[i] * t[which(band == j)])))^2
      
    }
    
  }
  
  Sj <- cbind(freq, Sj)
  
  topfreq <- matrix(0, nrow = m, ncol = b)
  
  for (k in 1:b){
    
    topfreq[,k] <- Sj[order(Sj[,(k+1)], decreasing = TRUE),1]
    
  }
  
  nost <- as.vector(topfreq[m:1,])
  
  if (dim(topfreq)[1] > nt){
    
    topfreq <- topfreq[1:nt,]
    
  }
  
  freq <- sort(unique(as.vector(topfreq[findInterval(topfreq, c(0.9999, 1.0001)) != 1L])))
  
  for (z in 1:length(filt)){
    
    freq <- freq[which(!(abs((1 / freq) - (1 / filt[z])) < eps*(1 / filt[z])))]
    
    freq <- freq[which(!(((1 / freq) > (1 / filt[z])) & (abs(((1 / freq) / (1 / filt[z])) - floor(((1 / freq) / (1 / filt[z])))) < eps)))]
    
    freq <- freq[which(!(((1 / freq) < (1 / filt[z])) & (abs(((1 / filt[z]) / (1 / freq)) - floor(((1 / filt[z]) / (1 / freq)))) < eps)))]
    
  }
  
  m <- length(freq)
  
  if (m < 1){
    
    freq <- 1 / (max(t) - min(t))
    
    m <- length(freq)
    
  }
  
  print(paste("Number of trial periods extracted: ", m, ".", sep=""))
  
  pwr <- rep(0, m)
  
  y <- y - mean(y)
  
  t <- t - min(t)
  
  ts <- cbind(t, y, x[,3])
  
  tsin <- matrix(c(seq(from = min(t), to = max(t), length.out = 15768), rep(mean(y), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, ts)
  
  tsin <- tsin[order(tsin[,1]),]
  
  for (k in 1:m) {
    
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
  
  pwr <- max(pwr) - pwr
  
  # Fit a exponential function to cancel out the red noise
  
  rfit <- 10^seq(log10(from), log10(to), length.out = 51)
  
  redf <- rep(0, 50)
  
  for (j in 1:(length(rfit)-1)){
    
    while (redf[j] == 0 | any(freq == redf[j])){
      
      redf[j] <- runif(1, rfit[j], rfit[j+1])
      
    }
    
  }
  
  redn <- rep(0, length(redf))
  
  for (k in 1:length(redf)) {
    
    model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*redf[k]*tsin[,1]) + cos(2*pi*redf[k]*tsin[,1]) +
                  sin(4*pi*redf[k]*tsin[,1]) + cos(4*pi*redf[k]*tsin[,1]) + sin(6*pi*redf[k]*tsin[,1]) +
                  cos(6*pi*redf[k]*tsin[,1]) + sin(8*pi*redf[k]*tsin[,1]) + cos(8*pi*redf[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    coeff <- model$coefficients
    
    coeff[is.na(coeff)] <- 0
    
    yp <- y - (coeff[2]*t + coeff[3]*sin(2*pi*redf[k]*t) + coeff[4]*cos(2*pi*redf[k]*t) +
                 coeff[5]*sin(4*pi*redf[k]*t) + coeff[6]*cos(4*pi*redf[k]*t) + coeff[7]*sin(6*pi*redf[k]*t) +
                 coeff[8]*cos(6*pi*redf[k]*t) + coeff[9]*sin(8*pi*redf[k]*t) + coeff[10]*cos(8*pi*redf[k]*t))
    
    redn[k] <- var(yp) / var(y)
    
  }
  
  redn <- max(redn) - redn
  
  arang <- seq(from = 0, to = 10, by = 0.01)
  
  mseold <- 10000
  
  for (k in 1:length(arang)){
    
    rprof <- lm(redn ~ exp(-1.0*arang[k]*redf), data = as.data.frame(cbind(redf, redn)))
    
    mse <- sum((rprof$residuals)^2)
    
    if (mse < mseold){
      
      mseold <- mse
      
      besta <- arang[k]
      
      rnco <- rprof$coefficients
      
      rnco[is.na(rnco)] <- 0
      
    }
    
  }
  
  print("Red Noise fitted using exponential profile, k*exp(-At) + C.")
  
  print(paste("k = ", rnco[1], ", A = ", besta, ", C = ", rnco[2], ".", sep=""))
  
  pwr <- pwr - (rnco[1] + rnco[2] * exp(-1.0*besta*freq))
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  if (plot == TRUE){
    
    plot(freq, pwr, pch=19, type="n", main="GLS", xlab = "Frequency (Cycles/day)", ylab = "Spectral Power", ylim = c(min(pwr), max(pwr)))
    
    lines(freq, pwr, type="l")
    
  }
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}






VRP2 <- function(x, ofac = 1, from = 0.001, to = 20, nt = 150, b = 3, plot = TRUE) {
  
  type = "frequency"
  alpha = 0.01
  level = 0
  
  eps = 1e-3
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  yb <- y
  
  n <- length(y)
  tspan <- t[n] - t[1]
  
  df <- 1 / (tspan * ofac)
  
  print(paste("Starting trial period estimation. # of trials per bin: ", nt, " and # of bins: ", b, ".", sep=""))
  
  d <- rep(0, n-1)
  
  for (k in 1:(n-1)){
    
    d[k] <- (y[k+1] - y[k]) / (t[k+1] - t[k])
    
  }
  
  sortedmags <- sort(y)
  
  E1 <- sortedmags[ceiling(0.1 * n)]
  
  E2 <- sortedmags[ceiling(0.2 * n)]
  
  E3 <- sortedmags[ceiling(0.3 * n)]
  
  E4 <- sortedmags[ceiling(0.4 * n)]
  
  E5 <- sortedmags[ceiling(0.5 * n)]
  
  E6 <- sortedmags[ceiling(0.6 * n)]
  
  E7 <- sortedmags[ceiling(0.7 * n)]
  
  E8 <- sortedmags[ceiling(0.8 * n)]
  
  E9 <- sortedmags[ceiling(0.9 * n)]
  
  Dj <- rep(0, 10)
  
  band <- rep(0, (n-1))
  
  for (i in 1:(n-1)){
    
    if ((y[i] >= min(y)) & (y[i] < E1)){
      
      Dj[1] <- Dj[1] + abs(d[k])
      
      band[i] <- 1
      
    }
    else if ((y[i] >= E1) & (y[i] < E2)){
      
      Dj[2] <- Dj[2] + abs(d[k])
      
      band[i] <- 2
      
    }
    else if ((y[i] >= E2) & (y[i] < E3)){
      
      Dj[3] <- Dj[3] + abs(d[k])
      
      band[i] <- 3
      
    }
    else if ((y[i] >= E3) & (y[i] < E4)){
      
      Dj[4] <- Dj[4] + abs(d[k])
      
      band[i] <- 4
      
    }
    else if ((y[i] >= E4) & (y[i] < E5)){
      
      Dj[5] <- Dj[5] + abs(d[k])
      
      band[i] <- 5
      
    }
    else if ((y[i] >= E5) & (y[i] < E6)){
      
      Dj[6] <- Dj[6] + abs(d[k])
      
      band[i] <- 6
      
    }
    else if ((y[i] >= E6) & (y[i] < E7)){
      
      Dj[7] <- Dj[7] + abs(d[k])
      
      band[i] <- 7
      
    }
    else if ((y[i] >= E7) & (y[i] < E8)){
      
      Dj[8] <- Dj[8] + abs(d[k])
      
      band[i] <- 8
      
    }
    else if ((y[i] >= E8) & (y[i] < E9)){
      
      Dj[9] <- Dj[9] + abs(d[k])
      
      band[i] <- 9
      
    }
    else if ((y[i] >= E9) & (y[i] <= max(y))){
      
      Dj[10] <- Dj[10] + abs(d[k])
      
      band[i] <- 10
      
    }
    
  }
  
  Dj <- as.matrix(cbind(Dj, c(1:10)))
  
  Dj <- Dj[order(Dj[,1], decreasing = T),]
  
  Dj <- Dj[(1:b),]
  
  freq <- seq(from, to, by = df)
  
  m = length(freq)
  
  Sj <- matrix(0, nrow = m, ncol = b)
  
  for (i in 1:m){
    
    #gaup <- (1 / (sqrt(2 * pi) * sigt)) * exp(-1.0 * (2.0 * (sin(pi * freq[i] * (dift)))^2 / (sigt^2)))
    
    #gaups <- sum(sum(gaup))
    
    for (j in 1:b){
      
      Sj[i,j] <- abs(sum(exp(1i * 2.0 * pi * freq[i] * t[which(band == j)])))^2
      
    }
    
  }
  
  Sj <- cbind(freq, Sj)
  
  topfreq <- matrix(0, nrow = m, ncol = b)
  
  for (k in 1:b){
    
    topfreq[,k] <- Sj[order(Sj[,(k+1)], decreasing = TRUE),1]
    
  }
  
  nost <- as.vector(topfreq[m:1,])
  
  if (dim(topfreq)[1] > nt){
    
    topfreq <- topfreq[1:nt,]
    
  }
  
  freq <- sort(unique(as.vector(topfreq)))
  
  #freq <- sort(unique(as.vector(topfreq[findInterval(topfreq, c(0.9999, 1.0001)) != 1L])))
  
  #for (z in 1:length(filt)){
    
  #  freq <- freq[which(!(abs((1 / freq) - (1 / filt[z])) < eps*(1 / filt[z])))]
    
  #  freq <- freq[which(!(((1 / freq) > (1 / filt[z])) & (abs(((1 / freq) / (1 / filt[z])) - floor(((1 / freq) / (1 / filt[z])))) < eps)))]
    
  #  freq <- freq[which(!(((1 / freq) < (1 / filt[z])) & (abs(((1 / filt[z]) / (1 / freq)) - floor(((1 / filt[z]) / (1 / freq)))) < eps)))]
    
  #}
  
  m <- length(freq)
  
  if (m < 1){
    
    freq <- 1 / (max(t) - min(t))
    
    m <- length(freq)
    
  }
  
  print(paste("Number of trial periods extracted: ", m, ".", sep=""))
  
  pwr <- rep(0, m)
  
  y <- y - mean(y)
  
  t <- t - min(t)
  
  ts <- cbind(t, y, x[,3])
  
  tsin <- matrix(c(seq(from = min(t), to = max(t), length.out = 15768), rep(mean(y), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, ts)
  
  tsin <- tsin[order(tsin[,1]),]
  
  for (k in 1:m) {
    
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
  
  pwr <- max(pwr) - pwr
  
  # Fit a exponential function to cancel out the red noise
  
  rfit <- 10^seq(log10(from), log10(to), length.out = 51)
  
  redf <- rep(0, 50)
  
  for (j in 1:(length(rfit)-1)){
    
    while (redf[j] == 0 | any(freq == redf[j])){
      
      redf[j] <- runif(1, rfit[j], rfit[j+1])
      
    }
    
  }
  
  redn <- rep(0, length(redf))
  
  for (k in 1:length(redf)) {
    
    model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*redf[k]*tsin[,1]) + cos(2*pi*redf[k]*tsin[,1]) +
                  sin(4*pi*redf[k]*tsin[,1]) + cos(4*pi*redf[k]*tsin[,1]) + sin(6*pi*redf[k]*tsin[,1]) +
                  cos(6*pi*redf[k]*tsin[,1]) + sin(8*pi*redf[k]*tsin[,1]) + cos(8*pi*redf[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    coeff <- model$coefficients
    
    coeff[is.na(coeff)] <- 0
    
    yp <- y - (coeff[2]*t + coeff[3]*sin(2*pi*redf[k]*t) + coeff[4]*cos(2*pi*redf[k]*t) +
                 coeff[5]*sin(4*pi*redf[k]*t) + coeff[6]*cos(4*pi*redf[k]*t) + coeff[7]*sin(6*pi*redf[k]*t) +
                 coeff[8]*cos(6*pi*redf[k]*t) + coeff[9]*sin(8*pi*redf[k]*t) + coeff[10]*cos(8*pi*redf[k]*t))
    
    redn[k] <- var(yp) / var(y)
    
  }
  
  redn <- max(redn) - redn
  
  arang <- seq(from = 0, to = 10, by = 0.01)
  
  mseold <- 10000
  
  for (k in 1:length(arang)){
    
    rprof <- lm(redn ~ exp(-1.0*arang[k]*redf), data = as.data.frame(cbind(redf, redn)))
    
    mse <- sum((rprof$residuals)^2)
    
    if (mse < mseold){
      
      mseold <- mse
      
      besta <- arang[k]
      
      rnco <- rprof$coefficients
      
      rnco[is.na(rnco)] <- 0
      
    }
    
  }
  
  print("Red Noise fitted using exponential profile, k*exp(-At) + C.")
  
  print(paste("k = ", rnco[1], ", A = ", besta, ", C = ", rnco[2], ".", sep=""))
  
  pwr <- pwr - (rnco[1] + rnco[2] * exp(-1.0*besta*freq))
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  if (plot == TRUE){
    
    plot(freq, pwr, pch=19, type="n", main="GLS", xlab = "Frequency (Cycles/day)", ylab = "Spectral Power", ylim = c(min(pwr), max(pwr)))
    
    lines(freq, pwr, type="l")
    
  }
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}