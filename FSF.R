FSF <- function(x, filt, melc, amplc, ofac = 1, from = 0.001, to = 20, nt = 150, b = 3) {
  
  type = "frequency"
  alpha = 0.01
  level = 0
  
  eps = 1e-3
  rang = 4
  
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
  
  plot(0, type="n", main = "Bands method colour-coded into ten bands", xlab="Modified Julian Date (days)", ylab="Magnitude", xlim=c(min(t), max(t)), ylim=c(max(y), min(y)))
  
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
      
      points(t[i], y[i], type = "p", pch = 19, col = "#FF0000")
      
    }
    else if ((y[i] >= E1) & (y[i] < E2)){
      
      Dj[2] <- Dj[2] + abs(d[k])
      
      band[i] <- 2
      
      points(t[i], y[i], type = "p", pch = 19, col = "#E2571E")
      
    }
    else if ((y[i] >= E2) & (y[i] < E3)){
      
      Dj[3] <- Dj[3] + abs(d[k])
      
      band[i] <- 3
      
      points(t[i], y[i], type = "p", pch = 19, col = "#FF7F00")
      
    }
    else if ((y[i] >= E3) & (y[i] < E4)){
      
      Dj[4] <- Dj[4] + abs(d[k])
      
      band[i] <- 4
      
      points(t[i], y[i], type = "p", pch = 19, col = "#FFFF00")
      
    }
    else if ((y[i] >= E4) & (y[i] < E5)){
      
      Dj[5] <- Dj[5] + abs(d[k])
      
      band[i] <- 5
      
      points(t[i], y[i], type = "p", pch = 19, col = "#00FF00")
      
    }
    else if ((y[i] >= E5) & (y[i] < E6)){
      
      Dj[6] <- Dj[6] + abs(d[k])
      
      band[i] <- 6
      
      points(t[i], y[i], type = "p", pch = 19, col = "#00FFFF")
      
    }
    else if ((y[i] >= E6) & (y[i] < E7)){
      
      Dj[7] <- Dj[7] + abs(d[k])
      
      band[i] <- 7
      
      points(t[i], y[i], type = "p", pch = 19, col = "#0000FF")
      
    }
    else if ((y[i] >= E7) & (y[i] < E8)){
      
      Dj[8] <- Dj[8] + abs(d[k])
      
      band[i] <- 8
      
      points(t[i], y[i], type = "p", pch = 19, col = "#4B0082")
      
    }
    else if ((y[i] >= E8) & (y[i] < E9)){
      
      Dj[9] <- Dj[9] + abs(d[k])
      
      band[i] <- 9
      
      points(t[i], y[i], type = "p", pch = 19, col = "#8B00FF")
      
    }
    else if ((y[i] >= E9) & (y[i] <= max(y))){
      
      Dj[10] <- Dj[10] + abs(d[k])
      
      band[i] <- 10
      
      points(t[i], y[i], type = "p", pch = 19, col = "#FF1493")
      
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
  
  plot(freq, Sj[,1], pch=19, type="n", main="Band Method power spectrum", xlab = "Frequency (Cycles/days)", ylab = "Spectral Power", ylim = c(min(as.vector(Sj)), max(as.vector(Sj))))
  
  lines(freq, Sj[,1], type="l", col="red")
  
  lines(freq, Sj[,2], type="l", col="blue")
  
  lines(freq, Sj[,3], type="l", col="green")
  
  Sj <- cbind(freq, Sj)
  
  topfreq <- matrix(0, nrow = m, ncol = b)
  
  for (k in 1:b){
    
    topfreq[,k] <- Sj[order(Sj[,(k+1)], decreasing = TRUE),1]
    
  }
  
  write.table(topfreq, file = "M:/topfreq.csv", sep = ",")
  
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
  
  x[,2] <- x[,2] - melc
  
  x[,2] <- x[,2] / amplc
  
  x[,3] <- x[,3] / amplc 
  
  x <- x[which(x[,2] < 1),]
  
  x <- x[which(x[,2] > -1),]
  
  xb <- x
  
  ansmat <- matrix(0, nrow = m, ncol = 24)
  
  for (k in 1:m) {
    
    x <- xb
    
    x1 <- xb
    
    p <- (((x[,1] - x[which.min(x[,2]),1]) / (1 / freq[k])) + 0.25) %% 1
    
    x[, 1] <- p
    
    x <- x[order(p),]
    
    fold <- x[,1] - 1
    
    folded <- c(fold, x[,1])
    
    foldts <- c(x[,2], x[,2])
    
    magerr <- c(x[,3], x[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    ftsin <- matrix(c(seq(from = -1, to = 1, length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    ftsin <- rbind(ftsin, fts)
    
    ftsin <- ftsin[order(ftsin[,1]),]
    
    #model <- lm(ftsin[,2] ~ 0 + sin(2*pi*ftsin[,1]) + cos(2*pi*ftsin[,1]) +
    #              sin(4*pi*ftsin[,1]) + cos(4*pi*ftsin[,1]) + sin(6*pi*ftsin[,1]) +
    #              cos(6*pi*ftsin[,1]) + sin(8*pi*ftsin[,1]) + cos(8*pi*ftsin[,1]) +
    #              sin(10*pi*ftsin[,1]) + cos(10*pi*ftsin[,1]) + sin(12*pi*ftsin[,1]) + cos(12*pi*ftsin[,1]) +
    #              sin(14*pi*ftsin[,1]) + cos(14*pi*ftsin[,1]) + sin(16*pi*ftsin[,1]) + cos(16*pi*ftsin[,1]), data = as.data.frame(ftsin[,1:2]), weights = 1/ftsin[,3])
    
    #coeff <- model$coefficients
    
    #coeff[is.na(coeff)] <- 0
    
    Z <- cbind(sin(2*pi*ftsin[,1]), cos(2*pi*ftsin[,1]), sin(4*pi*ftsin[,1]), cos(4*pi*ftsin[,1]), sin(6*pi*ftsin[,1]), cos(6*pi*ftsin[,1]), sin(8*pi*ftsin[,1]), cos(8*pi*ftsin[,1]),
               sin(10*pi*ftsin[,1]), cos(10*pi*ftsin[,1]), sin(12*pi*ftsin[,1]), cos(12*pi*ftsin[,1]), sin(14*pi*ftsin[,1]), cos(14*pi*ftsin[,1]), sin(16*pi*ftsin[,1]), cos(16*pi*ftsin[,1]))
    
    Z <- as.matrix(Z)
    
    lambda <- 10^(-3)
    
    M <- t(Z)%*%Z + diag(rep(16,1))*lambda
    
    coeff <- solve(M)%*%t(Z)%*%(ftsin[,2])
    
    #if (abs((1 / freq[k]) - 1.00) < 0.01){
    
    #plot(fts[,1], fts[,2], main = paste("", 1/freq[k], " Folded Light Curve spread from -1.0 to 1.0", sep=""), xlab = paste("Phase at ", (1 / freq[k]), " days", sep=""),
    #   ylab = "Scaled Magnitude", pch=19, ylim = c(1, -1))
    
    #lines(sort(ftsin[,1]), model$fitted.values, col = 2)
    
    #(paste("Period : ", 1/freq[k], " days", sep = ""))
    
    #print(coeff)
    
    #}
    
    A1 <- (coeff[1]^2.0 + coeff[2]^2.0)^0.5
    
    A2 <- (coeff[3]^2.0 + coeff[4]^2.0)^0.5
    
    A3 <- (coeff[5]^2.0 + coeff[6]^2.0)^0.5
    
    A4 <- (coeff[7]^2.0 + coeff[8]^2.0)^0.5
    
    A5 <- (coeff[9]^2.0 + coeff[10]^2.0)^0.5
    
    A6 <- (coeff[11]^2.0 + coeff[12]^2.0)^0.5
    
    A7 <- (coeff[13]^2.0 + coeff[14]^2.0)^0.5
    
    A8 <- (coeff[15]^2.0 + coeff[16]^2.0)^0.5
    
    NPH1 <- atan2(coeff[2], coeff[1])
    
    NPH2 <- atan2(coeff[4], coeff[3]) - NPH1
    
    NPH3 <- atan2(coeff[6], coeff[5]) - NPH1
    
    NPH4 <- atan2(coeff[8], coeff[7]) - NPH1
    
    NPH5 <- atan2(coeff[10], coeff[9]) - NPH1
    
    NPH6 <- atan2(coeff[12], coeff[11]) - NPH1
    
    NPH7 <- atan2(coeff[14], coeff[13]) - NPH1
    
    NPH8 <- atan2(coeff[16], coeff[15]) - NPH1
    
    PH1 <- 0
    
    PH2 <- atan2(sin(NPH2), cos(NPH2))
    
    PH3 <- atan2(sin(NPH3), cos(NPH3))
    
    PH4 <- atan2(sin(NPH4), cos(NPH4))
    
    PH5 <- atan2(sin(NPH5), cos(NPH5))
    
    PH6 <- atan2(sin(NPH6), cos(NPH6))
    
    PH7 <- atan2(sin(NPH7), cos(NPH7))
    
    PH8 <- atan2(sin(NPH8), cos(NPH8))
    
    ansmat[k,] <- c(freq[k], as.numeric(A1), as.numeric(A2), as.numeric(A3), as.numeric(A4),
                    as.numeric(A5), as.numeric(A6), as.numeric(A7), as.numeric(A8),
                    as.numeric(NPH1), as.numeric(NPH2), as.numeric(NPH3), as.numeric(NPH4),
                    as.numeric(NPH5), as.numeric(NPH6), as.numeric(NPH7), as.numeric(NPH8),
                    as.numeric(PH2), as.numeric(PH3), as.numeric(PH4), as.numeric(PH5),
                    as.numeric(PH6), as.numeric(PH7), as.numeric(PH8))
    
  }
  
  ansmat <- as.data.frame(ansmat)
  
  colnames(ansmat) <- c("Period", "A1", "A2", "A3", "A4", "A5", "A6", "A7",
                        "A8", "NPH1", "NPH2", "NPH3", "NPH4", "NPH5", "NPH6",
                        "NPH7", "NPH8", "PH2", "PH3", "PH4",
                        "PH5", "PH6", "PH7", "PH8")
  
  write.table(ansmat, file = "M:/ansmat.csv", sep = ",")
  
  pwr <- pwr * ((n - 1) / (2 * n))
  
  pwr <- 1 - pwr
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}