IFA <- function(x, filt, ofac = 1, from = 0.001, to = 20, nt = 150, b = 3, pixsize = 32) {
  
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
  
  plot(freq, Sj[,1], pch=19, type="n", main="Band Method power spectrum", xlab = "Frequency (Cycles/days)", ylab = "Spectral Power", ylim = c(min(as.vector(Sj)), max(as.vector(Sj))))
  
  lines(freq, Sj[,1], type="l", col="red")
  
  lines(freq, Sj[,2], type="l", col="blue")
  
  lines(freq, Sj[,3], type="l", col="green")
  
  Sj <- cbind(freq, Sj)
  
  topfreq <- matrix(0, nrow = m, ncol = b)
  
  for (k in 1:b){
    
    topfreq[,k] <- Sj[order(Sj[,(k+1)], decreasing = TRUE),1]
    
  }
  
  write.table(topfreq, file = "C:/Users/Ross/Desktop/topfreq.csv", sep = ",")
  
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
  
  for (k in 1:m) {
    
    imag <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    x <- xb
    
    x1 <- xb
    
    p <- (t / (1 / freq[k])) %% 1
    
    x[, 1] <- p
    
    x <- x[order(p),]
    
    fold <- x[,1] - 1
    
    folded <- c(fold, x[,1])
    
    foldts <- c(x[,2], x[,2])
    
    pos <- weighted.mean(x[,2], (1 / x[,3])) + rang * sd(x[,2])
    
    neg <- weighted.mean(x[,2], (1 / x[,3])) - rang * sd(x[,2])
    
    magerr <- c(x[,3], x[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    fts <- fts[which(fts[,2] <= pos & fts[,2] >= neg),]
    
    splitx <- seq(from = -1, to = 1, length.out = pixsize+1)
    
    splity <- seq(from = min(fts[,2]), to = max(fts[,2]), length.out = pixsize+1)
    
    for (j in 1:(length(splitx)-1)){
      
      imgcol <- matrix(fts[which(fts[,1] > splitx[j] & fts[,1] <= splitx[j+1]),], ncol = 3)
      
      if (as.numeric(dim(imgcol)[1]) != 0){
        
        for (l in 1:(length(splity)-1)){
          
          imgrow <- matrix(imgcol[which(imgcol[,2] > splity[l] & imgcol[,2] <= splity[l+1]),], ncol = 3)
          
          imag[j,l] <- dim(imgrow)[1]
          
        }
        
      }
      else{
        
        imag[j,] <- 0
        
      }
      
    }
    
    imgcol <- matrix(fts[which(fts[,1] > splitx[1] & fts[,1] <= splitx[2]),], ncol = 3)
    
    imgrow <- matrix(imgcol[which(imgcol[,2] > splity[1] & imgcol[,2] <= splity[2]),], ncol = 3)
    
    imag[1,1] <- as.numeric(dim(imgrow)[1])
    
    is.na(imag) <- 0
    
    imag <- imag / max(imag)
    
    imag <- imag[, rev(seq_len(ncol(imag)))]
    
    png(paste("F:/Documents/PhD/lcdump/", as.character(1/freq[k]), ".png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
  }
  
  pwr <- pwr * ((n - 1) / (2 * n))
  
  pwr <- 1 - pwr
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}