GENR <- function(x, filt, melc, amplc, ofac = 1, from = 0.001, to = 20, nt = 150, b = 3) {
  
  maxgen <- 2
  fdif <- 0.6
  pairups <- 200
  crossover <- 0.65
  mutation <- 0.003
  
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
  
  pwr <- rep(0, 3)
  
  x[,2] <- x[,2] - melc
  
  x[,2] <- x[,2] / amplc
  
  x[,3] <- x[,3] / amplc
  
  x <- x[which(x[,2] < 1),]
  
  x <- x[which(x[,2] > -1),]
  
  xb <- x
  
  # Now we build the genetic algorithm.
  
  # Initially we must test the Chi-Squared value of the initial frequencies.
  
  chi2 <- rep(0, m)
  
  tsin <- matrix(c(seq(from = min(x[,1]), to = max(x[,1]), length.out = 15768), rep(weighted.mean(x[,2], (1 / x[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, x)
  
  tsin <- tsin[order(tsin[,1]),]
  
  for (k in 1:m){
    
    SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*freq[k]*tsin[,1]) + cos(2*pi*freq[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    co <- SSTlm$coefficients
    
    co[is.na(co)] <- 0
    
    chi2[k] <- sum((((co[1] + co[2]*x[,1] + co[3]*sin(2*pi*freq[k]*x[,1]) + co[4]*cos(2*pi*freq[k]*x[,1])) - x[,2]) / x[,3]) ^ 2.0)
    
  }
  
  # First we need to scale the periods between 0 and 1.
  
  freq <- (freq - from) / (to - from)
  
  # Now encode the initial periods in string form.
  
  chroms <- as.numeric(format(round(freq, 10), nsmall = 10)) * 10^10
  
  chroms <- gsub(" ", "0", sprintf("%10s", chroms))
  
  for (k in 2:maxgen){
    
    print(paste("Generation : ", k, ".", sep=""))
    
    # Next we breed the next generation.
    
    chroms <- chroms[order(chi2, decreasing = FALSE)]
    
    rank <- m:1
    
    breedprob <- (1/m) + fdif*((m + 1 - 2*rank) / (m * (m + 1)))
    
    probsum <- sum(breedprob)
    
    coups <- matrix(0, nrow = pairups, ncol = 2)
    
    for (i in 1:pairups){
      
      val <- runif(1)
      
      for (j in 1:m){
        
        val <- val - breedprob[j]
        
        if (val <= 0){
          
          coups[i,1] <- j
          
          break
          
        }
        
        coups[i,1] <- 1
        
      }
      
      val <- runif(1)
      
      for (j in 1:m){
        
        val <- val - breedprob[j]
        
        if (val <= 0){
          
          coups[i,2] <- j
          
          break
          
        }
        
        coups[i,2] <- 1
        
      }
      
    }
    
    print(chroms)
    
    print(coups)
    
    # Now we perform crossover and mutation if the probability allows for each couple...
    
    newchroms <- matrix("", nrow = pairups, ncol = 2)
    
    oldchroms <- newchroms
    
    for (i in 1:pairups){
      
      oldchroms[i,1] <- chroms[coups[i,1]]
      
      oldchroms[i,2] <- chroms[coups[i,2]]
      
      if (runif(1) <= crossover){
        
        print(paste("Pairup ", i, " procced a crossover!", sep=""))
        
        choosegene <- (floor(runif(1, min = 0, max = 10)) %% 10) + 1
        
        chrom1 <- chroms[coups[i,1]]
        
        chrom2 <- chroms[coups[i,2]]
        
        if (choosegene == 1){
          
          newchroms[i,1] <- chrom2
          
          newchroms[i,2] <- chrom1
          
        }
        else if (choosegene == 10){
          
          newchroms[i,1] <- paste(substr(chrom1, 1, 9), substr(chrom2, 10, 10), sep="")
          
          newchroms[i,2] <- paste(substr(chrom2, 1, 9), substr(chrom1, 10, 10), sep="")
          
        }
        else{
          
          newchroms[i,1] <- paste(substr(chrom1, 1, choosegene), substr(chrom2, choosegene+1, 10), sep="")
          
          newchroms[i,2] <- paste(substr(chrom2, 1, choosegene), substr(chrom1, choosegene+1, 10), sep="")
          
        }
        
      }
      else{
        
        newchroms[i,1] <- chroms[coups[i,1]]
        
        newchroms[i,2] <- chroms[coups[i,2]]
        
      }
      
      # First we test for mutations in the first child.
      
      didwemut <- (runif(10) <= mutation)
      
      muthere <- which(didwemut == TRUE, arr.ind = TRUE)
      
      for (j in muthere){
        
        print(paste("Child 1 from Pairup ", i, " position ", j, " procced a mutation!", sep=""))
        
        chooseval <- (floor(runif(1, min = 0, max = 10)) %% 10)
        
        if (j == 1){
          
          newchroms[i,1] <- paste(chooseval, substr(newchroms[i,1], 2, 10), sep="")
          
        }
        else if (j == 10){
          
          newchroms[i,1] <- paste(substr(newchroms[i,1], 1, 9), chooseval, sep="")
          
        }
        else{
          
          newchroms[i,1] <- paste(substr(newchroms[i,1], 1, j-1), chooseval, substr(newchroms[i,1], j+1, 10), sep="")
          
        }
        
      }
      
      # Then we test for mutations in the second child.
      
      didwemut <- (runif(10) <= mutation)
      
      muthere <- which(didwemut == TRUE, arr.ind = TRUE)
      
      for (j in muthere){
        
        print(paste("Child 2 from Pairup ", i, " position ", j, " procced a mutation!", sep=""))
        
        chooseval <- (floor(runif(1, min = 0, max = 10)) %% 10)
        
        if (j == 1){
          
          newchroms[i,2] <- paste(chooseval, substr(newchroms[i,2], 2, 10), sep="")
          
        }
        else if (j == 10){
          
          newchroms[i,2] <- paste(substr(newchroms[i,2], 1, 9), chooseval, sep="")
          
        }
        else{
          
          newchroms[i,2] <- paste(substr(newchroms[i,2], 1, j-1), chooseval, substr(newchroms[i,2], j+1, 10), sep="")
          
        }
        
      }
      
    }
    
    print(oldchroms)
    
    print(newchroms)
    
  }
  
  pwr <- pwr * ((n - 1) / (2 * n))
  
  pwr <- 1 - pwr
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}