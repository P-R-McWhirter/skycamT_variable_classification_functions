LCQGEN <- function(x, spurs, fit = "lsp", melc, amplc, from = 0.001, to = 20, seedno = 20, pop = 1000, pairups = 200, nogen = 100, crossover = 0.65, mutation = 0.03, fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20) {
  
  # This is a modification to the sinusoid and harmonic sinusoid method through direct incorporation of astronomical periodogram methods.
  
  set.seed(seedno)
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
  
  #print(paste("Initiating Genetic Algorithm with Population: ", pop, " and # of pairups: ", pairups, ".", sep=""))
  
  #print(paste("Algorithm will run ", nogen, " total generations.", sep=""))
  
  logpop <- floor(pop/2)
  
  freq1 <- 10^runif(logpop, min = log10(from), max = log10(to))
  
  freq2 <- 1.0 / runif((pop-logpop), min = (1/to), max = (1/from))
  
  freq <- c(freq1, freq2)
  
  freq <- sort(freq)
  
  #freq <- runif(pop, min = 1/to, max = 1/from)
  
  m = length(freq)
  
  pwr <- rep(0, 3)
  
  x[,2] <- x[,2] - melc
  
  x[,2] <- x[,2] / amplc
  
  x[,3] <- x[,3] / amplc
  
  xo <- x
  
  prerem <- length(x[,1])
  
  x <- x[which(x[,2] < 1), , drop = F]
  
  x <- x[which(x[,2] > -1), , drop = F]
  
  postrem <- length(x[,1])
  
  if ((postrem / prerem) < clipcut){
    
    x <- xo
    
  }
  
  xb <- x
  
  # Now we build the genetic algorithm.
  
  # Initially we must test the Chi-Squared value of the initial frequencies.
  
  #print(paste("Generation : 1."))
  
  fstat <- rep(0, m)
  
  tsin <- matrix(c(seq(from = min(x[,1]), to = max(x[,1]), length.out = 15768), rep(weighted.mean(x[,2], (1 / x[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, x)
  
  tsin <- tsin[order(tsin[,1]),]
  
  norm <- 1 / (2 * var(x[,2]))
  
  tot <- length(x[,1])
  
  if (fit == "bgls"){
    
    err2 <- x[,3]*x[,3]
    
    err2 <- rep(10, length(x[,3]))
    
    w <- 1 / err2
    
    W <- sum(w)
    
    bigY <- sum(w*x[,2])
    
    p <- rep(0, m)
    
    constants <- rep(0, m)
    
    exponents <- rep(0, m)
    
  }
  
  for (k in 1:m){
    
    if (fit == "lsp"){
      
      wi <- 2 * pi * freq[k]
      
      tau <- 0.5 * atan2(sum(sin(wi * x[,1])), sum(cos(wi * x[,1])))/wi
      
      arg <- wi * (x[,1] - tau)
      
      cs <- cos(arg)
      
      sn <- sin(arg)
      
      A <- (sum(x[,2] * cs))^2
      
      B <- sum(cs * cs)
      
      C <- (sum(x[,2] * sn))^2
      
      D <- sum(sn * sn)
      
      PN <- norm * (A/B + C/D)
      
      fstat[k] <- -PN
      
      for (z in 1:length(spurs)){
        
        fstat[k] <- fstat[k] * (1 - 0.5*exp(-1.0 * ((freq[k] - (spurs[z]))^2.0 / 2.0*(50/freq[k])^2.0)))
        
      }
      
    }
    else if (fit == "bgls"){
      
      wi <- 2 * pi * freq[k]
      
      theta <- 0.5 * atan2(sum(w*sin(2*wi*x[,1])), sum(w*cos(2*wi*x[,1])))
      
      st <- (wi * x[,1]) - theta
      
      cosx <- cos(st)
      sinx <- sin(st)
      wcosx <- w*cosx
      wsinx <- w*sinx
      
      C <- sum(wcosx)
      S <- sum(wsinx)
      
      YCh <- sum(x[,2]*wcosx)
      YSh <- sum(x[,2]*wsinx)
      CCh <- sum(wcosx*cosx)
      SSh <- sum(wsinx*sinx)
      
      if (CCh != 0 & SSh != 0){
        
        K <- (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
        
        L <- (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
        
        M <- (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
        
        constants[k] <- 1/(sqrt(CCh*SSh*abs(K)))
        
      }
      else if (CCh == 0){
        
        K <- (S*S - W*SSh)/(2*SSh)
        
        L <- (bigY*SSh - S*YSh)/(SSh)
        
        M <- (YSh*YSh)/(2*SSh)
        
        constants[k] <- 1/(sqrt(SSh*abs(K)))
        
      }
      else if (SSh == 0){
        
        K <- (C*C - W*CCh)/(2*CCh)
        
        L <- (bigY*CCh - C*YCh)/(CCh)
        
        M <- (YCh*YCh)/(2*CCh)
        
        constants[k] <- 1/(sqrt(CCh*abs(K)))
        
      }
      
      if (K > 0){
        
        print("K is positive, This should not happen.")
        
      }
      
      exponents[k] <- (M - ((L*L)/(4*K)))
      
      logp <- log10(constants[k]) + (exponents[k] * log10(exp(1)))
      
      p <- 10^(as.numeric(logp))
      
      if (p < 0.00001){
        
        p <- 0
        
      }
      
      fstat[k] <- -p
      
      for (z in 1:length(spurs)){
        
        fstat[k] <- fstat[k] * (1 - 0.5*exp(-1.0 * ((freq[k] - (spurs[z]))^2.0 / 2.0*(50/freq[k])^2.0)))
        
      }
      
    }
    else if (fit == "ce"){
      
      p <- matrix(0, nrow = bins, ncol = bins)
      
      xts <- x
      
      xts[,1] <- (xts[,1] * freq[k]) %% 1
      
      xts[,2] <- (xts[,2] - min(xts[,2]))
      
      xts[,2] <- xts[,2] / (max(xts[,2]) - min(xts[,2]))
      
      memg <- as.numeric(tapply(xts[,2], cut(xts[,1], bins), mean))
      
      meph <- seq(0, 1, length.out = bins+1)
      
      meph <- meph[1:bins] + diff(meph)/2
      
      xts <- cbind(meph, memg)
      
      xts <- xts[complete.cases((xts)),]
      
      ab = matrix(c(0, 0, 1, 1), 2, 2)
      
      datbin <- bin2(cbind(xts[,1], xts[,2]), ab = ab, nbin = c(bins, bins))
      
      p <- datbin$nc / tot
      
      pj <- colSums(p)
      
      po <- log(pj / p)
      
      po[which(po == Inf | po == -Inf | is.nan(po))] <- 0
      
      Hc <- sum(rowSums(p * po))
      
      fstat[k] <- Hc
      
    }
    else if (fit == "ml"){
      
      fitlc <- lcsplinefit(AAVSO_expset[4355,], spur = T, pero = 1/freq[k], quiet = T)
      
      fitlc <- c(1/freq[k], fitlc[c(5, 8:106)])
      
      intlc1 <- t(as.matrix(as.numeric(spline.pca$rotation[,1]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc2 <- t(as.matrix(as.numeric(spline.pca$rotation[,2]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc3 <- t(as.matrix(as.numeric(spline.pca$rotation[,3]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc4 <- t(as.matrix(as.numeric(spline.pca$rotation[,4]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc5 <- t(as.matrix(as.numeric(spline.pca$rotation[,5]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc6 <- t(as.matrix(as.numeric(spline.pca$rotation[,6]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc7 <- t(as.matrix(as.numeric(spline.pca$rotation[,7]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc8 <- t(as.matrix(as.numeric(spline.pca$rotation[,8]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
      
      intlc <- as.data.frame(c(fitlc[1:2], intlc1, intlc2, intlc3, intlc4, intlc5, intlc6, intlc7, intlc8))
      
      colnames(intlc) <- c("Period", "Amplitude", paste("PCA", 1:8, sep = ""))
      
      pred <- predict(DCEPsplinemodel, intlc, type = "prob")
      
      fstat[k] <- -1.0*as.numeric(pred[,1])
      
    }
    
  }
  
  # First we need to scale the periods between 0 and 1.
  
  intfreq <- freq
  
  freq <- (freq - from) / (to - from)
  
  # Now encode the initial periods in string form.
  
  chroms <- as.numeric(format(round(freq, 10), nsmall = 10)) * 10^10
  
  chroms <- gsub(" ", "0", sprintf("%10s", chroms))
  
  evo <- matrix(0, nrow = nogen, ncol = 3)
  
  evo[,1] <- 1:nogen
  
  evo[1,2] <- min(as.numeric(fstat))
  
  evo[1,3] <- median(as.numeric(fstat))
  
  posit <- matrix(0, nrow = nogen, ncol = pop)
  
  posit[1,] <- sort(1 / intfreq)
  
  rank <- 1:m
  
  breedprob <- (1/m) + fdif*((m + 1 - 2*rank) / (m * (m + 1)))
  
  probsum <- sum(breedprob)
  
  for (k in 2:nogen){
    
    if (k %% 10 == 0){
      
      #print(paste("Generation : ", k, ".", sep=""))
      
    }
    
    # Next we breed the next generation.
    
    chroms <- chroms[order(as.numeric(fstat), decreasing = FALSE)]
    
    fstat <- fstat[order(as.numeric(fstat), decreasing = FALSE)]
    
    coups <- matrix(0, nrow = pairups, ncol = 2)
    
    for (i in 1:pairups){
      
      val <- runif(1)
      
      for (j in 1:m){
        
        val <- val - breedprob[j]
        
        if (val <= 0){
          
          coups[i,1] <- j
          
          break
          
        }
        
        coups[i,1] <- m
        
      }
      
      val <- runif(1)
      
      for (j in 1:m){
        
        val <- val - breedprob[j]
        
        if (val <= 0){
          
          coups[i,2] <- j
          
          break
          
        }
        
        coups[i,2] <- m
        
      }
      
    }
    
    #print(chroms)
    
    #print(coups)
    
    # Now we perform crossover and mutation if the probability allows for each couple...
    
    newchroms <- matrix("", nrow = pairups, ncol = 2)
    
    oldchroms <- newchroms
    
    for (i in 1:pairups){
      
      oldchroms[i,1] <- chroms[coups[i,1]]
      
      oldchroms[i,2] <- chroms[coups[i,2]]
      
      if (runif(1) <= crossover[k]){
        
        choosegene <- (floor(runif(1, min = 0, max = 10)) %% 10) + 1
        
        #print(paste("Pairup ", i, " procced a crossover at gene ", choosegene, "!", sep=""))
        
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
      
      didwemut <- (runif(10) <= mutation[k])
      
      muthere <- which(didwemut == TRUE, arr.ind = TRUE)
      
      for (j in muthere){
        
        #print(paste("Child 1 from Pairup ", i, ", gene ", j, " procced a mutation!", sep=""))
        
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
      
      didwemut <- (runif(10) <= mutation[k])
      
      muthere <- which(didwemut == TRUE, arr.ind = TRUE)
      
      for (j in muthere){
        
        #print(paste("Child 2 from Pairup ", i, ", gene ", j, " procced a mutation!", sep=""))
        
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
    
    #print(oldchroms)
    
    #print(newchroms)
    
    newchroms <- unique(as.vector(newchroms))
    
    for (i in 1:length(chroms)){
      
      newchroms <- newchroms[which(newchroms != chroms[i])]
      
    }
    
    fitness <- cbind(chroms, fstat)
    
    testchroms <- ((as.numeric(newchroms) / 10^10) * (to - from)) + from
    
    testpop <- length(testchroms)
    
    fstatn <- rep(0, testpop)
    
    if (testpop > 0){
      
      for (i in 1:testpop){
        
        if (fit == "bgls"){
          
          p <- rep(0, testpop)
          
          constants <- rep(0, testpop)
          
          exponents <- rep(0, testpop)
          
        }
        
        if (fit == "lsp"){
          
          wi <- 2 * pi * testchroms[i]
          
          tau <- 0.5 * atan2(sum(sin(wi * x[,1])), sum(cos(wi * x[,1])))/wi
          
          arg <- wi * (x[,1] - tau)
          
          cs <- cos(arg)
          
          sn <- sin(arg)
          
          A <- (sum(x[,2] * cs))^2
          
          B <- sum(cs * cs)
          
          C <- (sum(x[,2] * sn))^2
          
          D <- sum(sn * sn)
          
          PN <- norm * (A/B + C/D)
          
          fstatn[i] <- -PN
          
          for (z in 1:length(spurs)){
            
            fstatn[i] <- fstatn[i] * (1 - 0.5*exp(-1.0 * ((testchroms[i] - (spurs[z]))^2.0 / 2.0*(50/testchroms[i])^2.0)))
            
          }
          
        }
        else if (fit == "bgls"){
          
          wi <- 2 * pi * testchroms[i]
          
          theta <- 0.5 * atan2(sum(w*sin(2*wi*x[,1])), sum(w*cos(2*wi*x[,1])))
          
          st <- (wi * x[,1]) - theta
          
          cosx <- cos(st)
          sinx <- sin(st)
          wcosx <- w*cosx
          wsinx <- w*sinx
          
          C <- sum(wcosx)
          S <- sum(wsinx)
          
          YCh <- sum(x[,2]*wcosx)
          YSh <- sum(x[,2]*wsinx)
          CCh <- sum(wcosx*cosx)
          SSh <- sum(wsinx*sinx)
          
          if (CCh != 0 & SSh != 0){
            
            K <- (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
            
            L <- (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
            
            M <- (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
            
            constants[i] <- 1/(sqrt(CCh*SSh*abs(K)))
            
          }
          else if (CCh == 0){
            
            K <- (S*S - W*SSh)/(2*SSh)
            
            L <- (bigY*SSh - S*YSh)/(SSh)
            
            M <- (YSh*YSh)/(2*SSh)
            
            constants[i] <- 1/(sqrt(SSh*abs(K)))
            
          }
          else if (SSh == 0){
            
            K <- (C*C - W*CCh)/(2*CCh)
            
            L <- (bigY*CCh - C*YCh)/(CCh)
            
            M <- (YCh*YCh)/(2*CCh)
            
            constants[i] <- 1/(sqrt(CCh*abs(K)))
            
          }
          
          if (K > 0){
            
            print("K is positive, This should not happen.")
            
          }
          
          exponents[i] <- (M - L*L/(4*K))
          
          logp <- log10(constants[i]) + (exponents[i] * log10(exp(1)))
          
          p <- 10^(as.numeric(logp))
          
          fstatn[i] <- -p
          
          for (z in 1:length(spurs)){
            
            fstatn[i] <- fstatn[i] * (1 - 0.5*exp(-1.0 * ((testchroms[i] - (spurs[z]))^2.0 / 2.0*(50/testchroms[i])^2.0)))
            
          }
          
        }
        else if (fit == "ce"){
          
          p <- matrix(0, nrow = bins, ncol = bins)
          
          xts <- x
          
          xts[,1] <- (xts[,1] * testchroms[i]) %% 1
          
          xts[,2] <- (xts[,2] - min(xts[,2]))
          
          xts[,2] <- xts[,2] / (max(xts[,2]) - min(xts[,2]))
          
          memg <- tapply(xts[,2], cut(xts[,1], bins), mean)
          
          meph <- seq(0, 1, length.out = bins+1)
          
          meph <- meph[1:bins] + diff(meph)/2
          
          xts <- cbind(meph, memg)
          
          xts <- xts[complete.cases(xts),]
          
          ab = matrix(c(0, 0, 1, 1), 2, 2)
          
          datbin <- bin2(cbind(xts[,1], xts[,2]), ab = ab, nbin = c(bins, bins))
          
          p <- datbin$nc / tot
          
          pj <- colSums(p)
          
          po <- log(pj / p)
          
          po[which(po == Inf | po == -Inf | is.nan(po))] <- 0
          
          Hc <- sum(rowSums(p * po))
          
          fstatn[i] <- Hc
          
        }
        else if (fit == "ml"){
          
          fitlc <- lcsplinefit(asas_train[1,], spur = T, pero = 1/testchroms[i], quiet = T)
          
          fitlc <- c(1/testchroms[i], fitlc[c(5, 8:106)])
          
          intlc1 <- t(as.matrix(as.numeric(spline.pca$rotation[,1]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc2 <- t(as.matrix(as.numeric(spline.pca$rotation[,2]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc3 <- t(as.matrix(as.numeric(spline.pca$rotation[,3]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc4 <- t(as.matrix(as.numeric(spline.pca$rotation[,4]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc5 <- t(as.matrix(as.numeric(spline.pca$rotation[,5]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc6 <- t(as.matrix(as.numeric(spline.pca$rotation[,6]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc7 <- t(as.matrix(as.numeric(spline.pca$rotation[,7]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc8 <- t(as.matrix(as.numeric(spline.pca$rotation[,8]))) %*% (as.matrix(as.numeric(fitlc[3:101])))
          
          intlc <- as.data.frame(c(fitlc[1:2], intlc1, intlc2, intlc3, intlc4, intlc5, intlc6, intlc7, intlc8))
          
          colnames(intlc) <- c("Period", "Amplitude", paste("PCA", 1:8, sep = ""))
          
          pred <- predict(DCEPsplinemodel, intlc, type = "prob")
          
          fstatn[i] <- -1.0*as.numeric(pred[,1])
          
        }
        
      }
      
    }
    
    period <- 1.0 / (((as.numeric(c(chroms, newchroms)) / 10^10) * (to - from)) + from)
    
    newfit <- cbind(newchroms, fstatn)
    
    fitness <- rbind(fitness, newfit)
    
    fitness <- cbind(period, fitness)
    
    fitness <- fitness[order(as.numeric(fitness[,3]), decreasing = FALSE),]
    
    randfrac <- floor(dfrac * pop)
    
    dsamp <- sample((randfrac+1):length(fitness[,2]), pop-randfrac, replace=F)
    
    if (randfrac <= 0){
      
      fitness <- fitness[dsamp,]
      
    }
    else if (randfrac >= 1){
      
      fitness <- fitness[1:pop,]
      
    }
    else{
      
      fitness <- fitness[c(1:randfrac, dsamp),]
      
    }
    
    # Let's try to make a better population update method
    
    #fitlen <- length(fitness[,2])-1
    
    #rank <- fitlen:1
    
    #deathprob <- (1/fitlen) + ddif*((fitlen + 1 - 2*rank) / (fitlen * (fitlen + 1)))
    
    #probsum <- sum(deathprob)
    
    #deathpop <- fitlen - pop
    
    #cull <- c()
    
    #while (length(unique(cull)) <= deathpop){
    
    #val <- runif(1)
    
    #for (j in 1:fitlen){
    
    #val <- val - deathprob[j]
    
    #if (val <= 0){
    
    #cull <- c(cull, j)
    
    #cull <- unique(cull)
    
    #break
    
    #}
    
    #cull <- c(cull, fitlen)
    
    #}
    
    #}
    
    #fitness <- fitness[-(cull+1),]
    
    evo[k,2] <- min(as.numeric(fitness[,3]))
    
    evo[k,3] <- median(as.numeric(fitness[,3]))
    
    #print(paste("Chi-squared of fittest member of generation ", k, " is : ", evo[k,2], ".", sep=""))
    
    #print(paste("Median Chi-squared of generation ", k, " is : ", evo[k,3], ".", sep=""))
    
    chroms <- fitness[,2]
    
    fstat <- fitness[,3]
    
    posit[k,] <- sort(fitness[,1])
    
  }
  
  finfreq <- as.numeric(1/as.numeric(fitness[1,1]))
  
  plot(evo[,1], evo[,3], pch=19, type="n", main="Evolution of Genetic Algorithm", xlab = "Generation", ylab = "Chi-Squared", ylim = c(min(evo[,2]), max(evo[,3])), xlim = c(1, nogen))
  
  lines(evo[,1], evo[,2], type="l", col = "red")
  
  lines(evo[,1], evo[,3], type="l", col = "blue")
  
  write.table(posit, file = "F:/Documents/PhD/posit.csv", sep = ",")
  
  matplot(posit, pch=4)
  
  print(paste("Completed... Period found: ", 1/finfreq, ".", sep=""))
  
  SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*finfreq*tsin[,1]) + cos(2*pi*finfreq*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  onemodel <- lm(tsin[,2] ~ 1, data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  mul2model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*finfreq*tsin[,1]*(1/2)) + cos(2*pi*finfreq*tsin[,1]*(1/2)), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  mul3model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*finfreq*tsin[,1]*(1/3)) + cos(2*pi*finfreq*tsin[,1]*(1/3)), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  smul2model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*finfreq*tsin[,1]*2) + cos(2*pi*finfreq*tsin[,1]*2), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  smul3model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*finfreq*tsin[,1]*3) + cos(2*pi*finfreq*tsin[,1]*3), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  al1model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*abs(finfreq + 1)*tsin[,1]) + cos(2*pi*abs(finfreq + 1)*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  aln1model <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*abs(finfreq - 1)*tsin[,1]) + cos(2*pi*abs(finfreq - 1)*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  alhmodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*abs(finfreq + 2)*tsin[,1]) + cos(2*pi*abs(finfreq + 2)*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  alnhmodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*abs(finfreq - 2)*tsin[,1]) + cos(2*pi*abs(finfreq - 2)*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  othmodel <- list(mul2model, mul3model, smul2model, smul3model, al1model, aln1model, alhmodel, alnhmodel)
  
  othper <- 1 / c((1/2)*finfreq, (1/3)*finfreq, 2*finfreq, 3*finfreq, abs(finfreq+1), abs(finfreq-1), abs(finfreq+2), abs(finfreq-2))
  
  sp.out <- list(period = as.character(1/finfreq), fstat = as.character(fitness[1,3]), model = SSTlm, onemodel = onemodel, othmodel = othmodel, othper = othper)
  
}