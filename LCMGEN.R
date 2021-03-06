LCMGEN <- function(x, fit = "lsp", melc, amplc, from = 0.001, to = 0.1, jit = 0, noimod = 1, seedno = 20, redo = 5, pop = 1000, pairups = 200, nogen = 100, crossover = 0.65, mutation = 0.03, fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20) {
  
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
  
  if (fit == "ml"){
    
    compmodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]) + cos(2*pi*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
  }
  
  norm <- 1 / (2 * var(x[,2]))
  
  tot <- length(x[,1])
  
  if (fit == "bgls"){
    
    err2 <- noimod*(x[,3]*x[,3]) + rep(jit*jit, length(x[,3]))
    
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
      
    }
    else if (fit == "ce"){
      
      p <- matrix(0, nrow = bins, ncol = bins)
      
      xts <- x
      
      xts[,1] <- (xts[,1] * freq[k]) %% 1
      
      xts[,2] <- (xts[,2] - min(xts[,2]))
      
      xts[,2] <- xts[,2] / (max(xts[,2]) - min(xts[,2]))
      
      #memg <- as.numeric(tapply(xts[,2], cut(xts[,1], bins), mean))
      
      #meph <- seq(0, 1, length.out = bins+1)
      
      #meph <- meph[1:bins] + diff(meph)/2
      
      #xts <- cbind(meph, memg)
      
      #xts <- xts[complete.cases((xts)),]
      
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
      
      chkmodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*freq[k]*tsin[,1]) + cos(2*pi*freq[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      vuongres <- as.numeric(vuongtest(chkmodel, compmodel)$LRTstat)
      
      fstat[k] <- -1.0*vuongres
      
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
  
  tplot <- rep(0, nogen-1)
  
  clset <- matrix(0, nrow = 0, ncol = 3)
  
  clset2 <- matrix(0, nrow = 0, ncol = 3)
  
  clset3 <- matrix(0, nrow = 0, ncol = 3)
  
  clset4 <- matrix(0, nrow = 0, ncol = 3)
  
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
      
      filtchroms <- as.numeric(newchroms)
      
      newchroms <- newchroms[which(!is.na(filtchroms))]
      
    }
    
    fitness <- cbind(chroms, fstat)
    
    testchroms <- ((as.numeric(newchroms) / 10^10) * (to - from)) + from
    
    testchroms <- testchroms[which(testchroms >= from & testchroms <= to)]
    
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
          
        }
        else if (fit == "ce"){
          
          p <- matrix(0, nrow = bins, ncol = bins)
          
          xts <- x
          
          xts[,1] <- (xts[,1] * testchroms[i]) %% 1
          
          xts[,2] <- (xts[,2] - min(xts[,2]))
          
          xts[,2] <- xts[,2] / (max(xts[,2]) - min(xts[,2]))
          
          #memg <- tapply(xts[,2], cut(xts[,1], bins), mean)
          
          #meph <- seq(0, 1, length.out = bins+1)
          
          #meph <- meph[1:bins] + diff(meph)/2
          
          #xts <- cbind(meph, memg)
          
          #xts <- xts[complete.cases(xts),]
          
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
          
          chkmodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*testchroms[i]*tsin[,1]) + cos(2*pi*testchroms[i]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
          
          vuongres <- as.numeric(vuongtest(chkmodel, compmodel)$LRTstat)
          
          fstatn[i] <- -1.0*vuongres
          
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
    
    fitness[which(as.numeric(fitness[,1]) <= 1/tspan),1] <- 1/tspan
    
    #print(cbind((1/as.numeric(fitness[,1])), scale(1/as.numeric(fitness[,1]))))
    
    kset <- try(kmeans(scale(1/as.numeric(fitness[,1])), centers = redo, nstart = 20), TRUE)
    
    if (class(kset) == "try-error"){
    
      ## We seem to have converged to a single period. Dont bother clustering, it will not work.
      ## Output redo clusters and break the generation loop.
      
      finclus <- 1/as.numeric(fitness[,1])
      
      kmn <- rep(finclus, redo)
      
      ksd <- rep(0, redo)
      
    }
    else{
      
      kmn <- rep(0, redo)
      
      ksd <- rep(0, redo)
      
      for (l in 1:redo){
        
        kmn[l] <- mean(1/as.numeric(fitness[which(kset$cluster == l),1]))
        
        ksd[l] <- sd(1/as.numeric(fitness[which(kset$cluster == l),1]))
        
      }
      
    }
    
    #kdr <- matrix(0, nrow = redo, ncol = redo)
    
    #for (l in 1:redo){
      
    #  for (m in 1:redo){
        
    #    if (l == m){
          
    #      kdr[l,m] <- 2
          
    #    }
    #    else{
          
    #      kdr[l,m] <- (separability(1/as.numeric(fitness[which(kset$cluster == l),1]), 1/as.numeric(fitness[which(kset$cluster == m),1]), plot = F)$JM)
          
    #    }
        
    #  }
      
    #}
    
    #print(min(kdr))
    
    #tplot[k-1] <- min(kdr, na.rm = T)
    
    st <- cbind(rep(k, redo), kmn, ksd)
    
    st2 <- st[which(st[,3] < 0.00001),]
    
    ## For LPVs sometimes we need to relax this constraint. Create a new set of these clusters.
    
    st3 <- st[which(st[,3] < 0.0001),]
    
    st4 <- st[which(st[,3] < 0.001),]
    
    st5 <- st[which(st[,3] < 0.01),]
    
    clset <- rbind(clset, st2)
    
    clset2 <- rbind(clset2, st3)
    
    clset3 <- rbind(clset3, st4)
    
    clset4 <- rbind(clset4, st5)
    
    tplot[k-1] <- sd(st[,2])
    
    if (sd(st[,2]) < 0.01){
      
      brkpt <- k
      
      break
    
    }
    else{
      
      brkpt <- nogen
      
    }
    
  }
  
  clset <- clset[order(clset[,2]),, drop = F]
  
  clset2 <- clset2[order(clset2[,2]),, drop = F]
  
  clset3 <- clset3[order(clset3[,2]),, drop = F]
  
  clset4 <- clset4[order(clset4[,2]),, drop = F]
  
  clset[,2] <- 1/as.numeric(round(1/as.numeric(as.vector(clset[,2])), digits = 2))
  
  #clset <- subset(clset, !duplicated(clset[,2]))
  
  clset2[,2] <- 1/as.numeric(round(1/as.numeric(as.vector(clset2[,2])), digits = 2))
  
  #clset2 <- subset(clset2, !duplicated(clset2[,2]))
  
  clset3[,2] <- 1/as.numeric(round(1/as.numeric(as.vector(clset3[,2])), digits = 2))
  
  clset4[,2] <- 1/as.numeric(round(1/as.numeric(as.vector(clset4[,2])), digits = 2))

  kset <- try(kmeans(scale(1/as.numeric(clset[,2])), centers = redo-1, nstart = 20), TRUE)
  
  if (class(kset) == "try-error"){
    
    ## Here we apply the looser constaint set as the previous kmeans likely failed due to this.
    
    clset <- clset2
    
    try(kset <- kmeans(scale(1/as.numeric(clset[,2])), centers = redo-1, nstart = 20), TRUE)
    
  }
  
  if (class(kset) == "try-error"){
    
    ## Here we apply the even looser constaint set as the previous kmeans likely failed due to this.
    
    clset <- clset3
    
    kset <- try(kmeans(scale(1/as.numeric(clset[,2])), centers = redo-1, nstart = 20), TRUE)
    
  }
  
  if (class(kset) == "try-error"){
    
    ## Here we apply the even looser again constaint set as the previous kmeans likely failed due to this.
    
    clset <- clset4
    
    kset <- try(kmeans(scale(1/as.numeric(clset[,2])), centers = redo-1, nstart = 20), TRUE)
    
  }
  
  if (class(kset) == "try-error"){
    
    ## At this point it probably will not converge, now we just output the best clusters we have.
    
    clset <- subset(clset, !duplicated(clset[,2]),drop=F)
    
    if (is.null(nrow(clset))){
    
      clset <- clset[nrow(clset):1,]
    
      reptime <- ceiling((redo-1)/nrow(clset))
    
      clset <- clset[rep(seq_len(nrow(clset)), reptime),,drop = F][1:(redo-1),,drop = F]
    
      finfreq <- as.numeric(1/as.numeric(fitness[1,1]))
    
      kmn <- c(finfreq, clset[,2])
    
    }
    else{
      
      finfreq <- as.numeric(1/as.numeric(fitness[1,1]))
      
      kmn <- rep(finfreq, redo)
      
    }
    
    perran <- cbind(1/as.numeric(kmn), 0.9*(1/as.numeric(kmn)), 1.1*(1/as.numeric(kmn)))
    
  }
  else{
  
    kmn <- rep(0, redo)
  
    ksd <- rep(0, redo)
  
    for (l in 2:redo){
    
      kmn[l] <- mean(as.numeric(clset[which(kset$cluster == l-1),2]))
    
      ksd[l] <- sd(as.numeric(clset[which(kset$cluster == l-1),2]))
    
    }
  
    finfreq <- as.numeric(1/as.numeric(fitness[1,1]))
  
    kmn[1] <- finfreq
  
    ksd[1] <- 0
  
    perran <- cbind(1/as.numeric(kmn), 0.9*(1/as.numeric(kmn)), 1.1*(1/as.numeric(kmn)))
  
  }
  
  #print(paste(min(evo[,2]), max(evo[,2]), min(evo[,3]), max(evo[,3])))
  
  if (-1.0*max(evo[,3]) > 1e+12){
    
    close(channel)
    
    stop("Invalid Convergence")
    
  }
  
  evo <- evo[1:brkpt,]
  
  plot(2:nogen, tplot, type = "l")
  
  dots <- which((1:length(evo[,1])) %% 3 == 0)
  
  plot(evo[,1], evo[,2], type="l", col = "red", main="Evolution of Genetic Algorithm", xlab = "Generation", ylab = "Fitness Statistic", ylim = c(min(evo[,2]), max(evo[,3])), xlim = c(1, brkpt))
  
  lines(evo[dots,1], evo[dots,2], type="p", col = "red", pch = 19)
  
  lines(evo[,1], evo[,3], type="l", col = "blue")
  
  lines(evo[dots,1], evo[dots,3], type="p", col = "blue", pch = 17)
  
  legend("bottomleft", inset = .05, cex = 1, title = "Fitness", c("Minimum", "Median"), lty = c(1,1), lwd = c(2,2), col = c("red", "blue"), pch = c(19, 17), bg = "grey96")
  
  matplot(posit, pch=4, main = "Period Optimisation", xlab = "Generation", ylab = "Period (days)", xlim = c(1, brkpt))
  
  perran[which(perran[,1] == 0),] <- perran[which(perran[,1] != 0),][1,]
  
  for (i in 1:redo){
  
    lines(-14:110, rep(perran[i,1], 125), lwd = 3)
  
  }
  
  sp.out <- list(periods = perran[,1], permin = perran[,2], permax = perran[,3])
  
}