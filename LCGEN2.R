LCGEN2 <- function(x, spurs, fit = "sin", melc, amplc, from = 0.001, to = 20, seedno = 20, pop = 1000, pairups = 200, nogen = 100, crossover = 0.65, mutation = 0.03, fdif = 0.6, dfrac = 0.7) {
  
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
  
  #freq <- runif(pop, min = from, max = to)
  
  m = length(freq)
  
  pwr <- rep(0, 3)
  
  x[,2] <- x[,2] - melc
  
  x[,2] <- x[,2] / amplc
  
  x[,3] <- x[,3] / amplc
  
  x <- x[which(x[,2] < 1), , drop = F]
  
  x <- x[which(x[,2] > -1), , drop = F]
  
  xb <- x
  
  # Now we build the genetic algorithm.
  
  # Initially we must test the Chi-Squared value of the initial frequencies.
  
  #print(paste("Generation : 1."))
  
  chi2 <- rep(0, m)
  
  tsin <- matrix(c(seq(from = min(x[,1]), to = max(x[,1]), length.out = 15768), rep(weighted.mean(x[,2], (1 / x[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, x)
  
  tsin <- tsin[order(tsin[,1]),]
  
  for (k in 1:m){
    
    if (fit == "harm"){
      
      SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*freq[k]*tsin[,1]) + cos(2*pi*freq[k]*tsin[,1]) + 
                    sin(4*pi*freq[k]*tsin[,1]) + cos(4*pi*freq[k]*tsin[,1]) + sin(6*pi*freq[k]*tsin[,1]) + cos(6*pi*freq[k]*tsin[,1]) + 
                    sin(8*pi*freq[k]*tsin[,1]) + cos(8*pi*freq[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      co <- SSTlm$coefficients
      
      co[is.na(co)] <- 0
      
      chi2[k] <- sum(((x[,2] - (co[1] + co[2]*x[,1] + co[3]*sin(2*pi*freq[k]*x[,1]) + co[4]*cos(2*pi*freq[k]*x[,1]) + 
                                  co[5]*sin(4*pi*freq[k]*x[,1]) + co[6]*cos(4*pi*freq[k]*x[,1]) + co[7]*sin(6*pi*freq[k]*x[,1]) + co[8]*cos(6*pi*freq[k]*x[,1]) + 
                                  co[9]*sin(8*pi*freq[k]*x[,1]) + co[10]*cos(8*pi*freq[k]*x[,1]))) / x[,3]) ^ 2.0)
      
    }
    
    else if (fit == "sin"){
      
      SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*freq[k]*tsin[,1]) + cos(2*pi*freq[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      co <- SSTlm$coefficients
      
      co[is.na(co)] <- 0
      
      chi2[k] <- sum(((x[,2] - (co[1] + co[2]*x[,1] + co[3]*sin(2*pi*freq[k]*x[,1]) + co[4]*cos(2*pi*freq[k]*x[,1]))) / x[,3]) ^ 2.0)
      
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
  
  evo[1,2] <- min(as.numeric(chi2))
  
  evo[1,3] <- median(as.numeric(chi2))
  
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
    
    chroms <- chroms[order(as.numeric(chi2), decreasing = FALSE)]
    
    chi2 <- chi2[order(as.numeric(chi2), decreasing = FALSE)]
    
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
    
    fitness <- cbind(chroms, chi2)
    
    testchroms <- ((as.numeric(newchroms) / 10^10) * (to - from)) + from
    
    testpop <- length(testchroms)
    
    chi2n <- rep(0, testpop)
    
    if (testpop > 0){
      
      for (i in 1:testpop){
        
        if (fit == "sin"){
          
          SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*testchroms[i]*tsin[,1]) + cos(2*pi*testchroms[i]*tsin[,1]) + 
                        sin(4*pi*testchroms[i]*tsin[,1]) + cos(4*pi*testchroms[i]*tsin[,1]) + sin(6*pi*testchroms[i]*tsin[,1]) + cos(6*pi*testchroms[i]*tsin[,1]) + 
                        sin(8*pi*testchroms[i]*tsin[,1]) + cos(8*pi*testchroms[i]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
          
          co <- SSTlm$coefficients
          
          co[is.na(co)] <- 0
          
          chi2n[i] <- sum((((co[1] + co[2]*x[,1] + co[3]*sin(2*pi*testchroms[i]*x[,1]) + co[4]*cos(2*pi*testchroms[i]*x[,1]) + 
                               co[5]*sin(4*pi*testchroms[i]*x[,1]) + co[6]*cos(4*pi*testchroms[i]*x[,1]) + co[7]*sin(6*pi*testchroms[i]*x[,1]) + co[8]*cos(6*pi*testchroms[i]*x[,1]) +
                               co[9]*sin(8*pi*testchroms[i]*x[,1]) + co[10]*cos(8*pi*testchroms[i]*x[,1])) - x[,2]) / x[,3]) ^ 2.0)
          
        }
        else if (fit == "harm"){
          
          SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*testchroms[i]*tsin[,1]) + cos(2*pi*testchroms[i]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
          
          co <- SSTlm$coefficients
          
          co[is.na(co)] <- 0
          
          chi2n[i] <- sum((((co[1] + co[2]*x[,1] + co[3]*sin(2*pi*testchroms[i]*x[,1]) + co[4]*cos(2*pi*testchroms[i]*x[,1])) - x[,2]) / x[,3]) ^ 2.0)
          
        }
        
      }
      
    }
    
    period <- 1.0 / (((as.numeric(c(chroms, newchroms)) / 10^10) * (to - from)) + from)
    
    newfit <- cbind(newchroms, chi2n)
    
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
    
    chi2 <- fitness[,3]
    
    posit[k,] <- sort(fitness[,1])
    
  }
  
  #print(fitness[1:10,])
  
  plot(evo[,1], evo[,3], pch=19, type="n", main="Evolution of Genetic Algorithm", xlab = "Generation", ylab = "Chi-Squared", ylim = c(min(evo[,2]), max(evo[,3])), xlim = c(1, nogen))
  
  lines(evo[,1], evo[,2], type="l", col = "red")
  
  lines(evo[,1], evo[,3], type="l", col = "blue")
  
  matplot(posit, pch=4)
  
  print(paste("Completed... Period found: ", fitness[1,1], ".", sep=""))
  
  finfreq <- as.numeric(1/as.numeric(fitness[1,1]))
  
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
  
  sp.out <- list(period = as.character(fitness[1,1]), chi2 = as.character(fitness[1,3]), model = SSTlm, onemodel = onemodel, othmodel = othmodel, othper = othper)
  
}