polygen <- function(x, bins = 100, lambda = 0.001, eps = 1e-3, delta = 1e-3, seedno = 20, pop = 1000, pairups = 200, nogen = 100, crossover = 0.65, mutation = 0.03, fdif = 0.6, dfrac = 0.7) {
  
  set.seed(seedno)
  
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
  
  # First we generate the initial population.
  
  medy <- median(y)
  
  knots1 <- runif(pop, -0.5, 0.5)
  
  knots2 <- runif(pop, -0.5, 0.5)
  
  knots3 <- runif(pop, -0.5, 0.5)
  
  knots4 <- runif(pop, -0.5, 0.5)
  
  knotsall <- cbind(knots1, knots2, knots3, knots4)
  
  m <- pop
  
  # TEMP PLEASE PUT IN ACTUAL CONVERSIONS
  
  crossover <- rep(crossover, nogen)
  
  #mutation <- rep(mutation, nogen)
  
  for (z in 1:pop){
    
    knotsall[z,] <- sort(knotsall[z,])
    
  }
  
  cost <- rep(0, pop)
  
  # Then we test them for goodness-of-fit.
  
  for (z in 1:m){
    
    knots <- sort(knotsall[z,])
    
    set1 <- matrix(x[which(x[,1] >= knots[1]),], ncol = 3)
    
    set1 <- matrix(set1[which(set1[,1] < knots[2]),], ncol = 3)
    
    set2 <- matrix(x[which(x[,1] >= knots[2]),], ncol = 3)
    
    set2 <- matrix(set2[which(set2[,1] < knots[3]),], ncol = 3)
    
    set3 <- matrix(x[which(x[,1] >= knots[3]),], ncol = 3)
    
    set3 <- matrix(set3[which(set3[,1] < knots[4]),], ncol = 3)
    
    set4 <- matrix(x[which(x[,1] >= knots[4]),], ncol = 3)
    
    setloop <- matrix(x[which(x[,1] < knots[1]),], ncol = 3)
    
    setloop[,1] <- setloop[,1] + 1
    
    set4 <- rbind(set4, setloop)
    
    if (length(set1[,1]) >= 3 & length(set2[,1]) >= 2 & length(set3[,1]) >= 2 & length(set4[,1]) >= 2){
      
      Z1 <- as.matrix(cbind(1, (set1[,1]-knots[1]), (set1[,1]-knots[1])^2))
      
      W1 <- as.matrix(diag(set1[,3]^-1))
      
      M1 <- t(Z1)%*%W1%*%Z1 + diag(rep(3,1))*lambda
      
      coeff1 <- solve(M1)%*%t(Z1)%*%W1%*%(set1[,2])
      
      a2 <- coeff1[3]*(knots[2] - knots[1])^2 + coeff1[2]*(knots[2] - knots[1]) + coeff1[1]
      
      Z2 <- as.matrix(cbind((set2[,1]-knots[2]), (set2[,1]-knots[2])^2))
      
      W2 <- as.matrix(diag(set2[,3]^-1))
      
      M2 <- t(Z2)%*%W2%*%Z2 + diag(rep(2,1))*lambda
      
      coeff2 <- solve(M2)%*%t(Z2)%*%W2%*%(set2[,2]-a2)
      
      coeff2 <- c(a2, coeff2)
      
      a3 <- coeff2[3]*(knots[3] - knots[2])^2 + coeff2[2]*(knots[3] - knots[2]) + coeff2[1]
      
      Z3 <- as.matrix(cbind((set3[,1]-knots[3]), (set3[,1]-knots[3])^2))
      
      W3 <- as.matrix(diag(set3[,3]^-1))
      
      M3 <- t(Z3)%*%W3%*%Z3 + diag(rep(2,1))*lambda
      
      coeff3 <- solve(M3)%*%t(Z3)%*%W3%*%(set3[,2]-a3)
      
      coeff3 <- c(a3, coeff3)
      
      a4 <- coeff3[3]*(knots[4] - knots[3])^2 + coeff3[2]*(knots[4] - knots[3]) + coeff3[1]
      
      Z4 <- as.matrix(rbind(cbind((set4[,1]-knots[4]), (set4[,1]-knots[4])^2), c((knots[1]+1-knots[4]), (knots[1]+1-knots[4])^2)))
      
      W4 <- as.matrix(diag(c((set4[,3]^-1), 1e+12)))
      
      M4 <- t(Z4)%*%W4%*%Z4 + diag(rep(2,1))*lambda
      
      coeff4 <- solve(M4)%*%t(Z4)%*%W4%*%(c((set4[,2]-a4), coeff1[1]-a4))
      
      coeff4 <- c(a4, coeff4)
      
      chi21 <- sum((coeff1[1] + coeff1[2]*(set1[,1] - knots[1]) + coeff1[3]*(set1[,1] - knots[1])^2 - set1[,2])^2)
      
      chi22 <- sum((coeff2[1] + coeff2[2]*(set2[,1] - knots[2]) + coeff2[3]*(set2[,1] - knots[2])^2 - set2[,2])^2)
      
      chi23 <- sum((coeff3[1] + coeff3[2]*(set3[,1] - knots[3]) + coeff3[3]*(set3[,1] - knots[3])^2 - set3[,2])^2)
      
      chi24 <- sum((coeff4[1] + coeff4[2]*(set4[,1] - knots[4]) + coeff4[3]*(set4[,1] - knots[4])^2 - set4[,2])^2)
      
      chi2 <- chi21 + chi22 + chi23 + chi24
      
      knotrep <- eps*((knots[2] - knots[1])^-2 + (knots[3] - knots[2])^-2 + (knots[4] - knots[3])^-2 + (knots[1] + 1 - knots[4])^-2)
      
      medat <- delta*((coeff1[1] - medy)^2 + (coeff2[1] - medy)^2 + (coeff3[1] - medy)^2 + (coeff4[1] - medy)^2)
      
      cost[z] <- chi2 + knotrep + medat
      
    }
    else{
      
      cost[z] <- 1e+12
      
    }
    
  }
  
  # We now initiate the genetic algorithm. Scale the knots between 0 and 1.
  
  knotsall <- knotsall + 0.5
  
  # Now encode the knots into chromosome strings.
  
  chroms1 <- as.numeric(format(round(knotsall[,1], 5), nsmall = 5)) * 10^5
  
  chroms1 <- gsub(" ", "0", sprintf("%5s", chroms1))
  
  chroms2 <- as.numeric(format(round(knotsall[,2], 5), nsmall = 5)) * 10^5
  
  chroms2 <- gsub(" ", "0", sprintf("%5s", chroms2))
  
  chroms3 <- as.numeric(format(round(knotsall[,3], 5), nsmall = 5)) * 10^5
  
  chroms3 <- gsub(" ", "0", sprintf("%5s", chroms3))
  
  chroms4 <- as.numeric(format(round(knotsall[,4], 5), nsmall = 5)) * 10^5
  
  chroms4 <- gsub(" ", "0", sprintf("%5s", chroms4))
  
  chroms <- paste(chroms1, chroms2, chroms3, chroms4, sep = "")
  
  evo <- matrix(0, nrow = nogen, ncol = 3)
  
  evo[,1] <- 1:nogen
  
  evo[1,2] <- min(as.numeric(cost))
  
  evo[1,3] <- median(as.numeric(cost))
  
  posit <- matrix(0, nrow = nogen, ncol = pop)
  
  posit[1,] <- sort(knotsall[,1])
  
  rank <- 1:m
  
  breedprob <- (1/m) + fdif*((m + 1 - 2*rank) / (m * (m + 1)))
  
  probsum <- sum(breedprob)
  
  for (k in 2:nogen){
    
    if (k %% 10 == 0){
      
      #print(paste("Generation : ", k, ".", sep=""))
      
    }
    
    # Next we breed the next generation.
    
    chroms <- chroms[order(as.numeric(cost), decreasing = FALSE)]
    
    cost <- cost[order(as.numeric(cost), decreasing = FALSE)]
    
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
        
        choosegene <- (floor(runif(1, min = 0, max = 20)) %% 20) + 1
        
        #print(paste("Pairup ", i, " procced a crossover at gene ", choosegene, "!", sep=""))
        
        chrom1 <- chroms[coups[i,1]]
        
        chrom2 <- chroms[coups[i,2]]
        
        if (choosegene == 1){
          
          newchroms[i,1] <- chrom2
          
          newchroms[i,2] <- chrom1
          
        }
        else if (choosegene == 20){
          
          newchroms[i,1] <- paste(substr(chrom1, 1, 19), substr(chrom2, 20, 20), sep="")
          
          newchroms[i,2] <- paste(substr(chrom2, 1, 19), substr(chrom1, 20, 20), sep="")
          
        }
        else{
          
          newchroms[i,1] <- paste(substr(chrom1, 1, choosegene), substr(chrom2, choosegene+1, 20), sep="")
          
          newchroms[i,2] <- paste(substr(chrom2, 1, choosegene), substr(chrom1, choosegene+1, 20), sep="")
          
        }
        
      }
      else{
        
        newchroms[i,1] <- chroms[coups[i,1]]
        
        newchroms[i,2] <- chroms[coups[i,2]]
        
      }
      
      # First we test for mutations in the first child.
      
      didwemut <- (runif(20) <= mutation[k])
      
      muthere <- which(didwemut == TRUE, arr.ind = TRUE)
      
      for (j in muthere){
        
        #print(paste("Child 1 from Pairup ", i, ", gene ", j, " procced a mutation!", sep=""))
        
        chooseval <- (floor(runif(1, min = 0, max = 10)) %% 10)
        
        if (j == 1){
          
          newchroms[i,1] <- paste(chooseval, substr(newchroms[i,1], 2, 20), sep="")
          
        }
        else if (j == 20){
          
          newchroms[i,1] <- paste(substr(newchroms[i,1], 1, 19), chooseval, sep="")
          
        }
        else{
          
          newchroms[i,1] <- paste(substr(newchroms[i,1], 1, j-1), chooseval, substr(newchroms[i,1], j+1, 20), sep="")
          
        }
        
      }
      
      # Then we test for mutations in the second child.
      
      didwemut <- (runif(20) <= mutation[k])
      
      muthere <- which(didwemut == TRUE, arr.ind = TRUE)
      
      for (j in muthere){
        
        #print(paste("Child 2 from Pairup ", i, ", gene ", j, " procced a mutation!", sep=""))
        
        chooseval <- (floor(runif(1, min = 0, max = 10)) %% 10)
        
        if (j == 1){
          
          newchroms[i,2] <- paste(chooseval, substr(newchroms[i,2], 2, 20), sep="")
          
        }
        else if (j == 20){
          
          newchroms[i,2] <- paste(substr(newchroms[i,2], 1, 19), chooseval, sep="")
          
        }
        else{
          
          newchroms[i,2] <- paste(substr(newchroms[i,2], 1, j-1), chooseval, substr(newchroms[i,2], j+1, 20), sep="")
          
        }
        
      }
      
    }
    
    #print(oldchroms)
    
    #print(newchroms)
    
    newchroms <- unique(as.vector(newchroms))
    
    for (i in 1:length(chroms)){
      
      newchroms <- newchroms[which(newchroms != chroms[i])]
      
    }
    
    fitness <- cbind(chroms, cost)
    
    testchroms1 <- (as.numeric(substr(newchroms, 1, 5)) / 10^5) - 0.5
    
    testchroms2 <- (as.numeric(substr(newchroms, 6, 10)) / 10^5) - 0.5
    
    testchroms3 <- (as.numeric(substr(newchroms, 11, 15)) / 10^5) - 0.5
    
    testchroms4 <- (as.numeric(substr(newchroms, 16, 20)) / 10^5) - 0.5
    
    testchroms <- cbind(testchroms1, testchroms2, testchroms3, testchroms4)
    
    testpop <- length(testchroms[,1])
    
    for (z in 1:testpop){
      
      testchroms[z,] <- sort(testchroms[z,])
      
    }
    
    #print(testchroms)
    
    costn <- rep(0, testpop)
    
    if (testpop > 0){
      
      for (i in 1:testpop){
        
        #if (fit == "sin"){
        
        #SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*testchroms[i]*tsin[,1]) + cos(2*pi*testchroms[i]*tsin[,1]) + 
        #              sin(4*pi*testchroms[i]*tsin[,1]) + cos(4*pi*testchroms[i]*tsin[,1]) + sin(6*pi*testchroms[i]*tsin[,1]) + cos(6*pi*testchroms[i]*tsin[,1]) + 
        #              sin(8*pi*testchroms[i]*tsin[,1]) + cos(8*pi*testchroms[i]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
        #co <- SSTlm$coefficients
        
        #co[is.na(co)] <- 0
        
        #chi2n[i] <- sum((((co[1] + co[2]*x[,1] + co[3]*sin(2*pi*testchroms[i]*x[,1]) + co[4]*cos(2*pi*testchroms[i]*x[,1]) + 
        #                     co[5]*sin(4*pi*testchroms[i]*x[,1]) + co[6]*cos(4*pi*testchroms[i]*x[,1]) + co[7]*sin(6*pi*testchroms[i]*x[,1]) + co[8]*cos(6*pi*testchroms[i]*x[,1]) +
        #                     co[9]*sin(8*pi*testchroms[i]*x[,1]) + co[10]*cos(8*pi*testchroms[i]*x[,1])) - x[,2]) / x[,3]) ^ 2.0)
        
        #}
        #else if (fit == "harm"){
        
        #SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*testchroms[i]*tsin[,1]) + cos(2*pi*testchroms[i]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
        #co <- SSTlm$coefficients
        
        #co[is.na(co)] <- 0
        
        #chi2n[i] <- sum((((co[1] + co[2]*x[,1] + co[3]*sin(2*pi*testchroms[i]*x[,1]) + co[4]*cos(2*pi*testchroms[i]*x[,1])) - x[,2]) / x[,3]) ^ 2.0)
        
        #}
        
        knots <- sort(testchroms[i,])
        
        set1 <- matrix(x[which(x[,1] >= knots[1]),], ncol = 3)
        
        set1 <- matrix(set1[which(set1[,1] < knots[2]),], ncol = 3)
        
        set2 <- matrix(x[which(x[,1] >= knots[2]),], ncol = 3)
        
        set2 <- matrix(set2[which(set2[,1] < knots[3]),], ncol = 3)
        
        set3 <- matrix(x[which(x[,1] >= knots[3]),], ncol = 3)
        
        set3 <- matrix(set3[which(set3[,1] < knots[4]),], ncol = 3)
        
        set4 <- matrix(x[which(x[,1] >= knots[4]),], ncol = 3)
        
        setloop <- matrix(x[which(x[,1] < knots[1]),], ncol = 3)
        
        setloop[,1] <- setloop[,1] + 1
        
        set4 <- rbind(set4, setloop)
        
        if (length(set1[,1]) >= 3 & length(set2[,1]) >= 2 & length(set3[,1]) >= 2 & length(set4[,1]) >= 2){
          
          Z1 <- as.matrix(cbind(1, (set1[,1]-knots[1]), (set1[,1]-knots[1])^2))
          
          W1 <- as.matrix(diag(set1[,3]^-1))
          
          M1 <- t(Z1)%*%W1%*%Z1 + diag(rep(3,1))*lambda
          
          coeff1 <- solve(M1)%*%t(Z1)%*%W1%*%(set1[,2])
          
          a2 <- coeff1[3]*(knots[2] - knots[1])^2 + coeff1[2]*(knots[2] - knots[1]) + coeff1[1]
          
          Z2 <- as.matrix(cbind((set2[,1]-knots[2]), (set2[,1]-knots[2])^2))
          
          W2 <- as.matrix(diag(set2[,3]^-1))
          
          M2 <- t(Z2)%*%W2%*%Z2 + diag(rep(2,1))*lambda
          
          coeff2 <- solve(M2)%*%t(Z2)%*%W2%*%(set2[,2]-a2)
          
          coeff2 <- c(a2, coeff2)
          
          a3 <- coeff2[3]*(knots[3] - knots[2])^2 + coeff2[2]*(knots[3] - knots[2]) + coeff2[1]
          
          Z3 <- as.matrix(cbind((set3[,1]-knots[3]), (set3[,1]-knots[3])^2))
          
          W3 <- as.matrix(diag(set3[,3]^-1))
          
          M3 <- t(Z3)%*%W3%*%Z3 + diag(rep(2,1))*lambda
          
          coeff3 <- solve(M3)%*%t(Z3)%*%W3%*%(set3[,2]-a3)
          
          coeff3 <- c(a3, coeff3)
          
          a4 <- coeff3[3]*(knots[4] - knots[3])^2 + coeff3[2]*(knots[4] - knots[3]) + coeff3[1]
          
          Z4 <- as.matrix(rbind(cbind((set4[,1]-knots[4]), (set4[,1]-knots[4])^2), c((knots[1]+1-knots[4]), (knots[1]+1-knots[4])^2)))
          
          W4 <- as.matrix(diag(c((set4[,3]^-1), 1e+12)))
          
          M4 <- t(Z4)%*%W4%*%Z4 + diag(rep(2,1))*lambda
          
          coeff4 <- solve(M4)%*%t(Z4)%*%W4%*%(c((set4[,2]-a4), coeff1[1]-a4))
          
          coeff4 <- c(a4, coeff4)
          
          chi21 <- sum((coeff1[1] + coeff1[2]*(set1[,1] - knots[1]) + coeff1[3]*(set1[,1] - knots[1])^2 - set1[,2])^2)
          
          chi22 <- sum((coeff2[1] + coeff2[2]*(set2[,1] - knots[2]) + coeff2[3]*(set2[,1] - knots[2])^2 - set2[,2])^2)
          
          chi23 <- sum((coeff3[1] + coeff3[2]*(set3[,1] - knots[3]) + coeff3[3]*(set3[,1] - knots[3])^2 - set3[,2])^2)
          
          chi24 <- sum((coeff4[1] + coeff4[2]*(set4[,1] - knots[4]) + coeff4[3]*(set4[,1] - knots[4])^2 - set4[,2])^2)
          
          chi2 <- chi21 + chi22 + chi23 + chi24
          
          knotrep <- eps*((knots[2] - knots[1])^-2 + (knots[3] - knots[2])^-2 + (knots[4] - knots[3])^-2 + (knots[1] + 1 - knots[4])^-2)
          
          medat <- delta*((coeff1[1] - medy)^2 + (coeff2[1] - medy)^2 + (coeff3[1] - medy)^2 + (coeff4[1] - medy)^2)
          
          costn[i] <- chi2 + knotrep + medat
          
        }
        else{
          
          costn[i] <- 1e+12
          
        }
        
      }
      
    }
    
    #print(c(chroms, newchroms))
    
    #knotn <- (as.numeric(c(chroms, newchroms)) / 10^5) - 0.5
    
    knotn <- c(chroms, newchroms)
    
    newfit <- cbind(newchroms, costn)
    
    fitness <- rbind(fitness, newfit)
    
    fitness <- cbind(knotn, fitness)
    
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
    
    cost <- fitness[,3]
    
    posit[k,] <- sort(fitness[,1])
    
  }
  
  #print(fitness[1:10,])
  
  plot(evo[,1], evo[,3], pch=19, type="n", main="Evolution of Genetic Algorithm", xlab = "Generation", ylab = "Chi-Squared", ylim = c(min(evo[,2]), max(evo[,3])), xlim = c(1, nogen))
  
  lines(evo[,1], evo[,2], type="l", col = "red")
  
  lines(evo[,1], evo[,3], type="l", col = "blue")
  
  #matplot(posit, pch=4)
  
  #print(paste("Completed... Period found: ", fitness[1,1], ".", sep=""))
  
  finfreq <- as.numeric(as.numeric(fitness[1,1]))
  
  knots1 <- (as.numeric(substr(fitness[1,1], 1, 5)) / 10^5) - 0.5
  
  knots2 <- (as.numeric(substr(fitness[1,1], 6, 10)) / 10^5) - 0.5
  
  knots3 <- (as.numeric(substr(fitness[1,1], 11, 15)) / 10^5) - 0.5
  
  knots4 <- (as.numeric(substr(fitness[1,1], 16, 20)) / 10^5) - 0.5
  
  knots <- c(knots1, knots2, knots3, knots4)
  
  knots <- sort(knots)
  
  sp.out <- list(knots = knots, cost = fitness[1,3])
  
}