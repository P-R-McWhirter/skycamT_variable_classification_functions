briergen <- function(probs, trainclass, seedno = 20, pop = 1000, pairups = 200, nogen = 100, crossover = 0.65, mutation = 0.03, fdif = 0.6, dfrac = 0.7) {
  
  set.seed(seedno)
  
  A <- runif(pop, -20, 20)
  
  B <- runif(pop, -20, 20)
  
  m <- pop
  
  # TEMP PLEASE PUT IN ACTUAL CONVERSIONS
  
  crossover <- rep(crossover, nogen)
  
  #mutation <- rep(mutation, nogen)
  
  cost <- rep(0, pop)
  
  # Then we test them for goodness-of-fit.
  
  for (z in 1:m){
    
    problargest <- as.numeric(as.vector(apply(probs, 1, max)))
    
    prob2ndlargest <- as.numeric(as.vector(apply(probs, 1, function(x) sort(x)[dim(probs)[2]-1])))
    
    probmarg <- problargest - prob2ndlargest
    
    r <- 1 / (1 + exp(A[z]*probmarg + B[z]))
    
    corprobs <- probs*(1-r)
    
    maxprobpos <- as.numeric(as.vector(match(colnames(probs)[max.col(probs,ties.method="first")], names(probs))))
    
    for (i in 1:nrow(probs)){
      
      corprobs[i,maxprobpos[i]] <- probs[i,maxprobpos[i]] + r[i]*(1-probs[i,maxprobpos[i]])
      
    }
    
    cost[z] <- (1/nrow(corprobs)) * sum(sum((trainclass - corprobs)^2.0))
    
  }
  
  # We now initiate the genetic algorithm. Scale A and B between 0 and 1.
  
  A <- (A + 20) / 40
  
  B <- (B + 20) / 40
  
  # Now encode the knots into chromosome strings.
  
  chroms1 <- as.numeric(format(round(A, 10), nsmall = 10)) * 10^10
  
  chroms1 <- gsub(" ", "0", sprintf("%10s", chroms1))
  
  chroms2 <- as.numeric(format(round(B, 10), nsmall = 10)) * 10^10
  
  chroms2 <- gsub(" ", "0", sprintf("%10s", chroms2))
  
  chroms <- paste(chroms1, chroms2, sep = "")
  
  evo <- matrix(0, nrow = nogen, ncol = 3)
  
  evo[,1] <- 1:nogen
  
  evo[1,2] <- min(as.numeric(cost))
  
  evo[1,3] <- median(as.numeric(cost))
  
  posit <- matrix(0, nrow = nogen, ncol = pop)
  
  posit[1,] <- A
  
  rank <- 1:m
  
  breedprob <- (1/m) + fdif*((m + 1 - 2*rank) / (m * (m + 1)))
  
  probsum <- sum(breedprob)
  
  for (k in 2:nogen){
    
    if (k %% 10 == 0){
      
      print(paste("Generation : ", k, ".", sep=""))
      
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
    
    testchroms1 <- ((as.numeric(substr(newchroms, 1, 10)) / 10^10) * 40) - 20
    
    testchroms2 <- ((as.numeric(substr(newchroms, 11, 20)) / 10^10) * 40) - 20
    
    testchroms <- cbind(testchroms1, testchroms2)
    
    testpop <- length(testchroms[,1])
    
    if (testpop > 0){
    
      for (z in 1:testpop){
      
        testchroms[z,] <- sort(testchroms[z,])
      
      }
      
    }
    
    #print(testchroms)
    
    costn <- rep(0, testpop)
    
    if (testpop > 0){
      
      for (i in 1:testpop){
        
        testA <- testchroms[i,1]
        
        testB <- testchroms[i,2]
        
        r <- 1 / (1 + exp(testA*probmarg + testB))
        
        corprobs <- probs*(1-r)
        
        for (j in 1:nrow(probs)){
          
          corprobs[j,maxprobpos[j]] <- probs[j,maxprobpos[j]] + r[j]*(1-probs[j,maxprobpos[j]])
          
        }
        
        costn[i] <- (1/nrow(corprobs)) * sum(sum((trainclass - corprobs)^2.0))
        
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
  
  dots <- which((1:length(evo[,1])) %% 3 == 0)
  
  plot(evo[,1], evo[,2], type="l", col = "red", main="Evolution of Genetic Algorithm", xlab = "Generation", ylab = "Fitness Statistic", ylim = c(min(evo[,2]), max(evo[,3])), xlim = c(1, nogen))
  
  lines(evo[dots,1], evo[dots,2], type="p", col = "red", pch = 19)
  
  lines(evo[,1], evo[,3], type="l", col = "blue")
  
  lines(evo[dots,1], evo[dots,3], type="p", col = "blue", pch = 17)
  
  legend("topright", inset = .05, cex = 1, title = "Fitness", c("Minimum", "Median"), lty = c(1,1), lwd = c(2,2), col = c("red", "blue"), pch = c(19, 17), bg = "grey96")
  
  #matplot(posit, pch=4)
  
  #print(paste("Completed... Period found: ", fitness[1,1], ".", sep=""))
  
  finfreq <- as.numeric(as.numeric(fitness[1,1]))
  
  A <- ((as.numeric(substr(fitness[1,1], 1, 10)) / 10^10) * 40) - 20
  
  B <- ((as.numeric(substr(fitness[1,1], 11, 20)) / 10^10) * 40) - 20
  
  sp.out <- list(A = A, B = B, cost = fitness[1,3])
  
}