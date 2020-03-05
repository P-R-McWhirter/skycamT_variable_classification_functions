lcgpolyfit <- function(data, lim = 100, detrend = FALSE, fixint = FALSE, spur = FALSE, pop = 100, pairups = 20, nogen = 100, crossover = 0.6, mutation = 0.03, fdif = 0.6, dfrac = 0.7, eps = 0.01, lambda = 0.001, radii = 0.1, maxp = 20000, pero = 1.0, seedno = 10) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = 18)
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    keep <- sigclip(info$Rcat, sig = 3, tol = 0.000001)
    
    info <- info[keep,]
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    if (detrend == TRUE){
      
      fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts), weights = 1/ts[,3])
      
      coeff <- coefficients(fit)
      
      ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
      
    }
    
    allp <- length(ts[,1])
    
    print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
    
    if (maxp < length(ts[,1])){
      
      print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    #amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    fts <- ts
    
    fts[,1] <- ((((fts[,1] - fts[,1][which.max(abs(fts[,2]))]) / (per[i])) + 0.5) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (fixint == TRUE){
      
      meanmag <- weighted.mean(fts[,2], (1/fts[,3]))
      
      ftsin <- matrix(c(seq(from = min(fts[,1]), to = max(fts[,1]), length.out = 100), rep(meanmag, times = 100), rep(1/0, times = 100)), nrow = 100)
      
      fts <- rbind(fts, ftsin)
      
      fts <- fts[order(fts[,1]),]
      
    }
    
    knots <- as.vector(c(-0.4, -0.1, 0.1, 0.4))
    
    knotsold <- knots
    
    cost <- 1e+12
    
    costold <- cost
    
    okay <- F
    
    for (z in 1:maxit){
      
      if (cost >= costold){
        
        knots <- knotsold
        
      }
      else{
        
        okay <- F
        
      }
      
      costold <- cost
      
      if (z > 1){
        
        knotsold <- knots
        
        while (okay == F){
          
          knots <- sort(((knots + 0.5 + rnorm(knots, delta)) %% 1) - 0.5)
          
          set1 <- matrix(fts[which(fts[,1] >= knots[1]),], ncol = 3)
          
          set1 <- matrix(set1[which(set1[,1] < knots[2]),], ncol = 3)
          
          set2 <- matrix(fts[which(fts[,1] >= knots[2]),], ncol = 3)
          
          set2 <- matrix(set2[which(set2[,1] < knots[3]),], ncol = 3)
          
          set3 <- matrix(fts[which(fts[,1] >= knots[3]),], ncol = 3)
          
          set3 <- matrix(set3[which(set3[,1] < knots[4]),], ncol = 3)
          
          set4 <- matrix(fts[which(fts[,1] >= knots[4]),], ncol = 3)
          
          setloop <- matrix(fts[which(fts[,1] < knots[1]),], ncol = 3)
          
          setloop[,1] <- setloop[,1] + 1
          
          set4 <- rbind(set4, setloop)
          
          if (length(set1[,1]) >= 3 & length(set2[,1]) >= 2 & length(set3[,1]) >= 2 & length(set4[,1]) >= 2){
            
            okay <- T
            
          }
          
        }
        
      }
      
      set1 <- matrix(fts[which(fts[,1] >= knots[1]),], ncol = 3)
      
      set1 <- matrix(set1[which(set1[,1] < knots[2]),], ncol = 3)
      
      set2 <- matrix(fts[which(fts[,1] >= knots[2]),], ncol = 3)
      
      set2 <- matrix(set2[which(set2[,1] < knots[3]),], ncol = 3)
      
      set3 <- matrix(fts[which(fts[,1] >= knots[3]),], ncol = 3)
      
      set3 <- matrix(set3[which(set3[,1] < knots[4]),], ncol = 3)
      
      set4 <- matrix(fts[which(fts[,1] >= knots[4]),], ncol = 3)
      
      setloop <- matrix(fts[which(fts[,1] < knots[1]),], ncol = 3)
      
      setloop[,1] <- setloop[,1] + 1
      
      set4 <- rbind(set4, setloop)
      
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
      
      cost <- chi2 + knotrep
      
    }
    
    fts <- fts[which(fts[,3] != 1/0),]
    
    plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
         ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
    
    set1in <- matrix(c(seq(from = knots[1], to = knots[2], length.out = 1000), rep(mean(set1[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
    
    set2in <- matrix(c(seq(from = knots[2], to = knots[3], length.out = 1000), rep(mean(set2[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
    
    set3in <- matrix(c(seq(from = knots[3], to = knots[4], length.out = 1000), rep(mean(set3[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
    
    set4in <- matrix(c(seq(from = knots[4], to = (knots[1]+1), length.out = 1000), rep(mean(set4[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
    
    set1in <- rbind(set1in, set1)
    
    set1in <- set1in[order(set1in[,1]),]
    
    set2in <- rbind(set2in, set2)
    
    set2in <- set2in[order(set2in[,1]),]
    
    set3in <- rbind(set3in, set3)
    
    set3in <- set3in[order(set3in[,1]),]
    
    set4in <- rbind(set4in, set4)
    
    set4in <- set4in[order(set4in[,1]),]
    
    set4in1 <- set4in[which(set4in[,1] <= 0.5),]
    
    set4in2 <- set4in[which(set4in[,1] > 0.5),]
    
    preset1in <- coeff1[1] + coeff1[2]*(set1in[,1] - knots[1]) + coeff1[3]*(set1in[,1] - knots[1])^2
    
    preset2in <- coeff2[1] + coeff2[2]*(set2in[,1] - knots[2]) + coeff2[3]*(set2in[,1] - knots[2])^2
    
    preset3in <- coeff3[1] + coeff3[2]*(set3in[,1] - knots[3]) + coeff3[3]*(set3in[,1] - knots[3])^2
    
    preset4in1 <- coeff4[1] + coeff4[2]*(set4in1[,1] - knots[4]) + coeff4[3]*(set4in1[,1] - knots[4])^2
    
    preset4in2 <- coeff4[1] + coeff4[2]*(set4in2[,1] - knots[4]) + coeff4[3]*(set4in2[,1] - knots[4])^2
    
    set4in2[,1] <- set4in2[,1] - 1
    
    lines(set1in[,1], preset1in, col='red', lwd=2)
    
    lines(set2in[,1], preset2in, col='red', lwd=2)
    
    lines(set3in[,1], preset3in, col='red', lwd=2)
    
    lines(set4in1[,1], preset4in1, col='red', lwd=2)
    
    lines(set4in2[,1], preset4in2, col='red', lwd=2)
    
    #feats <- predict(lo, seq(-0.5, 0.5, length.out = len)) - weighted.mean(fts[,2], (1 / fts[,3]))
    
    print(paste("Final knots: ", knots[1], " ", knots[2], " ", knots[3], " ", knots[4], ".", sep = ""))
    
    feats <- c(allp, knots, coeff1, coeff2, coeff3, coeff4, chi2)
    
    fullfeats[i,] <- feats
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, fullfeats, log10(lambda), bins)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "#.Observations", "Knot.1", "Knot.2", "Knot.3", "Knot.4", "Poly.1.co.1", "Poly.1.co.2", "Poly.1.co.3",
                           "Poly.2.co.1", "Poly.2.co.2", "Poly.2.co.3", "Poly.3.co.1", "Poly.3.co.2", "Poly.3.co.3", "Poly.4.co.1", "Poly.4.co.2", "Poly.4.co.3", "Chi.Squared", "Log.Regularisation", "No.Bins")
  
  fullfeats
  
}





lcgbinpolyfit <- function(data, lim = 100, detrend = FALSE, zeroed = TRUE, fixint = FALSE, spur = FALSE, pop = 100, pairups = 20, nogen = 100, crossover = 0.65, mutation = 0.03, fdif = 0.6, dfrac = 0.7, eps = 0.01, delta = 0.01, lambda = 0.001, radii = 0.1, maxp = 20000, pero = 1.0, bins = 100, seedno = 10) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = 18)
  
  clipped <- rep(1, length(data$Name))
  
  lcamp <- rep(0, length(data$Name))
  
  mnmag <- rep(0, length(data$Name))
  
  sdmag <- rep(0, length(data$Name))
  
  cordist <- rep(0, length(data$Name))
  
  pzpm <- rep(0, length(data$Name))
  
  pzpr <- rep(0, length(data$Name))
  
  pfzp <- rep(0, length(data$Name))
  
  actbins <- rep(0, length(data$Name))
  
  phase <- seq(-0.5, 0.5, by = 0.001)
  
  fitmag <- rep(0, length(phase))
  
  intamp <- rep(0, length(data$Name))
  
  normag <- rep(0, length(data$Name))
  
  intphase <- seq(-0.49, 0.49, by = 0.01)
  
  intmag <- matrix(0, ncol = length(intphase), nrow = length(data$Name))
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    cordist[i] <- min(geodist)*60*60
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    keep <- sigclip(info$Rcat, sig = 3, tol = 0.000001)
    
    info <- info[keep,]
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    if (detrend == TRUE){
      
      fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts), weights = 1/ts[,3])
      
      coeff <- coefficients(fit)
      
      ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
      
    }
    
    allp <- length(ts[,1])
    
    print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
    
    if (maxp < length(ts[,1])){
      
      print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    lcamp[i] <- amplitude
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    ts <- as.matrix(ts, ncol = 3)
    
    fts <- ts
    
    prerem <- length(ts[,1])
    
    fts[,1] <- ((((fts[,1] - fts[,1][which.max(abs(fts[,2] - weighted.mean(fts[,2], (1 / fts[,3]))))]) / (per[i])) + 0.0) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (bins < 2){
      
      print("Number of bins must be 2 or more, setting number of bins to 2.")
      
      bins <- 2
      
    }
    
    splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
    
    tery <- rep(0, (bins-1))
    
    tere <- tery
    
    tert <- tery
    
    for (k in 1:(bins-1)){
      
      tery[k] <- mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2, drop = F])
      
      tere[k] <- mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F])
      
      #tere[k] <- abs((max(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2]) - min(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2])) / 2)
      
      tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
      
    }
    
    #tere[which(tere == 0)] <- min(tere[which(tere != 0)])
    
    fts <- cbind(tert, tery, tere)
    
    row.has.na <- apply(fts, 1, function(x){any(is.na(x))})
    
    fts <- fts[!row.has.na, , drop = F]
    
    fts <- fts[order(fts[,1]), , drop = F]
    
    actbins[i] <- length(fts[,1])
    
    mnmag[i] <- weighted.mean(fts[,2], (1 / fts[,3]))
    
    sdmag[i] <- sd(fts[,2])
    
    if (zeroed == TRUE){
      
      fts[,2] <- fts[,2] - weighted.mean(fts[,2], (1 / fts[,3]))
      
    }
    
    # Here we reset the 0.0 phase to the minimum brightness.
    
    # We will calculate the mean phase of the maximum 5th percentile of the binned light curve.
    
    maxcut <- median(tail(sort(fts[,2]), ceiling(length(fts[,2]) * 0.05)))
    
    # Now use mean the phase bins with magnitudes below this to determine the minimum phase.
    
    cutfts <- fts[which(fts[,2] >= maxcut),]
    
    if (is.null(dim(cutfts))){
      
      print("5th Percentile cut not sufficient for this light curve, performing a 20th percentile cut...")
      
      maxcut <- median(tail(sort(fts[,2]), ceiling(length(fts[,2]) * 0.2)))
      
    }
    
    if (is.null(dim(cutfts))){
      
      print("20th Percentile cut not sufficient for this light curve, choosing minimum point as phase...")
      
      mopr <- fts[which.max(fts[,2]),1]
      
      pzpm[i] <- mopr
      
      pzpr[i] <- 1
      
    }
    else{
      
      mopr <- mean(cutfts[,1])
      
      pzpm[i] <- mopr
      
      pzpr[i] <- abs(max(cutfts[,1]) - min(cutfts[,1]))
      
    }
    
    # Now we adjust the folded time series so that this phase is 0.0.
    
    fts[,1] <- fts[,1] - mopr
    
    fts[which(fts[,1] < -0.5),1] <- fts[which(fts[,1] < -0.5),1] + 1
    
    fts[which(fts[,1] > 0.5),1] <- fts[which(fts[,1] > 0.5),1] - 1
    
    # And we are done!
    
    medy <- median(fts[,2])
    
    if (fixint == TRUE){
      
      meanmag <- weighted.mean(fts[,2], (1/fts[,3]))
      
      ftsin <- matrix(c(seq(from = min(fts[,1]), to = max(fts[,1]), length.out = 100), rep(meanmag, times = 100), rep(1/0, times = 100)), nrow = 100)
      
      fts <- rbind(fts, ftsin)
      
      fts <- fts[order(fts[,1]),]
      
    }
    
    genres <- polygen(fts, bins = bins, lambda = lambda, eps = eps, delta = delta, seedno = seedno, pop = pop, pairups = pairups, nogen = nogen, crossover = crossover, mutation = mutation, fdif = fdif, dfrac = dfrac)
    
    knots <- genres$knots
    
    bestcost <- genres$cost
    
    if (bestcost == 1e+12){
      
      print("Algorithm failed to converge.")
      
      knots <- c(NA, NA, NA, NA)
      
      coeff1 <- c(NA, NA, NA)
      
      coeff2 <- coeff1
      
      coeff3 <- coeff1
      
      coeff4 <- coeff1
      
      feats <- c(allp, knots, coeff1, coeff2, coeff3, coeff4, NA)
      
      intmag[i,] <- NA
      
      fullfeats[i,] <- feats
      
    }
    else{
      
      set1 <- matrix(fts[which(fts[,1] >= knots[1]),], ncol = 3)
      
      set1 <- matrix(set1[which(set1[,1] < knots[2]),], ncol = 3)
      
      set2 <- matrix(fts[which(fts[,1] >= knots[2]),], ncol = 3)
      
      set2 <- matrix(set2[which(set2[,1] < knots[3]),], ncol = 3)
      
      set3 <- matrix(fts[which(fts[,1] >= knots[3]),], ncol = 3)
      
      set3 <- matrix(set3[which(set3[,1] < knots[4]),], ncol = 3)
      
      set4 <- matrix(fts[which(fts[,1] >= knots[4]),], ncol = 3)
      
      setloop <- matrix(fts[which(fts[,1] < knots[1]),], ncol = 3)
      
      setloop[,1] <- setloop[,1] + 1
      
      set4 <- rbind(set4, setloop)
      
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
      
      cost <- chi2 + knotrep + medat
      
      
      # Extract coefficients and knots from parameters
      
      # Use polyfit to perform one final zero phase alignment.
      
      fitmag[which(phase < knots[1])] <- coeff4[1] +
        coeff4[2]*(phase[which(phase < knots[1])] - knots[4] + 1) +
        coeff4[3]*(phase[which(phase < knots[1])] - knots[4] + 1)^2
      
      fitmag[which(phase >= knots[1] & phase < knots[2])] <- coeff1[1] +
        coeff1[2]*(phase[which(phase >= knots[1] & phase < knots[2])] - knots[1]) +
        coeff1[3]*(phase[which(phase >= knots[1] & phase < knots[2])] - knots[1])^2
      
      fitmag[which(phase >= knots[2] & phase < knots[3])] <- coeff2[1] +
        coeff2[2]*(phase[which(phase >= knots[2] & phase < knots[3])] - knots[2]) +
        coeff2[3]*(phase[which(phase >= knots[2] & phase < knots[3])] - knots[2])^2
      
      fitmag[which(phase >= knots[3] & phase < knots[4])] <- coeff3[1] +
        coeff3[2]*(phase[which(phase >= knots[3] & phase < knots[4])] - knots[3]) +
        coeff3[3]*(phase[which(phase >= knots[3] & phase < knots[4])] - knots[3])^2
      
      fitmag[which(phase >= knots[4])] <- coeff4[1] +
        coeff4[2]*(phase[which(phase >= knots[4])] - knots[4]) +
        coeff4[3]*(phase[which(phase >= knots[4])] - knots[4])^2
      
      zerophase <- phase[which(fitmag == max(fitmag), arr.ind = T)]
      
      pfzp[i] <- zerophase
      
      knots <- knots - zerophase
      
      fts[,1] <- fts[,1] - zerophase
      
      fts[which(fts[,1] < -0.5),1] <- fts[which(fts[,1] < -0.5),1] + 1
      
      fts[which(fts[,1] > 0.5),1] <- fts[which(fts[,1] > 0.5),1] - 1
      
      fts <- fts[order(fts[,1]),]
      
      knots <- ((knots + 0.5) %% 1) - 0.5
      
      knotord <- order(knots)
      
      coeff1b <- coeff1
      
      coeff2b <- coeff2
      
      coeff3b <- coeff3
      
      coeff4b <- coeff4
      
      if (knotord[1] == 2){
        
        coeff1 <- coeff2b
        
      }
      else if (knotord[1] == 3){
        
        coeff1 <- coeff3b
        
      }
      else if (knotord[1] == 4){
        
        coeff1 <- coeff4b
        
      }
      
      if (knotord[2] == 1){
        
        coeff2 <- coeff1b
        
      }
      else if (knotord[2] == 3){
        
        coeff2 <- coeff3b
        
      }
      else if (knotord[2] == 4){
        
        coeff2 <- coeff4b
        
      }
      
      if (knotord[3] == 1){
        
        coeff3 <- coeff1b
        
      }
      else if (knotord[3] == 2){
        
        coeff3 <- coeff2b
        
      }
      else if (knotord[3] == 4){
        
        coeff3 <- coeff4b
        
      }
      
      if (knotord[4] == 1){
        
        coeff4 <- coeff1b
        
      }
      else if (knotord[4] == 2){
        
        coeff4 <- coeff2b
        
      }
      else if (knotord[4] == 3){
        
        coeff4 <- coeff3b
        
      }
      
      knots <- knots[order(knots)]
      
      
      set1 <- matrix(fts[which(fts[,1] >= knots[1]),], ncol = 3)
      
      set1 <- matrix(set1[which(set1[,1] < knots[2]),], ncol = 3)
      
      set2 <- matrix(fts[which(fts[,1] >= knots[2]),], ncol = 3)
      
      set2 <- matrix(set2[which(set2[,1] < knots[3]),], ncol = 3)
      
      set3 <- matrix(fts[which(fts[,1] >= knots[3]),], ncol = 3)
      
      set3 <- matrix(set3[which(set3[,1] < knots[4]),], ncol = 3)
      
      set4 <- matrix(fts[which(fts[,1] >= knots[4]),], ncol = 3)
      
      setloop <- matrix(fts[which(fts[,1] < knots[1]),], ncol = 3)
      
      setloop[,1] <- setloop[,1] + 1
      
      set4 <- rbind(set4, setloop)
      
      # Now we intepolate 99 phases to use as features.
      
      intmag[i,which(intphase < knots[1])] <- coeff4[1] +
        coeff4[2]*(intphase[which(intphase < knots[1])] - knots[4] + 1) +
        coeff4[3]*(intphase[which(intphase < knots[1])] - knots[4] + 1)^2
      
      intmag[i,which(intphase >= knots[1] & intphase < knots[2])] <- coeff1[1] +
        coeff1[2]*(intphase[which(intphase >= knots[1] & intphase < knots[2])] - knots[1]) +
        coeff1[3]*(intphase[which(intphase >= knots[1] & intphase < knots[2])] - knots[1])^2
      
      intmag[i,which(intphase >= knots[2] & intphase < knots[3])] <- coeff2[1] +
        coeff2[2]*(intphase[which(intphase >= knots[2] & intphase < knots[3])] - knots[2]) +
        coeff2[3]*(intphase[which(intphase >= knots[2] & intphase < knots[3])] - knots[2])^2
      
      intmag[i,which(intphase >= knots[3] & intphase < knots[4])] <- coeff3[1] +
        coeff3[2]*(intphase[which(intphase >= knots[3] & intphase < knots[4])] - knots[3]) +
        coeff3[3]*(intphase[which(intphase >= knots[3] & intphase < knots[4])] - knots[3])^2
      
      intmag[i,which(intphase >= knots[4])] <- coeff4[1] +
        coeff4[2]*(intphase[which(intphase >= knots[4])] - knots[4]) +
        coeff4[3]*(intphase[which(intphase >= knots[4])] - knots[4])^2
      
      # We want to normalise this between 0 and 1 and determine the normalisation factor.
      
      normfac <- abs(max(fitmag) - min(fitmag))
      
      intamp[i] <- normfac/2.0
      
      normag[i] <- min(fitmag)
      
      intmag[i,] <- ((intmag[i,] - min(fitmag)) / normfac)
      
      fts <- fts[which(fts[,3] != 1/0),]
      
      if (zeroed == TRUE){
        
        plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
             ylab = "Mean-Zeroed Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
        
      }
      else{
        
        plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
             ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
        
      }
      set1in <- matrix(c(seq(from = knots[1], to = knots[2], length.out = 1000), rep(mean(set1[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
      
      set2in <- matrix(c(seq(from = knots[2], to = knots[3], length.out = 1000), rep(mean(set2[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
      
      set3in <- matrix(c(seq(from = knots[3], to = knots[4], length.out = 1000), rep(mean(set3[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
      
      set4in <- matrix(c(seq(from = knots[4], to = (knots[1]+1), length.out = 1000), rep(mean(set4[,2]), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
      
      set1in <- rbind(set1in, set1)
      
      set1in <- set1in[order(set1in[,1]),]
      
      set2in <- rbind(set2in, set2)
      
      set2in <- set2in[order(set2in[,1]),]
      
      set3in <- rbind(set3in, set3)
      
      set3in <- set3in[order(set3in[,1]),]
      
      set4in <- rbind(set4in, set4)
      
      set4in <-set4in[order(set4in[,1]),]
      
      set4in1 <- matrix(set4in[which(set4in[,1] <= 0.5),], ncol = 3)
      
      set4in2 <- matrix(set4in[which(set4in[,1] > 0.5),], ncol = 3)
      
      preset1in <- coeff1[1] + coeff1[2]*(set1in[,1] - knots[1]) + coeff1[3]*(set1in[,1] - knots[1])^2
      
      preset2in <- coeff2[1] + coeff2[2]*(set2in[,1] - knots[2]) + coeff2[3]*(set2in[,1] - knots[2])^2
      
      preset3in <- coeff3[1] + coeff3[2]*(set3in[,1] - knots[3]) + coeff3[3]*(set3in[,1] - knots[3])^2
      
      preset4in1 <- coeff4[1] + coeff4[2]*(set4in1[,1] - knots[4]) + coeff4[3]*(set4in1[,1] - knots[4])^2
      
      preset4in2 <- coeff4[1] + coeff4[2]*(set4in2[,1] - knots[4]) + coeff4[3]*(set4in2[,1] - knots[4])^2
      
      set4in2[,1] <- set4in2[,1] - 1
      
      lines(set1in[,1], preset1in, col='red', lwd=2)
      
      lines(set2in[,1], preset2in, col='red', lwd=2)
      
      lines(set3in[,1], preset3in, col='red', lwd=2)
      
      lines(set4in1[,1], preset4in1, col='red', lwd=2)
      
      lines(set4in2[,1], preset4in2, col='red', lwd=2)
      
      #feats <- predict(lo, seq(-0.5, 0.5, length.out = len)) - weighted.mean(fts[,2], (1 / fts[,3]))
      
      print(paste("Final knots: ", knots[1], " ", knots[2], " ", knots[3], " ", knots[4], ".", sep = ""))
      
      feats <- c(allp, knots, coeff1, coeff2, coeff3, coeff4, chi2)
      
      fullfeats[i,] <- feats
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, lcamp, cordist, mnmag, sdmag, fullfeats, pzpm, pzpr, pfzp, eps, delta, log10(lambda), bins, actbins, clipped, "...", intamp, normag, intmag)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "Amplitude", "Match.Distance", "Mean.Magnitude", "StD.Magnitude", "#.Observations", "Knot.1", "Knot.2", "Knot.3", "Knot.4", "Poly.1.co.1", "Poly.1.co.2", "Poly.1.co.3",
                           "Poly.2.co.1", "Poly.2.co.2", "Poly.2.co.3", "Poly.3.co.1", "Poly.3.co.2", "Poly.3.co.3", "Poly.4.co.1", "Poly.4.co.2", "Poly.4.co.3", "Chi.Squared", "Pre.Zero.Phase.Mean", "Pre.Zero.Phase.Range", "Polyfit.Zero.Phase.Fix",
                           "Knot.Phase.Repulsion", "Knot.Median.Mag.Pull", "Log.Regularisation", "No.Bins", "Act.Bins", "Clipped", "...", "Int.Amplitude", "Norm.Magnitude", paste("V", 1:99, sep = ""))
  
  fullfeats
  
}