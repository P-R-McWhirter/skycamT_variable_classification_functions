lcsplinefit <- function(data, lim = 100, detrend = FALSE, spur = FALSE, radii = 0.1, maxp = 20000, pero = 1.0, bins = 100, norm = TRUE, len = 100, seedno = 10, quiet = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  if (quiet == FALSE){
    
    print("Folding all objects with the AAVSO period...")
    
    print("Objects without an AAVSO period have been removed.")
    
  }
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = len)
  
  normamp <- rep(0, length(data$Name))
  
  goodness <- rep(0, length(data$Name))
  
  binnedratio <- rep(0, length(data$Name))
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
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
    
    if (quiet == FALSE){
      
      print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
      
    }
    
    if (maxp < length(ts[,1])){
      
      if (quiet == FALSE){
        
        print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
        
      }
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    fts <- ts
    
    fts[,1] <- ((((fts[,1] - fts[,1][which.max(fts[,2])]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (bins < 2){
      
      if (quiet == FALSE){
        
        print("Number of bins must be 2 or more, setting number of bins to 2.")
        
      }
      
      bins <- 2
      
    }
    
    splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
    
    tery <- rep(0, (bins-1))
    
    tere <- tery
    
    tert <- tery
    
    for (k in 1:(bins-1)){
      
      tery[k] <- median(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2])
      
      tere[k] <- median(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3])
      
      tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
      
    }
    
    fts <- cbind(tert, tery, tere)
    
    row.has.na <- apply(fts, 1, function(x){any(is.na(x))})
    
    fts <- fts[!row.has.na,,drop=F]
    
    if (nrow(fts) > 1){
      
      fts <- fts[order(fts[,1]),]
      
    }
    
    binnedratio[i] <- length(fts[,1]) / bins
    
    ftsst <- c(-0.5, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
    
    ftsend <- c(0.5, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
    
    ftsin <- rbind(ftsst, fts, ftsend)
    
    lo <- loess(ftsin[,2]~ftsin[,1], data = as.data.frame(ftsin[,1:2]), weights = 1/ftsin[,3], span = 0.5, degree = 2, model = TRUE)
    
    # Put New models here
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    predlo <- predict(lo, se = TRUE)
    
    if (quiet == FALSE){
      
      lines(ftsin[,1], predict(lo), col='red', lwd=2)
      
      lines(ftsin[,1], predlo$fit-1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      lines(ftsin[,1], predlo$fit+1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
    }
    
    feats <- predict(lo, seq(-0.5, 0.5, length.out = len)) - weighted.mean(fts[,2], (1 / fts[,3]))
    
    goodness[i] <- mean(predlo$se.fit)
    
    if (norm == TRUE){
      
      maxabsfeat <- which.max(feats)
      
      feats <- c(feats[maxabsfeat:len], feats[1:maxabsfeat-1])
      
    }
    
    normfac <- abs(max(feats) - min(feats))
    
    feats <- (feats - min(feats)) / normfac
    
    normamp[i] <- normfac / 2.0
    
    fullfeats[i,] <- feats
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, binnedratio, normamp, goodness, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "Binned.Ratio", "Amplitude", "Goodness.of.Fit", paste("V", 1:len, sep = ""))
  
  fullfeats
  
}






lcsplinesel <- function(data, lim = 100, detrend = FALSE, spur = FALSE, radii = 0.1, maxp = 20000, pero = 1.0, bins = 100, norm = TRUE, len = 100, seedno = 10, quiet = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  if (quiet == FALSE){
    
    print("Folding all objects with the AAVSO period...")
    
    print("Objects without an AAVSO period have been removed.")
    
  }
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = len)
  
  normamp <- rep(0, length(data$Name))
  
  goodness <- rep(0, length(data$Name))
  
  binnedratio <- rep(0, length(data$Name))
  
  result <- list()
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT * FROM obsdat100split WHERE usnoref = '", 
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
    
    if (quiet == FALSE){
      
      print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
      
    }
    
    if (maxp < length(ts[,1])){
      
      if (quiet == FALSE){
        
        print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
        
      }
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    #amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    fts <- ts
    
    fts[,1] <- ((((fts[,1]) / (per[i]))) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (bins < 2){
      
      if (quiet == FALSE){
        
        print("Number of bins must be 2 or more, setting number of bins to 2.")
        
      }
      
      bins <- 2
      
    }
    
    splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
    
    tery <- rep(0, (bins-1))
    
    tere <- tery
    
    tert <- tery
    
    for (k in 1:(bins-1)){
      
      tery[k] <- median(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2])
      
      tere[k] <- median(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3])
      
      tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
      
    }
    
    fts <- cbind(tert, tery, tere)
    
    row.has.na <- apply(fts, 1, function(x){any(is.na(x))})
    
    fts <- fts[!row.has.na,]
    
    fts <- fts[order(fts[,1]),]
    
    binnedratio[i] <- length(fts[,1]) / bins
    
    ftsst <- c(-0.5, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
    
    ftsend <- c(0.5, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
    
    ftsin <- rbind(ftsst, fts, ftsend)
    
    lo <- loess(ftsin[,2]~ftsin[,1], data = as.data.frame(ftsin[,1:2]), weights = 1/ftsin[,3], span = 0.5, degree = 2, model = TRUE)
    
    # Put New models here
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    predlo <- predict(lo, se = TRUE)
    
    if (quiet == FALSE){
      
      lines(ftsin[,1], predict(lo), col='red', lwd=2)
      
      lines(ftsin[,1], predlo$fit-1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      lines(ftsin[,1], predlo$fit+1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
    }
    
    feats <- predict(lo, seq(-0.5, 0.5, length.out = len)) - weighted.mean(fts[,2], (1 / fts[,3]))
    
    goodness[i] <- mean(predlo$se.fit)
    
    if (norm == TRUE){
      
      maxabsfeat <- which.max(feats)
      
      feats <- c(feats[maxabsfeat:len], feats[1:maxabsfeat-1])
      
    }
    
    normfac <- abs(max(feats) - min(feats))
    
    feats <- (feats - min(feats)) / normfac
    
    normamp[i] <- normfac
    
    fullfeats[i,] <- feats
    
    foldprox <- predict(lo, ((ts[,1] / (per[i])) %% 1)-0.5)
    
    foldprox <- abs(foldprox - ts[,2]) / normfac
    
    foldprox <- as.numeric(foldprox <= max(foldprox)/10)
    
    res <- as.data.frame(matrix(0, nrow = length(foldprox), ncol = 23))
    
    for (z in 1:length(ts[,1])){
      
      res[z,] <- info[which(info$MJD == ts[z,1]),]
      
    }
    
    res <- cbind(res, foldprox)
    
    colnames(res) <- c(colnames(info), "Fold.Prox")
    
    result[[i]] <- res
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, binnedratio, normamp, goodness, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "Binned.Ratio", "Amplitude", "Goodness.of.Fit", paste("V", 1:len, sep = ""))
  
  final <- matrix(0, ncol = 24, nrow = 0)
  
  for (i in 1:length(data$Name)){
   
    final <- rbind(final, result[[i]])
     
  }
  
  colnames(final) <- colnames(res)
  
  #final <- final[,c(24, 5, 6, 12, 14, 15, 16, 21)]
  
  final <- final[,c(24, 1:10, 12, 14:16, 19:21)]
  
  colnames(final)[1] <- "Type"
  
  final
  
}






lcpolyfit <- function(data, lim = 100, detrend = FALSE, fixint = FALSE, spur = FALSE, maxit = 5000, eps = 0.01, delta = 0.05, lambda = 0.001, radii = 0.1, maxp = 20000, pero = 1.0, seedno = 10) {
  
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





lcbinpolyfit <- function(data, lim = 100, detrend = FALSE, fixint = FALSE, spur = FALSE, maxit = 2000, eps = 0.01, delta = 0.05, lambda = 0.001, radii = 0.1, maxp = 20000, pero = 1.0, bins = 100, seedno = 10) {
  
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
    
    if (bins < 2){
      
      print("Number of bins must be 2 or more, setting number of bins to 2.")
      
      bins <- 2
      
    }
    
    splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
    
    tery <- rep(0, (bins-1))
    
    tere <- tery
    
    tert <- tery
    
    for (k in 1:(bins-1)){
      
      tery[k] <- median(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2])
      
      tere[k] <- median(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3])
      
      #tere[k] <- abs((max(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2]) - min(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2])) / 2)
      
      tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
      
    }
    
    #tere[which(tere == 0)] <- min(tere[which(tere != 0)])
    
    fts <- cbind(tert, tery, tere)
    
    row.has.na <- apply(fts, 1, function(x){any(is.na(x))})
    
    fts <- fts[!row.has.na,]
    
    fts <- fts[order(fts[,1]),]
    
    if (fixint == TRUE){
      
      meanmag <- weighted.mean(fts[,2], (1/fts[,3]))
      
      ftsin <- matrix(c(seq(from = min(fts[,1]), to = max(fts[,1]), length.out = 100), rep(meanmag, times = 100), rep(1/0, times = 100)), nrow = 100)
      
      fts <- rbind(fts, ftsin)
      
      fts <- fts[order(fts[,1]),]
      
    }
    
    knots <- sort(runif(4, min = -0.5, max = 0.5))
    
    #knots <- as.vector(c(-0.4, -0.1, 0.1, 0.4))
    
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
        
        stopchk <- 1
        
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
          else{
            
            stopchk <- stopchk + 1
            
          }
          
          if (stopchk >= 10000){
            
            stop("While loop infinite.")
            
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
      
      chi2 <- (1/bins) * (chi21 + chi22 + chi23 + chi24)
      
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




lcrawsplinefit <- function(data, lim = 100, radii = 0.1, maxp = 20000, norm = TRUE, len = 100, seedno = 10, quiet = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = len)
  
  normamp <- rep(0, length(data$Name))
  
  goodness <- rep(0, length(data$Name))
  
  binnedratio <- rep(0, length(data$Name))
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
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
    
    allp <- length(ts[,1])
    
    if (quiet == FALSE){
      
      print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
      
    }
    
    if (maxp < length(ts[,1])){
      
      if (quiet == FALSE){
        
        print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
        
      }
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    ts <- ts[order(ts[,1]),]
    
    tso <- ts
    
    tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 1000), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
    
    tsin <- rbind(tsin, ts)
    
    ts <- tsin[order(tsin[,1]),]
    
    yrran <- 365.25/(max(ts[,1] - min(ts[,1])))
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    binnedratio[i] <- 0
    
    lo <- loess(ts[,2]~ts[,1], data = as.data.frame(ts[,1:2]), weights = 1/ts[,3], span = 2*yrran, degree = 2, model = TRUE)
    
    # Put New models here
    
    if (quiet == FALSE){
      
      plot(tso[,1], tso[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Light Curve", sep=""), xlab = "Time (MJD)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), xlim = c(max(ts[,1]), min(ts[,1])))
      
    }
    
    predlo <- predict(lo, se = TRUE)
    
    if (quiet == FALSE){
      
      lines(ts[,1], predict(lo), col='red', lwd=2)
      
      #lines(ts[,1], predlo$fit-1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      #lines(ts[,1], predlo$fit+1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
    }
    
    feats <- predict(lo, seq(-0.5, 0.5, length.out = len)) - weighted.mean(ts[,2], (1 / ts[,3]))
    
    goodness[i] <- mean(predlo$se.fit)
    
    if (norm == TRUE){
      
      maxabsfeat <- which.max(feats)
      
      feats <- c(feats[maxabsfeat:len], feats[1:maxabsfeat-1])
      
    }
    
    normfac <- abs(max(feats) - min(feats))
    
    feats <- (feats - min(feats)) / normfac
    
    normamp[i] <- normfac / 2.0
    
    fullfeats[i,] <- feats
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, binnedratio, normamp, goodness, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "Binned.Ratio", "Amplitude", "Goodness.of.Fit", paste("V", 1:len, sep = ""))
  
  fullfeats
  
}




lcsysrem <- function(data, lim = 100, radii = 0.1, maxp = 20000, norm = TRUE, len = 100, seedno = 10, quiet = FALSE, peroverride = FALSE, perseas = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = 2*len)
  
  whiten <- rep(0, length(data$Name))
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
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
    
    allp <- length(ts[,1])
    
    if (quiet == FALSE){
      
      print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
      
    }
    
    if (maxp < length(ts[,1])){
      
      if (quiet == FALSE){
        
        print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
        
      }
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    ts <- ts[order(ts[,1]),]
    
    tso <- ts
    
    tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 1000), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 1000), rep(1/0, times = 1000)), nrow = 1000)
    
    tsstart <- c(min(ts[,1])-1, weighted.mean(ts[,2], (1 / ts[,3])), 0.0000001)
    
    tsend <- c(max(ts[,1])+1, weighted.mean(ts[,2], (1 / ts[,3])), 0.0000001)
    
    tsin <- rbind(tsstart, tsin, ts, tsend)
    
    ts <- tsin[order(tsin[,1]),]
    
    yrran <- 365.25/(max(ts[,1]) - min(ts[,1]))
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    lo <- loess(ts[,2]~ts[,1], data = as.data.frame(ts[,1:2]), weights = 1/ts[,3], span = 2*yrran, degree = 2, model = TRUE)
    
    # Put New models here
    
    if (is.na(per[i]) | peroverride == TRUE){
      
      per[i] <- max(ts[,1]) - min(ts[,1])
      
    }
    
    if (perseas == TRUE){
      
      per[i] <- 365.25/2
      
    }
    
    if (quiet == FALSE){
      
      plot(tso[,1], tso[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Light Curve", sep=""), xlab = "Time (MJD)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), xlim = c(min(ts[,1]), max(ts[,1])))
      
    }
    
    predlo <- predict(lo, se = TRUE)
    
    if (quiet == FALSE){
      
      lines(ts[,1], predict(lo), col='red', lwd=2)
      
    }
    
    fts <- tso
    
    maxpos <- which.max(fts[,2])
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    prefix <-  (1.0 / ((length(ts[,1]) - 1) * (sd(fts[,2])) ^ 2.0) * sum((diff(fts[,2])) ^ 2.0))
    
    mnmag <- weighted.mean(ts[,2], (1 / ts[,3]))
    
    predts <- predict(lo, tso[,1])
    
    tso[,2] <- tso[,2] - predts
    
    fts <- tso
    
    fts[,2] <- fts[,2] + mnmag
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Corrected Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    feats <- predict(lo, seq(55000, 56000, length.out = len)) - weighted.mean(ts[,2], (1 / ts[,3]))
    
    feats[which(is.na(feats))] <- 0
    
    plot(seq(55000, 56000, length.out = len), feats, type="l", ylim=c(amplitude, -amplitude))
    
    postfix <-  (1.0 / ((length(ts[,1]) - 1) * (sd(fts[,2])) ^ 2.0) * sum((diff(fts[,2])) ^ 2.0))
    
    print(paste("Pre-correction Psi-eta: ", prefix, sep = ""))
    
    print(paste("Post-correction Psi-eta: ", postfix, sep = ""))
    
    AC <- acf(feats, lag.max = len, type = "correlation", plot = T)
    
    racf <- as.numeric(AC$acf)
    
    if (postfix < prefix){
      
      whiten[i] <- 1
      
    }
    
    fullfeats[i,] <- c(feats, racf)
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, whiten, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "Whiten", paste("V", 1:len, sep = ""), paste("A", 1:len, sep = ""))
  
  fullfeats
  
}








lctrendrem <- function(data, lim = 100, radii = 0.1, maxp = 20000, sigk = 3, seedno = 10, quiet = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = 25)
  
  whiten <- rep(0, length(data$Name))
  
  whitenper <- rep(0, length(data$Name))
  
  whitenper2 <- rep(0, length(data$Name))
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    ts <- unique(ts)
    
    preav <- weighted.mean(ts[,2], (1 / ts[,3]))
    
    ts[,2] <- ts[,2] - preav
    
    tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    tsin <- rbind(tsin, ts)
    
    tsin <- tsin[order(tsin[,1]),]
    
    seastrend <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/365.24217) + cos(2*pi*tsin[,1]/365.24217) +
                      sin(4*pi*tsin[,1]/365.24217) + cos(4*pi*tsin[,1]/365.24217), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    seascoe <- seastrend$coefficients
    
    seascoe[is.na(seascoe)] <- 0
    
    seasamp <- sqrt(as.numeric(seascoe[3])^2.0 + as.numeric(seascoe[4])^2.0)
    
    seasamp2 <- sqrt(as.numeric(seascoe[5])^2.0 + as.numeric(seascoe[6])^2.0)
    
    PH1 <- atan2(seascoe[4], seascoe[3])
    
    PH2 <- atan2(seascoe[6], seascoe[5]) - 2*PH1
    
    fitparam <- c(seascoe[1], seascoe[2], seasamp, seasamp2, PH2)
    
    ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    ts <- unique(ts)
    
    ts[,2] <- ts[,2] - preav
    
    maxtime <- min(c(800, max(ts[,1]) - min(ts[,1])))
    
    if (maxtime > 80){
      
      jitamp <- 2*abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
      
      LPVsrh <- try(BGLS2(ts, plow = 80, phigh = maxtime, ofac = 5, dt = NULL, lent = NULL, spur = F, jit = jitamp, plot = F), TRUE)
      
      LPVstab <- try(is.na(max(LPVsrh$p)), TRUE)
      
      if (class(LPVsrh) == "try-error" | class(LPVstab) == "try-error" | LPVstab == TRUE){
        
        print("---------------------------------------")
        
        jitamp <- 2*jitamp
        
        LPVsrh <- try(BGLS2(ts, plow = 80, phigh = maxtime, ofac = 5, dt = NULL, lent = NULL, spur = F, jit = jitamp, plot = F), TRUE)
        
        LPVstab <- try(is.na(max(LPVsrh$p)), TRUE)
        
      }
      
      if (class(LPVsrh) == "try-error" | class(LPVstab) == "try-error" | LPVstab == TRUE){
        
        LPVguess <- 800
        
      }
      else{
        
        LPVguess <- 1 / LPVsrh$f[which.max(LPVsrh$p)]
        
      }
      
      print(LPVguess)
      
      LPVfit <- try(harm2fit(tsin, (1/LPVguess), lambda = 1e-2), TRUE)
      
      if (class(LPVfit) == "try-error"){
        
        LPVfitted <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/LPVguess) + cos(2*pi*tsin[,1]/LPVguess) +
                          sin(4*pi*tsin[,1]/LPVguess) + cos(4*pi*tsin[,1]/LPVguess), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
        LPVfit <- LPVfitted$coefficients
        
        LPVfit[is.na(LPVfit)] <- 0
        
      }
      else{
        
        LPVfitted <- list()
        
        LPVfitted$fitted.values <- LPVfit[1] + LPVfit[2]*tsin[,1] + LPVfit[3]*sin(2*pi*tsin[,1]/LPVguess) + LPVfit[4]*cos(2*pi*tsin[,1]/LPVguess) +
          LPVfit[5]*sin(4*pi*tsin[,1]/LPVguess) + LPVfit[6]*cos(4*pi*tsin[,1]/LPVguess)
        
      }
      
      PH1 <- atan2(LPVfit[4], LPVfit[3])
      
      PH2 <- atan2(LPVfit[6], LPVfit[5]) - 2*PH1
      
      trendparam <- c(LPVfit[1], LPVfit[2], sqrt(LPVfit[3]^2 + LPVfit[4]^2), sqrt(LPVfit[5]^2 + LPVfit[6]^2), PH2)
      
      LPVamp <- sqrt(as.numeric(LPVfit[3])^2.0 + as.numeric(LPVfit[4])^2.0)
      
      LPVamp2 <- sqrt(as.numeric(LPVfit[5])^2.0 + as.numeric(LPVfit[6])^2.0)
      
      maxval <- max(LPVfitted$fitted.values) - max(seastrend$fitted.values)
      
      minval <- min(LPVfitted$fitted.values) - min(seastrend$fitted.values)
      
      print(c(seasamp, seasamp2, LPVamp, LPVamp2))
      
      ampsim <- as.vector(abs(outer(c(seasamp, seasamp2), c(LPVamp, LPVamp2), "-")))
      
      print(ampsim)
      
      ampratio <- seasamp2 / seasamp
      
      ampratio2 <- LPVamp2 / LPVamp
      
      print(c(ampratio, ampratio2))
      
      ampcut <- ampratio / ampratio2
      
      print(ampcut)
      
      print(any(ampsim < 0.05))
      
      print(c(seasamp, seasamp2))
      
      trendcor <- cor(seastrend$fitted.values, LPVfitted$fitted.values)
      
      print(trendcor)
      
    }
    
    keep <- sigclip(info$Rcat, sig = sigk, tol = 0.000001)
    
    info <- info[keep,]
    
    ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    prelen <- length(ts[,2])
    
    ts <- unique(ts)
    
    plot(ts[,1:2], pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
    lines(tsin[,1], seastrend$fitted.values+preav, col = "red")
    
    if (maxtime > 80){
      
      lines(tsin[,1], LPVfitted$fitted.values+preav, col = "blue")
      
    }
    
    tso <- ts
    
    fts <- tso
    
    maxpos <- which.max(fts[,2])
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    prefix <-  (1.0 / ((length(ts[,1]) - 1) * (sd(fts[,2])) ^ 2.0) * sum((diff(fts[,2])) ^ 2.0))
    
    mnmag <- weighted.mean(tso[,2], (1 / tso[,3]))
    
    tso[,2] <- tso[,2] - (seascoe[1] + seascoe[2]*tso[,1] + seascoe[3]*sin(2*pi*tso[,1]/365.24217) + seascoe[4]*cos(2*pi*tso[,1]/365.24217) +
                            seascoe[5]*sin(4*pi*tso[,1]/365.24217) + seascoe[6]*cos(4*pi*tso[,1]/365.24217))
    
    fts <- tso
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Corrected Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    postfix <-  (1.0 / ((length(ts[,1]) - 1) * (sd(fts[,2])) ^ 2.0) * sum((diff(fts[,2])) ^ 2.0))
    
    print(paste("Pre-correction Psi-eta: ", prefix, sep = ""))
    
    print(paste("Post-correction Psi-eta: ", postfix, sep = ""))
    
    t <- ts[,1]
    y <- ts[,2]-mean(ts[,2])
    norm <- 1/(2 * var(y))
    
    wi <- 2 * pi * (1/per[i])
    tau <- 0.5 * atan2(sum(sin(wi * t)), sum(cos(wi * t)))/wi
    arg <- wi * (t - tau)
    cs <- cos(arg)
    sn <- sin(arg)
    A <- (sum(y * cs))^2
    B <- sum(cs * cs)
    C <- (sum(y * sn))^2
    D <- sum(sn * sn)
    prep <- norm*(A/B + C/D)
    
    
    t <- tso[,1]
    y <- tso[,2]-mean(tso[,2])
    norm <- 1/(2 * var(y))
    
    
    tau <- 0.5 * atan2(sum(sin(wi * t)), sum(cos(wi * t)))/wi
    arg <- wi * (t - tau)
    cs <- cos(arg)
    sn <- sin(arg)
    A <- (sum(y * cs))^2
    B <- sum(cs * cs)
    C <- (sum(y * sn))^2
    D <- sum(sn * sn)
    postp <- norm*(A/B + C/D)
    
    
    print(paste("Pre-correction P: ", prep, sep = ""))
    
    print(paste("Post-correction P: ", postp, sep = ""))
    
    if (postfix < prefix){
      
      whiten[i] <- 1
      
    }
    
    if (postp > prep){
      
      whitenper[i] <- 1
      
    }
    
    feats <- c(LPVguess, trendparam, fitparam, LPVamp, LPVamp2, maxval, minval, ampsim, ampratio, ampratio2, ampcut, trendcor, prefix, postfix)
    
    fullfeats[i,] <- feats
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, whiten, whitenper, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", "Whiten", "Whitenper", paste("V", 1:25, sep = ""))
  
  fullfeats
  
}














lcsysfit <- function(data, trend, trcut, lim = 100, radii = 0.1, maxp = 20000, seedno = 10, model = "spline", quiet = FALSE, peroverride = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  trendmod <- cbind(seq(54884,56079, length.out = length(trend$lctrendmod)), trend$lctrendmod, rep(1, length(trend$lctrendmod)))
  
  trendmod[c(1,nrow(trendmod)),3] <- 1/0
  
  if (model == "spline"){
    
    yrran <- 365.25/(max(trendmod[,1]) - min(trendmod[,1]))
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- loess(trendmod[,2]~trendmod[,1], data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3], span = 0.25*yrran, degree = 2, model = TRUE)
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- predict(lo, trin[,1])
    
    predlo <- predict(lo, se = TRUE)
    
  }
  else if (model == "harm"){
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- lm(trendmod[,2] ~ trendmod[,1] + sin(2*pi*trendmod[,1]/365.24217) + cos(2*pi*trendmod[,1]/365.24217) +
               sin(4*pi*trendmod[,1]/365.24217) + cos(4*pi*trendmod[,1]/365.24217), data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3])
    
    lo <- lo$coefficients
    
    lo[is.na(lo)] <- 0
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- lo[1] + lo[2]*trin[,1] + lo[3]*sin(2*pi*trin[,1]/365.24217) + lo[4]*cos(2*pi*trin[,1]/365.24217) +
      lo[5]*sin(4*pi*trin[,1]/365.24217) + lo[6]*cos(4*pi*trin[,1]/365.24217)
    
    predlo <- lo[1] + lo[2]*trendmod[,1] + lo[3]*sin(2*pi*trendmod[,1]/365.24217) + lo[4]*cos(2*pi*trendmod[,1]/365.24217) +
      lo[5]*sin(4*pi*trendmod[,1]/365.24217) + lo[6]*cos(4*pi*trendmod[,1]/365.24217)
    
  }
  else{
    
    print(paste("'", model, "' is an unrecognised model. Exiting software.", sep = ""))
    
    stop("Invalid model")
    
  }
  
  if (quiet == FALSE){
    
    plot(trendmod[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(0.5, -0.5))
    
    lines(trin[,1], prelo, col='red', lwd=2)
    
  }
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
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
    
    ts[which(ts[,1] < 55200),2] <- ts[which(ts[,1] < 55200),2] - 0.1
    
    allp <- length(ts[,1])
    
    if (quiet == FALSE){
      
      print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
      
    }
    
    if (maxp < length(ts[,1])){
      
      if (quiet == FALSE){
        
        print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
        
      }
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    ts <- ts[order(ts[,1]),]
    
    tso <- ts
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (is.na(per[i]) | peroverride == TRUE){
      
      per[i] <- max(ts[,1]) - min(ts[,1])
      
    }
    
    if (quiet == FALSE){
      
      plot(tso[,1], tso[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Light Curve", sep=""), xlab = "Time (MJD)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), xlim = c(min(ts[,1]), max(ts[,1])))
      
    }
    
    fts <- tso
    
    maxpos <- which.max(fts[,2])
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if(model == "spline"){
      
      predts <- predict(lo, tso[,1])
      
    }
    else if (model == "harm"){
      
      predts <- lo[1] + lo[2]*tso[,1] + lo[3]*sin(2*pi*tso[,1]/365.24217) + lo[4]*cos(2*pi*tso[,1]/365.24217) +
        lo[5]*sin(4*pi*tso[,1]/365.24217) + lo[6]*cos(4*pi*tso[,1]/365.24217)
      
    }
    
    tso[,2] <- tso[,2] - predts
    
    #lsptest <- lsp(tso[,1:2], from = 0.001, to = 10, type = "frequency", ofac = 10)
    
    #print(1/lsptest$scanned[order(lsptest$power, decreasing = T)][1:100])
    
    if (quiet == FALSE){
      
      plot(tso[,1], tso[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Corrected Light Curve", sep=""), xlab = "Time (MJD)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), xlim = c(min(ts[,1]), max(ts[,1])))
      
    }
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    fts <- tso
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Corrected Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  lo
  
}




lcsysfit2 <- function(data, trend, trcut, lim = 100, radii = 0.1, maxp = 20000, seedno = 10, model = "spline", quiet = FALSE, peroverride = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  trendmod <- cbind(seq(54884,56079, length.out = length(trend$lctrendmod)), trend$lctrendmod, rep(1, length(trend$lctrendmod)))
  
  trendmod[c(1,nrow(trendmod)),3] <- 1/0
  
  if (model == "spline"){
    
    yrran <- 365.25/(max(trendmod[,1]) - min(trendmod[,1]))
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- loess(trendmod[,2]~trendmod[,1], data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3], span = 0.25*yrran, degree = 2, model = TRUE)
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- predict(lo, trin[,1])
    
    predlo <- predict(lo, se = TRUE)
    
  }
  else if (model == "harm"){
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- lm(trendmod[,2] ~ trendmod[,1] + sin(2*pi*trendmod[,1]/365.24217) + cos(2*pi*trendmod[,1]/365.24217), data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3])
    
    lo <- lo$coefficients
    
    lo[is.na(lo)] <- 0
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- lo[1] + lo[2]*trin[,1] + lo[3]*sin(2*pi*trin[,1]/365.24217) + lo[4]*cos(2*pi*trin[,1]/365.24217)
    
    predlo <- lo[1] + lo[2]*trendmod[,1] + lo[3]*sin(2*pi*trendmod[,1]/365.24217) + lo[4]*cos(2*pi*trendmod[,1]/365.24217)
    
  }
  else{
    
    print(paste("'", model, "' is an unrecognised model. Exiting software.", sep = ""))
    
    stop("Invalid model")
    
  }
  
  if (quiet == FALSE){
    
    plot(trendmod[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(0.5, -0.5))
    
    lines(trin[,1], prelo, col='red', lwd=2)
    
  }
  
  for (i in 1:length(data$Name)){
    
    if (quiet == FALSE){
      
      print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
      
    }
    
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
    
    allp <- length(ts[,1])
    
    if (quiet == FALSE){
      
      print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
      
    }
    
    if (maxp < length(ts[,1])){
      
      if (quiet == FALSE){
        
        print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
        
      }
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    ts <- ts[order(ts[,1]),]
    
    tso <- ts
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (is.na(per[i]) | peroverride == TRUE){
      
      per[i] <- max(ts[,1]) - min(ts[,1])
      
    }
    
    if (quiet == FALSE){
      
      plot(tso[,1], tso[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Light Curve", sep=""), xlab = "Time (MJD)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), xlim = c(min(ts[,1]), max(ts[,1])))
      
    }
    
    fts <- tso
    
    maxpos <- which.max(fts[,2])
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if(model == "spline"){
      
      predts <- predict(lo, tso[,1])
      
    }
    else if (model == "harm"){
      
      predts <- lo[1] + lo[2]*tso[,1] + lo[3]*sin(2*pi*tso[,1]/365.24217) + lo[4]*cos(2*pi*tso[,1]/365.24217)
      
    }
    
    tso[,2] <- tso[,2] - predts
    
    lsptest <- lsp(tso[,1:2], from = 0.001, to = 10, type = "frequency", ofac = 10)
    
    print(1/lsptest$scanned[order(lsptest$power, decreasing = T)][1:100])
    
    if (quiet == FALSE){
      
      plot(tso[,1], tso[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Corrected Light Curve", sep=""), xlab = "Time (MJD)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), xlim = c(min(ts[,1]), max(ts[,1])))
      
    }
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
    fts <- tso
    
    fts[,1] <- ((((fts[,1] - fts[,1][maxpos]) / (per[i])) + 0.75) %% 1) - 0.5
    
    fts <- fts[order(fts[,1]),]
    
    if (quiet == FALSE){
      
      plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Corrected Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
    }
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  print(c(as.numeric(sqrt(lo[3]^2 + lo[4]^2)), as.numeric(atan2(lo[4], lo[3]))))
  
  lo
  
}



lcsysmodplot <- function(trend, trcut, seedno = 10, model = "spline", quiet = FALSE) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  trendmod <- cbind(seq(54884,56079, length.out = length(trend$lctrendmod)), trend$lctrendmod, rep(1, length(trend$lctrendmod)))
  
  trendmod[c(1,nrow(trendmod)),3] <- 1/0
  
  if (model == "spline"){
    
    yrran <- 365.25/(max(trendmod[,1]) - min(trendmod[,1]))
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- loess(trendmod[,2]~trendmod[,1], data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3], span = 0.25*yrran, degree = 2, model = TRUE)
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- predict(lo, trin[,1])
    
    predlo <- predict(lo, se = TRUE)
    
  }
  else if (model == "harm"){
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- lm(trendmod[,2] ~ sin(2*pi*trendmod[,1]/365.24217) + cos(2*pi*trendmod[,1]/365.24217), data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3])
    
    lo <- lo$coefficients
    
    lo[is.na(lo)] <- 0
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- lo[1] + lo[2]*sin(2*pi*trin[,1]/365.24217) + lo[3]*cos(2*pi*trin[,1]/365.24217)
    
    predlo <- lo[1] + lo[2]*sin(2*pi*trendmod[,1]/365.24217) + lo[3]*cos(2*pi*trendmod[,1]/365.24217)
    
  }
  else{
    
    print(paste("'", model, "' is an unrecognised model. Exiting software.", sep = ""))
    
    stop("Invalid model")
    
  }
  
  if (quiet == FALSE){
    
    plot(trendmod[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(0.5, -0.5))
    
    lines(trin[,1], prelo, col='red', lwd=2)
    
  }
  
  if (quiet == FALSE){
    
    print("Operation Complete.")
    
    print(Sys.time() - start)
    
  }
  
  close(channel)
  
  amp <- c(as.numeric(sqrt(lo[2]^2 + lo[3]^2)), as.numeric(atan2(lo[3], lo[2])))
  
  sp.out <- list(lo = lo, amp = amp)
  
  sp.out
  
}



lcsysampinfo <- function(n = 10){
  
  p1a <- lcsysmodplot(trend = list(lctrendmod = lctrend[1,], lcnum = lcnum[1,]), n, model = "harm")$amp
  
  p2a <- lcsysmodplot(trend = list(lctrendmod = lctrend[2,], lcnum = lcnum[2,]), n, model = "harm")$amp
  
  p3a <- lcsysmodplot(trend = list(lctrendmod = lctrend[3,], lcnum = lcnum[3,]), n, model = "harm")$amp
  
  p4a <- lcsysmodplot(trend = list(lctrendmod = lctrend[4,], lcnum = lcnum[4,]), n, model = "harm")$amp
  
  p5a <- lcsysmodplot(trend = list(lctrendmod = lctrend[5,], lcnum = lcnum[5,]), n, model = "harm")$amp
  
  p6a <- lcsysmodplot(trend = list(lctrendmod = lctrend[6,], lcnum = lcnum[6,]), n, model = "harm")$amp
  
  p7a <- lcsysmodplot(trend = list(lctrendmod = lctrend[7,], lcnum = lcnum[7,]), n, model = "harm")$amp
  
  p8a <- lcsysmodplot(trend = list(lctrendmod = lctrend[8,], lcnum = lcnum[8,]), n, model = "harm")$amp
  
  p1b <- lcsysmodplot(trend = list(lctrendmod = lctrend[9,], lcnum = lcnum[9,]), n, model = "harm")$amp
  
  p2b <- lcsysmodplot(trend = list(lctrendmod = lctrend[10,], lcnum = lcnum[10,]), n, model = "harm")$amp
  
  p3b <- lcsysmodplot(trend = list(lctrendmod = lctrend[11,], lcnum = lcnum[11,]), n, model = "harm")$amp
  
  p4b <- lcsysmodplot(trend = list(lctrendmod = lctrend[12,], lcnum = lcnum[12,]), n, model = "harm")$amp
  
  p5b <- lcsysmodplot(trend = list(lctrendmod = lctrend[13,], lcnum = lcnum[13,]), n, model = "harm")$amp
  
  p6b <- lcsysmodplot(trend = list(lctrendmod = lctrend[14,], lcnum = lcnum[14,]), n, model = "harm")$amp
  
  p7b <- lcsysmodplot(trend = list(lctrendmod = lctrend[15,], lcnum = lcnum[15,]), n, model = "harm")$amp
  
  p8b <- lcsysmodplot(trend = list(lctrendmod = lctrend[16,], lcnum = lcnum[16,]), n, model = "harm")$amp
  
  p1c <- lcsysmodplot(trend = list(lctrendmod = lctrend[17,], lcnum = lcnum[17,]), n, model = "harm")$amp
  
  p2c <- lcsysmodplot(trend = list(lctrendmod = lctrend[18,], lcnum = lcnum[18,]), n, model = "harm")$amp
  
  p3c <- lcsysmodplot(trend = list(lctrendmod = lctrend[19,], lcnum = lcnum[19,]), n, model = "harm")$amp
  
  p4c <- lcsysmodplot(trend = list(lctrendmod = lctrend[20,], lcnum = lcnum[20,]), n, model = "harm")$amp
  
  p5c <- lcsysmodplot(trend = list(lctrendmod = lctrend[21,], lcnum = lcnum[21,]), n, model = "harm")$amp
  
  p6c <- lcsysmodplot(trend = list(lctrendmod = lctrend[22,], lcnum = lcnum[22,]), n, model = "harm")$amp
  
  p7c <- lcsysmodplot(trend = list(lctrendmod = lctrend[23,], lcnum = lcnum[23,]), n, model = "harm")$amp
  
  p8c <- lcsysmodplot(trend = list(lctrendmod = lctrend[24,], lcnum = lcnum[24,]), n, model = "harm")$amp
  
  p1d <- lcsysmodplot(trend = list(lctrendmod = lctrend[25,], lcnum = lcnum[25,]), n, model = "harm")$amp
  
  p2d <- lcsysmodplot(trend = list(lctrendmod = lctrend[26,], lcnum = lcnum[26,]), n, model = "harm")$amp
  
  p3d <- lcsysmodplot(trend = list(lctrendmod = lctrend[27,], lcnum = lcnum[27,]), n, model = "harm")$amp
  
  p4d <- lcsysmodplot(trend = list(lctrendmod = lctrend[28,], lcnum = lcnum[28,]), n, model = "harm")$amp
  
  p5d <- lcsysmodplot(trend = list(lctrendmod = lctrend[29,], lcnum = lcnum[29,]), n, model = "harm")$amp
  
  p6d <- lcsysmodplot(trend = list(lctrendmod = lctrend[30,], lcnum = lcnum[30,]), n, model = "harm")$amp
  
  p7d <- lcsysmodplot(trend = list(lctrendmod = lctrend[31,], lcnum = lcnum[31,]), n, model = "harm")$amp
  
  p8d <- lcsysmodplot(trend = list(lctrendmod = lctrend[32,], lcnum = lcnum[32,]), n, model = "harm")$amp
  
  p1e <- lcsysmodplot(trend = list(lctrendmod = lctrend[33,], lcnum = lcnum[33,]), n, model = "harm")$amp
  
  p2e <- lcsysmodplot(trend = list(lctrendmod = lctrend[34,], lcnum = lcnum[34,]), n, model = "harm")$amp
  
  p3e <- lcsysmodplot(trend = list(lctrendmod = lctrend[35,], lcnum = lcnum[35,]), n, model = "harm")$amp
  
  p4e <- lcsysmodplot(trend = list(lctrendmod = lctrend[36,], lcnum = lcnum[36,]), n, model = "harm")$amp
  
  p5e <- lcsysmodplot(trend = list(lctrendmod = lctrend[37,], lcnum = lcnum[37,]), n, model = "harm")$amp
  
  p6e <- lcsysmodplot(trend = list(lctrendmod = lctrend[38,], lcnum = lcnum[38,]), n, model = "harm")$amp
  
  p7e <- lcsysmodplot(trend = list(lctrendmod = lctrend[39,], lcnum = lcnum[39,]), n, model = "harm")$amp
  
  p8e <- lcsysmodplot(trend = list(lctrendmod = lctrend[40,], lcnum = lcnum[40,]), n, model = "harm")$amp
  
  p1f <- lcsysmodplot(trend = list(lctrendmod = lctrend[41,], lcnum = lcnum[41,]), n, model = "harm")$amp
  
  p2f <- lcsysmodplot(trend = list(lctrendmod = lctrend[42,], lcnum = lcnum[42,]), n, model = "harm")$amp
  
  p3f <- lcsysmodplot(trend = list(lctrendmod = lctrend[43,], lcnum = lcnum[43,]), n, model = "harm")$amp
  
  p4f <- lcsysmodplot(trend = list(lctrendmod = lctrend[44,], lcnum = lcnum[44,]), n, model = "harm")$amp
  
  p5f <- lcsysmodplot(trend = list(lctrendmod = lctrend[45,], lcnum = lcnum[45,]), n, model = "harm")$amp
  
  p6f <- lcsysmodplot(trend = list(lctrendmod = lctrend[46,], lcnum = lcnum[46,]), n, model = "harm")$amp
  
  p7f <- lcsysmodplot(trend = list(lctrendmod = lctrend[47,], lcnum = lcnum[47,]), n, model = "harm")$amp
  
  p8f <- lcsysmodplot(trend = list(lctrendmod = lctrend[48,], lcnum = lcnum[48,]), n, model = "harm")$amp
  
  p1g <- lcsysmodplot(trend = list(lctrendmod = lctrend[49,], lcnum = lcnum[49,]), n, model = "harm")$amp
  
  p2g <- lcsysmodplot(trend = list(lctrendmod = lctrend[50,], lcnum = lcnum[50,]), n, model = "harm")$amp
  
  p3g <- lcsysmodplot(trend = list(lctrendmod = lctrend[51,], lcnum = lcnum[51,]), n, model = "harm")$amp
  
  p4g <- lcsysmodplot(trend = list(lctrendmod = lctrend[52,], lcnum = lcnum[52,]), n, model = "harm")$amp
  
  p5g <- lcsysmodplot(trend = list(lctrendmod = lctrend[53,], lcnum = lcnum[53,]), n, model = "harm")$amp
  
  p6g <- lcsysmodplot(trend = list(lctrendmod = lctrend[54,], lcnum = lcnum[54,]), n, model = "harm")$amp
  
  p7g <- lcsysmodplot(trend = list(lctrendmod = lctrend[55,], lcnum = lcnum[55,]), n, model = "harm")$amp
  
  p8g <- lcsysmodplot(trend = list(lctrendmod = lctrend[56,], lcnum = lcnum[56,]), n, model = "harm")$amp
  
  p1h <- lcsysmodplot(trend = list(lctrendmod = lctrend[57,], lcnum = lcnum[57,]), n, model = "harm")$amp
  
  p2h <- lcsysmodplot(trend = list(lctrendmod = lctrend[58,], lcnum = lcnum[58,]), n, model = "harm")$amp
  
  p3h <- lcsysmodplot(trend = list(lctrendmod = lctrend[59,], lcnum = lcnum[59,]), n, model = "harm")$amp
  
  p4h <- lcsysmodplot(trend = list(lctrendmod = lctrend[60,], lcnum = lcnum[60,]), n, model = "harm")$amp
  
  p5h <- lcsysmodplot(trend = list(lctrendmod = lctrend[61,], lcnum = lcnum[61,]), n, model = "harm")$amp
  
  p6h <- lcsysmodplot(trend = list(lctrendmod = lctrend[62,], lcnum = lcnum[62,]), n, model = "harm")$amp
  
  p7h <- lcsysmodplot(trend = list(lctrendmod = lctrend[63,], lcnum = lcnum[63,]), n, model = "harm")$amp
  
  p8h <- lcsysmodplot(trend = list(lctrendmod = lctrend[64,], lcnum = lcnum[64,]), n, model = "harm")$amp
  
  z <- matrix(c(p1a[1], p2a[1], p3a[1], p4a[1], p5a[1], p6a[1], p7a[1], p8a[1], p1b[1], p2b[1], p3b[1], p4b[1], p5b[1], p6b[1], p7b[1], p8b[1], 
                   p1c[1], p2c[1], p3c[1], p4c[1], p5c[1], p6c[1], p7c[1], p8c[1], p1d[1], p2d[1], p3d[1], p4d[1], p5d[1], p6d[1], p7d[1], p8d[1], 
                   p1e[1], p2e[1], p3e[1], p4e[1], p5e[1], p6e[1], p7e[1], p8e[1], p1f[1], p2f[1], p3f[1], p4f[1], p5f[1], p6f[1], p7f[1], p8f[1], 
                   p1g[1], p2g[1], p3g[1], p4g[1], p5g[1], p6g[1], p7g[1], p8g[1], p1h[1], p2h[1], p3h[1], p4h[1], p5h[1], p6h[1], p7h[1], p8h[1]), nrow = 8, ncol = 8)
  
  print(z)
  
  decs <- colMedians(z)
  
  ras <- rowMedians(z)
  
  persp(seq(22.5, 337.5, length.out = 8), seq(-32.5, 72.5, length.out = 8), z, phi = 50, theta = 40,
        xlab = "Right Ascension", ylab = "Declination", zlab = "Trend Amplitude (mag)",
        main = "Trend Amplitude"
  )
  
  filled.contour(z, y = seq(-32.5, 72.5, length.out = 8), x = seq(22.5, 337.5, length.out = 8), col=colorpanel(27, "white", "grey10"), main = "Contour plot of the trend amplitude", ylab = "Declination (degrees)", xlab = "Right Ascension (degrees)")
  
  plot(seq(22.5, 337.5, length.out = 8), ras, type= "l")
  
  plot(seq(-32.5, 72.5, length.out = 8), decs, type= "l")
  
  z <- matrix(c(p1a[2], p2a[2], p3a[2], p4a[2], p5a[2], p6a[2], p7a[2], p8a[2], p1b[2], p2b[2], p3b[2], p4b[2], p5b[2], p6b[2], p7b[2], p8b[2], 
                p1c[2], p2c[2], p3c[2], p4c[2], p5c[2], p6c[2], p7c[2], p8c[2], p1d[2], p2d[2], p3d[2], p4d[2], p5d[2], p6d[2], p7d[2], p8d[2], 
                p1e[2], p2e[2], p3e[2], p4e[2], p5e[2], p6e[2], p7e[2], p8e[2], p1f[2], p2f[2], p3f[2], p4f[2], p5f[2], p6f[2], p7f[2], p8f[2], 
                p1g[2], p2g[2], p3g[2], p4g[2], p5g[2], p6g[2], p7g[2], p8g[2], p1h[2], p2h[2], p3h[2], p4h[2], p5h[2], p6h[2], p7h[2], p8h[2]), nrow = 8, ncol = 8)
  
  print(z)
  
  decs <- colMedians(z)
  
  ras <- rowMedians(z)
  
  persp(seq(22.5, 337.5, length.out = 8), seq(-32.5, 72.5, length.out = 8), z, phi = 50, theta = 40,
        xlab = "Right Ascension", ylab = "Declination", zlab = "Trend Phase",
        main = "Trend Phase"
  )
  
  filled.contour(z, y = seq(-32.5, 72.5, length.out = 8), x = seq(22.5, 337.5, length.out = 8), col=colorpanel(14, "white", "grey10"), main = "Contour plot of the trend phase", ylab = "Declination (degrees)", xlab = "Right Ascension (degrees)")
  
  plot(seq(22.5, 337.5, length.out = 8), ras, type= "l", xlim = c(0, 360), lwd = 2, main = "Plot of Trend Phase against Right Ascension", ylab = "Trend Phase", xlab = "Right Ascension (degrees)")
  
  plot(seq(-32.5, 72.5, length.out = 8), decs, type= "l")
  
}







lcrpartfit2 <- function(ts, quiet = FALSE, cut = 0.2) {
  
  library(rpart)
  
  ts <- ts[order(ts[,1]),]
  
  preav <- weighted.mean(ts[,2], (1 / ts[,3]))
  
  ts[,2] <- ts[,2] - preav
  
  diffts <- diff(ts[,1])
  
  rangts <- max(ts[,1]) - min(ts[,1])
  
  parts <- ceiling(rangts/365.24217)-1
  
  copart <- length(diffts[which(diffts >= max(diffts)/(parts+1))])
  
  if (copart == parts | copart == parts+1){
    
    split <- sort(c(1, order(diffts, decreasing = T)[1:copart], length(ts[,1])))
    
  }
  else{
  
    split <- sort(c(1, order(diffts, decreasing = T)[1:parts], length(ts[,1])))
  
  }
  
  tr <- rep(0, length(split)-1)
  
  for (i in 1:length(tr)){
    
    tr[i] <- weighted.mean(ts[split[i]:split[i+1],2], (1 / ts[split[i]:split[i+1],3]))
    
  }
  
  if (quiet == FALSE){
    
    plot(ts[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
    for (i in 1:length(tr)){
    
      lines(ts[split[i]:split[i+1],1], rep(tr[i], length(ts[split[i]:split[i+1],1])), col='red', lwd=2)
      
    }
    
  }
  
  sepamt <- max(tr) - min(tr)
  
  if (sepamt <= cut){
  
    for (i in 1:length(tr)){
    
      ts[split[i]:split[i+1],2] <- ts[split[i]:split[i+1],2] - tr[i]
    
    }
    
  }
  
  if (quiet == FALSE){
    
    plot(ts[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
  }
  
  out <- list(ts = ts, predlo = predlo, lo = lo)
  
  out
  
}





lcrpartfit <- function(ts, quiet = FALSE, cut = 0.2) {
  
  library(rpart)
  
  ts <- ts[order(ts[,1]),]
  
  lo <- rpart(ts[,2]~ts[,1], data = as.data.frame(ts[,1:2]), weights = 1/ts[,3], model = F, control=rpart.control(minsplit=length(ts[,1])/3, cp=0.1))
  
  predlo <- predict(lo, se = TRUE)
  
  if (quiet == FALSE){
    
    plot(ts[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
    lines(ts[,1], predlo, col='red', lwd=2)
    
  }
  
  preav <- weighted.mean(ts[,2], (1 / ts[,3]))
  
  rang <- max(predlo) - min(predlo)
  
  if (rang <= cut){
    
    ts[,2] <- ts[,2] - predlo + preav
    
  }
  
  out <- list(ts = ts, predlo = predlo, lo = lo)
  
  out
  
}








lcsysco <- function(RA, DEC, trend, trcut, model = "spline", quiet = FALSE) {
  
  rasplit <- seq(0, 360, length.out = 9)
  
  decsplit <- seq(-40, 80, length.out = 9)
  
  raloc <- 1
  
  decloc <- 1
  
  for (i in 1:(length(rasplit)-1)){
    
    if (RA >= rasplit[i] & RA < rasplit[i+1]){
      
      raloc <- i
      
    }
    else if(RA >= rasplit[i+1]){
      
      raloc <- length(rasplit)-1
      
    }
    
    if (DEC >= decsplit[i] & DEC < decsplit[i+1]){
      
      decloc <- i
      
    }
    else if(DEC >= decsplit[i+1]){
      
      decloc <- length(rasplit)-1
      
    }
    
  }
  
  rowloc <- matrix(1:64, nrow = 8, ncol = 8)[raloc,decloc]
  
  trend$lctrendmod <- trend$lctrendmod[rowloc,]
  
  trend$lcnum <- trend$lcnum[rowloc,]
  
  trendmod <- cbind(seq(54884,56079, length.out = length(trend$lctrendmod)), trend$lctrendmod, rep(1, length(trend$lctrendmod)))
  
  trendmod[c(1,nrow(trendmod)),3] <- 1/0
  
  if (model == "spline"){
    
    yrran <- 365.25/(max(trendmod[,1]) - min(trendmod[,1]))
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- loess(trendmod[,2]~trendmod[,1], data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3], span = 0.25*yrran, degree = 2, model = TRUE)
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- predict(lo, trin[,1])
    
    predlo <- predict(lo, se = TRUE)
    
  }
  else if (model == "harm"){
    
    trendmod <- trendmod[which(trend$lcnum >= trcut),]
    
    lo <- lm(trendmod[,2] ~ trendmod[,1] + sin(2*pi*trendmod[,1]/365.24217) + cos(2*pi*trendmod[,1]/365.24217) +
               sin(4*pi*trendmod[,1]/365.24217) + cos(4*pi*trendmod[,1]/365.24217), data = as.data.frame(trendmod[,1:2]), weights = 1/trendmod[,3])
    
    lo <- lo$coefficients
    
    lo[is.na(lo)] <- 0
    
    trin <- matrix(c(seq(from = min(trendmod[,1]), to = max(trendmod[,1]), length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    trin <- rbind(trin, trendmod)
    
    trin <- trin[order(trin[,1]),]
    
    prelo <- lo[1] + lo[2]*trin[,1] + lo[3]*sin(2*pi*trin[,1]/365.24217) + lo[4]*cos(2*pi*trin[,1]/365.24217) +
      lo[5]*sin(4*pi*trin[,1]/365.24217) + lo[6]*cos(4*pi*trin[,1]/365.24217)
    
    predlo <- lo[1] + lo[2]*trendmod[,1] + lo[3]*sin(2*pi*trendmod[,1]/365.24217) + lo[4]*cos(2*pi*trendmod[,1]/365.24217) +
      lo[5]*sin(4*pi*trendmod[,1]/365.24217) + lo[6]*cos(4*pi*trendmod[,1]/365.24217)
    
  }
  else{
    
    print(paste("'", model, "' is an unrecognised model. Exiting software.", sep = ""))
    
    stop("Invalid model")
    
  }
  
  if (quiet == FALSE){
    
    plot(trendmod[,1:2], main = "Trend Light Curve", xlab = "Modified Julian Date", ylab = "Apparent Magnitude", pch=19, ylim = c(0.5, -0.5))
    
    lines(trin[,1], prelo, col='red', lwd=2)
    
  }
  
  lo
  
}