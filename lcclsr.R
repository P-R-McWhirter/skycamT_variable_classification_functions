lcclsr <- function(data, lim = 100, bins = 200) {
  
  start <- Sys.time()
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  aper <- as.numeric(as.vector(data$Period))
  
  per <- seq(from = 0.001, to = 20, length.out = 50000)
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = length(per))
  
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
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    meants <- weighted.mean(ts[,2], (1 / ts[,3]))
    
    pos <- meants + amplitude
    
    neg <- meants - amplitude
    
    ts <- ts[which(ts[,2] <= pos & ts[,2] >= neg),]
    
    feats <- rep(0, length(per))
    
    for (j in 1:length(per)){
      
      fts <- ts
      
      fts[,1] <- ((((fts[,1] - fts[,1][which.max(abs(fts[,2] - meants))]) / (per[j])) + 0.25) %% 1)
      
      if (bins < 2){
        
        print("Number of bins must be 2 or more, setting number of bins to 2.")
        
        bins <- 2
        
      }
      
      splitphase <- seq(from = 0, to = 1, length.out = bins)
      
      tery <- rep(0, (bins-1))
      
      tere <- tery
      
      tert <- tery
      
      for (k in 1:(bins-1)){
        
        tery[k] <- mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2])
        
        tere[k] <- mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3])
        
        tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
        
      }
      
      fts <- cbind(tert, tery, tere)
      
      row.has.na <- apply(fts, 1, function(x){any(is.na(x))})
      
      fts <- fts[!row.has.na,]
      
      #tsin <- matrix(c(seq(from = 0, to = 1, length.out = 20000), rep(0, times = 20000), rep(1/0, times = 20000)), nrow = 20000)
      
      #tsin <- rbind(tsin, fts)
      
      tsin <- fts
      
      #model <- lm(tsin[,2] ~ 0 + sin(2*pi*tsin[,1]) + cos(2*pi*tsin[,1]) +
      #              sin(4*pi*tsin[,1]) + cos(4*pi*tsin[,1]) + sin(6*pi*tsin[,1]) +
      #              cos(6*pi*tsin[,1]) + sin(8*pi*tsin[,1]) + cos(8*pi*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      #coeff <- model$coefficients
      
      #coeff[is.na(coeff)] <- 0
      
      Z <- cbind(sin(2*pi*tsin[,1]), cos(2*pi*tsin[,1]), sin(4*pi*tsin[,1]), cos(4*pi*tsin[,1]), sin(6*pi*tsin[,1]), cos(6*pi*tsin[,1]), sin(8*pi*tsin[,1]), cos(8*pi*tsin[,1]),
                 sin(10*pi*tsin[,1]), cos(10*pi*tsin[,1]), sin(12*pi*tsin[,1]), cos(12*pi*tsin[,1]), sin(14*pi*tsin[,1]), cos(14*pi*tsin[,1]), sin(16*pi*tsin[,1]), cos(16*pi*tsin[,1]),
                 sin(18*pi*tsin[,1]), cos(18*pi*tsin[,1]), sin(20*pi*tsin[,1]), cos(20*pi*tsin[,1]), sin(22*pi*tsin[,1]), cos(22*pi*tsin[,1]), sin(24*pi*tsin[,1]), cos(24*pi*tsin[,1]),
                 sin(26*pi*tsin[,1]), cos(26*pi*tsin[,1]), sin(28*pi*tsin[,1]), cos(28*pi*tsin[,1]), sin(30*pi*tsin[,1]), cos(30*pi*tsin[,1]), sin(32*pi*tsin[,1]), cos(32*pi*tsin[,1]),
                 sin(34*pi*tsin[,1]), cos(34*pi*tsin[,1]), sin(36*pi*tsin[,1]), cos(36*pi*tsin[,1]), sin(38*pi*tsin[,1]), cos(38*pi*tsin[,1]), sin(40*pi*tsin[,1]), cos(40*pi*tsin[,1]))
      
      Z <- as.matrix(Z)
      
      lambda <- 10^(3)
      
      W <- diag(1/tsin[,3])
      
      M <- t(Z)%*%W%*%Z + diag(rep(40,1))*lambda
      
      coeff <- solve(M)%*%t(Z)%*%W%*%(tsin[,2])
      
      A1 <- (coeff[1]^2.0 + coeff[2]^2.0)^0.5
      
      A2 <- (coeff[3]^2.0 + coeff[4]^2.0)^0.5
      
      A3 <- (coeff[5]^2.0 + coeff[6]^2.0)^0.5
      
      A4 <- (coeff[7]^2.0 + coeff[8]^2.0)^0.5
      
      A5 <- (coeff[9]^2.0 + coeff[10]^2.0)^0.5
      
      A6 <- (coeff[11]^2.0 + coeff[12]^2.0)^0.5
      
      A7 <- (coeff[13]^2.0 + coeff[14]^2.0)^0.5
      
      A8 <- (coeff[15]^2.0 + coeff[16]^2.0)^0.5
      
      A9 <- (coeff[17]^2.0 + coeff[18]^2.0)^0.5
      
      A10 <- (coeff[19]^2.0 + coeff[20]^2.0)^0.5
      
      A11 <- (coeff[21]^2.0 + coeff[22]^2.0)^0.5
      
      A12 <- (coeff[23]^2.0 + coeff[24]^2.0)^0.5
      
      A13 <- (coeff[25]^2.0 + coeff[26]^2.0)^0.5
      
      A14 <- (coeff[27]^2.0 + coeff[28]^2.0)^0.5
      
      A15 <- (coeff[29]^2.0 + coeff[30]^2.0)^0.5
      
      A16 <- (coeff[31]^2.0 + coeff[32]^2.0)^0.5
      
      A17 <- (coeff[33]^2.0 + coeff[34]^2.0)^0.5
      
      A18 <- (coeff[35]^2.0 + coeff[36]^2.0)^0.5
      
      A19 <- (coeff[37]^2.0 + coeff[38]^2.0)^0.5
      
      A20 <- (coeff[39]^2.0 + coeff[40]^2.0)^0.5
      
      A21 <- (coeff[41]^2.0 + coeff[42]^2.0)^0.5
      
      A22 <- (coeff[43]^2.0 + coeff[44]^2.0)^0.5
      
      A23 <- (coeff[45]^2.0 + coeff[46]^2.0)^0.5
      
      A24 <- (coeff[47]^2.0 + coeff[48]^2.0)^0.5
      
      A25 <- (coeff[49]^2.0 + coeff[50]^2.0)^0.5
      
      A26 <- (coeff[51]^2.0 + coeff[52]^2.0)^0.5
      
      A27 <- (coeff[53]^2.0 + coeff[54]^2.0)^0.5
      
      A28 <- (coeff[55]^2.0 + coeff[56]^2.0)^0.5
      
      A29 <- (coeff[57]^2.0 + coeff[58]^2.0)^0.5
      
      A30 <- (coeff[59]^2.0 + coeff[60]^2.0)^0.5
      
      A31 <- (coeff[61]^2.0 + coeff[62]^2.0)^0.5
      
      A32 <- (coeff[63]^2.0 + coeff[64]^2.0)^0.5
      
      A33 <- (coeff[65]^2.0 + coeff[66]^2.0)^0.5
      
      A34 <- (coeff[67]^2.0 + coeff[68]^2.0)^0.5
      
      A35 <- (coeff[69]^2.0 + coeff[70]^2.0)^0.5
      
      A36 <- (coeff[71]^2.0 + coeff[72]^2.0)^0.5
      
      A37 <- (coeff[73]^2.0 + coeff[74]^2.0)^0.5
      
      A38 <- (coeff[75]^2.0 + coeff[76]^2.0)^0.5
      
      A39 <- (coeff[77]^2.0 + coeff[78]^2.0)^0.5
      
      A40 <- (coeff[79]^2.0 + coeff[80]^2.0)^0.5
      
      amps <- c(A1, A2, A3, A4, A5, A6, A7, A8, A9, A10,
              A11, A12, A13, A14, A15, A16, A17, A18, A19, A20,
              A21, A22, A23, A24, A25, A26, A27, A28, A29, A30,
              A31, A32, A33, A34, A35, A36, A37, A38, A39, A40)
      
      amps[is.na(amps)] <- 0
      
      feats[j] <- max(amps)
        
      }
    
    fullfeats[i,] <- feats
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", paste("V", 1:length(per), sep = ""))
  
  fullfeats
  
}