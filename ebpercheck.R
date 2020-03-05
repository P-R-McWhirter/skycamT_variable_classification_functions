ebpercheck <- function(data, lim = 100, bins = 100, ep = 0.01, cull = 0.8, spur = FALSE, pero = 1) {
  
  start <- Sys.time()
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  if (spur == TRUE){
    
    per <- rep(pero, length(per))
     
  }
  
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
    
    tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    tsin <- rbind(tsin, ts)
    
    tsin <- tsin[order(tsin[,1]),]
    
    seastrend <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/365.24217) + cos(2*pi*tsin[,1]/365.24217) +
                      sin(4*pi*tsin[,1]/365.24217) + cos(4*pi*tsin[,1]/365.24217), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    seascoe <- seastrend$coefficients
    
    seascoe[is.na(seascoe)] <- 0
    
    LPVcut <- sqrt(as.numeric(seascoe[3])^2.0 + as.numeric(seascoe[4])^2.0)
    
    print(paste("Projected amplitude of 1 year seasonal trend: ", LPVcut, ".", sep=""))
    
    preav <- weighted.mean(ts[,2], (1 / ts[,3]))
    
    if (LPVcut < 0.5){
      
      ts[,2] <- ts[,2] - (seascoe[1] + seascoe[2]*ts[,1] + seascoe[3]*sin(2*pi*ts[,1]/365.24217) + seascoe[4]*cos(2*pi*ts[,1]/365.24217) +
                            seascoe[5]*sin(4*pi*ts[,1]/365.24217) + seascoe[6]*cos(4*pi*ts[,1]/365.24217) - preav)
      
    }
    
    
    
    
    lowper <- per[i] * (1 - ep)
    
    highper <- per[i] * (1 + ep)
    
    pers <- seq(lowper, highper, length.out = bins)
    
    fullset <- matrix(0, nrow = bins, ncol = length(ts[,1]))
    
    for (j in 1:bins){
      
      fullset[j,] <- ((((ts[,1] - ts[which.min(ts[,2]),1])/ pers[j]) + 0.5) %% 1) - 0.5
      
    }
    
    phasdif <- rep(0, length(ts[,1]))
    
    for (j in 1:length(ts[,1])){
      
      phasdif[j] <- max(fullset[,j]) - min(fullset[,j])
      
    }
    
    pos <- max(ts[,2])
    
    neg <- min(ts[,2])
    
    plot(ts[,1:2], pch=19, ylim = c(pos, neg), main = "Raw light curve", col = "red")
    
    cts <- ts[which(phasdif > cull),]
    
    points(cts[,1:2], pch=19, col = "black")
    
    fts <- ts
    
    phasmin <- fts[which.min(fts[,2]),1]
    
    fts[,1] <- ((((fts[,1] - phasmin)/ per[i]) + 0.5) %% 1) - 0.5
    
    plot(fts[,1:2], pch=19, xlim = c(-0.5, 0.5), ylim = c(pos, neg), main = "Folded light curve", col = "red")
    
    fcts <- cts
    
    fcts[,1] <- ((((fcts[,1] - phasmin)/ per[i]) + 0.5) %% 1) - 0.5
    
    points(fcts[,1:2], pch=19, col = "black")
    
    sllktest <- SLLK2(ts, from = 1/highper, to = 1/lowper, ofac = 100, plot = FALSE)
    
    sllktesteb <- SLLK2(cts, from = 1/highper, to = 1/lowper, ofac = 100)
    
    abline(v = 1/per[i], lwd = 2, col = "red")
    
    abline(v = sllktest$peak.at[1], lwd = 2, col = "darkgreen")
    
    abline(v = sllktesteb$peak.at[1], lwd = 2, col = "blue")
    
    print(paste("Reference Period: ", per[i], sep=""))
    
    print(paste("Initial period: ", sllktest$peak.at[2], seo=""))
    
    print(paste("Tuned period: ", sllktesteb$peak.at[2], sep=""))
    
    filt <- 1 - (length(cts[,1]) / length(ts[,1]))
    
    print(paste("Ratio of filtered datapoints: ", filt, sep=""))
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}