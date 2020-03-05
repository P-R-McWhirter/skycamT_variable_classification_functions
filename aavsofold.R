aavsofold <- function(data, lim = 100, spur = FALSE, pero = 1.0, errs = FALSE) {
  
  start <- Sys.time()
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  #data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
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
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    ts[,1] <- ((((ts[,1] - ts[,1][which.min(ts[,2])]) / (per[i])) + 0.25) %% 1)
    
    ts <- ts[order(ts[,1]),]
    
    fold <- ts[,1] - 1
    
    folded <- c(fold, ts[,1])
    
    foldts <- c(ts[,2], ts[,2])
    
    plot(folded, foldts, main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
         ylab = "Apparent Magnitude", pch=19, xlim = c(-1, 1), ylim = c(max(foldts), min(foldts)))
    
    magerr <- c(ts[,3], ts[,3])
    
    if (errs == TRUE){
      
      segments(folded, foldts-magerr, folded, foldts+magerr)
      
      epsilon = 0.01
      
      segments(folded-epsilon, foldts-magerr, folded+epsilon, foldts-magerr)
      
      segments(folded-epsilon, foldts+magerr, folded+epsilon, foldts+magerr)
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}







aavsofoldarea <- function(data, lim = 1, radii = 0.001, spur = FALSE, pero = 1.0, errs = FALSE) {
  
  start <- Sys.time()
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    ts <- matrix(0, nrow = 0, ncol = 3)
    
    for (j in 1:length(objects$usnoref)){
      
      info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                      objects$usnoref[j], "'", sep=""))
      
      ots <- cbind(info$MJD, info$Rcat, info$Rcaterr)
      
      ts <- rbind(ts, ots)
      
    }
    
    ts <- unique(ts)
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    ts[,1] <- ((((ts[,1] - ts[,1][which.min(ts[,2])]) / (per[i])) + 0.25) %% 1)
    
    ts <- ts[order(ts[,1]),]
    
    fold <- ts[,1] - 1
    
    folded <- c(fold, ts[,1])
    
    foldts <- c(ts[,2], ts[,2])
    
    pos <- weighted.mean(ts[,2], (1 / ts[,3])) + amplitude
    
    neg <- weighted.mean(ts[,2], (1 / ts[,3])) - amplitude
    
    plot(folded, foldts, main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
         ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
    
    magerr <- c(ts[,3], ts[,3])
    
    if (errs == TRUE){
      
      segments(folded, foldts-magerr, folded, foldts+magerr)
      
      epsilon = 0.01
      
      segments(folded-epsilon, foldts-magerr, folded+epsilon, foldts-magerr)
      
      segments(folded-epsilon, foldts+magerr, folded+epsilon, foldts+magerr)
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}





aavsofolddetrend <- function(data, lim = 100, spur = FALSE, pero = 1.0, errs = FALSE) {
  
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
    
    fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts), weights = 1/ts[,3])
    
    coeff <- coefficients(fit)
    
    ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
    
    #amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    ts[,1] <- ((((ts[,1] - ts[,1][which.min(ts[,2])]) / (per[i])) + 0.25) %% 1)
    
    ts <- ts[order(ts[,1]),]
    
    fold <- ts[,1] - 1
    
    folded <- c(fold, ts[,1])
    
    foldts <- c(ts[,2], ts[,2])
    
    plot(folded, foldts, main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0", sep=""), xlab = paste("Phase at a period of ", per[i], " days", sep=""),
         ylab = "Apparent Magnitude", pch=19)
    
    magerr <- c(ts[,3], ts[,3])
    
    if (errs == TRUE){
      
      segments(folded, foldts-magerr, folded, foldts+magerr)
      
      epsilon = 0.01
      
      segments(folded-epsilon, foldts-magerr, folded+epsilon, foldts-magerr)
      
      segments(folded-epsilon, foldts+magerr, folded+epsilon, foldts+magerr)
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}