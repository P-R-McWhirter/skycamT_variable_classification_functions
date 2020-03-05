lcmean <- function(data, lim = 100, radii = 0.1) {
  
  start <- Sys.time()
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  mag <- as.vector(data$Magnitude)
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = 6)
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + (RA[i] - objects$RA)^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
    
    sdmag <- sd(ts[,2])
    
    bright <- median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))
    
    faint <- median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))
    
    fullfeats[i,] <- c(sdmag, meanmag, bright, faint, abs(bright-faint), min(geodist)*60*60)
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(as.vector(data$Name), as.vector(data$Type), as.vector(data$Period), as.vector(data$Magnitude), fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "AAVSO.Period", "Magnitude", "SD.Mag", "Mean.Mag", "Bright.Mag", "Faint.Mag", "Range.Mag", "Min.Distance")
  
  fullfeats
  
}