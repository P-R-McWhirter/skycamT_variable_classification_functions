lcraws <- function(data, lim = 100, errs = FALSE) {
  
  start <- Sys.time()
  
  print("Plotting all objects...")
  
  #print("Objects without an AAVSO period have been removed.")
  
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
                                       DEC[i]+radii, "' AND entries >= '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    pos <- weighted.mean(ts[,2], (1 / ts[,3])) + amplitude
    
    neg <- weighted.mean(ts[,2], (1 / ts[,3])) - amplitude
    
    ts <- ts[which(ts[,2] <= pos & ts[,2] >= neg),]
    
    timdif <- max(ts[,1]) - min(ts[,1])
    
    plot(ts[,1], ts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Raw Light Curve", sep=""), xlab = "Modified Julian Date",
         ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
    magerr <- ts[,3]
    
    if (errs == TRUE){
      
      segments(ts[,1], ts[,2]-magerr, ts[,1], ts[,2]+magerr)
      
      epsilon = timdif/500
      
      segments(ts[,1]-epsilon, ts[,2]-magerr, ts[,1]+epsilon, ts[,2]-magerr)
      
      segments(ts[,1]-epsilon, ts[,2]+magerr, ts[,1]+epsilon, ts[,2]+magerr)
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}