writedata <- function(data, radii = 0.001, lim = 100) {
  
  start <- Sys.time()
  
  library(RODBC)
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  type <- as.character(as.vector(data$Type))
  
  Period <- as.numeric(as.vector(data$Period))
  
  radra <- radii / abs(cos(DEC*pi/180))
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra[i], "' AND '", RA[i]+radra[i], "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    write.table(ts, file = paste("F:/Documents/PhD/Writedata_3011/", data$Name[i], ".csv", sep=""), sep=",",
                row.names = FALSE, col.names = FALSE)
    
    #plot(ts[,1:2], pch=19)
    
  }
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
}