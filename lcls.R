lcls <- function(data, lim = 100, len = 5000) {
  
  start <- Sys.time()
  
  library(RODBC)
  
  print("Analysing all objects...")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  aper <- as.numeric(as.vector(data$Period))
  
  per <- seq(from = 0.01, to = 2, length.out = len)
  
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
    
    chi2 <- rep(0, length(per))
    
    for (j in 1:length(per)){
      
      if (j %% 100 == 0){
        
        print(paste("Loop ", j, "/", len, ".", sep = ""))
        
      }
      
      Z <- cbind(1, sin(2*pi*per[j]*ts[,1]), cos(2*pi*per[j]*ts[,1]))
      
      Z <- as.matrix(Z)
      
      lambda <- 10^(-3)
      
      W <- diag(1/ts[,3])
      
      M <- t(Z)%*%W%*%Z + diag(rep(3,1))*lambda
      
      coeff <- solve(M)%*%t(Z)%*%W%*%(ts[,2])
      
      chi2[j] <- sum(((ts[,2] - (coeff[1] + coeff[2] * sin(2*pi*per[j]*ts[,1]) + coeff[3] * cos(2*pi*per[j]*ts[,1])))/ts[,3])^2)
      
    }
    
    chi2 <- max(chi2) - chi2
    
    plot(per, chi2, type = "l")
    
    fullfeats[i,] <- chi2
    
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