imgfold <- function(data, lim = 100, rang = 2, pixsize = 100) {
  
  start <- Sys.time()
  
  imgset <- list()
  
  par(mfrow=c(1, 1))
  
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
    
    imag <- matrix(0, nrow = pixsize, ncol = pixsize)
    
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
    
    ts[,1] <- ((((ts[,1] - ts[,1][which.min(ts[,2])]) / (per[i])) + 0.25) %% 1)
    
    ts <- ts[order(ts[,1]),]
    
    fold <- ts[,1] - 1
    
    folded <- c(fold, ts[,1])
    
    foldts <- c(ts[,2], ts[,2])
    
    pos <- weighted.mean(ts[,2], (1 / ts[,3])) + rang * sd(ts[,2])
    
    neg <- weighted.mean(ts[,2], (1 / ts[,3])) - rang * sd(ts[,2])
    
    magerr <- c(ts[,3], ts[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    fts <- fts[which(fts[,2] <= pos & fts[,2] >= neg),]
    
    splitx <- seq(from = -1, to = 1, length.out = pixsize+1)
    
    splity <- seq(from = min(fts[,2]), to = max(fts[,2]), length.out = pixsize+1)
    
    for (j in 1:(length(splitx)-1)){
      
      imgcol <- matrix(fts[which(fts[,1] > splitx[j] & fts[,1] <= splitx[j+1]),], ncol = 3)
      
      if (as.numeric(dim(imgcol)[1]) != 0){
        
        for (k in 1:(length(splity)-1)){
          
          imgrow <- matrix(imgcol[which(imgcol[,2] > splity[k] & imgcol[,2] <= splity[k+1]),], ncol = 3)
          
          imag[j,k] <- dim(imgrow)[1]
          
        }
        
      }
      else{
        
        imag[j,] <- 0
        
      }
      
    }
    
    imgcol <- matrix(fts[which(fts[,1] > splitx[1] & fts[,1] <= splitx[2]),], ncol = 3)
    
    imgrow <- matrix(imgcol[which(imgcol[,2] > splity[1] & imgcol[,2] <= splity[2]),], ncol = 3)
    
    imag[1,1] <- as.numeric(dim(imgrow)[1])
    
    is.na(imag) <- 0
    
    imag <- imag / max(imag)
    
    imag <- imag[, rev(seq_len(ncol(imag)))]
    
    png(paste("F:/Documents/PhD/", data$Type[i], "-", data$Name[i], ".png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  imgset
  
}