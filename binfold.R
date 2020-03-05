binfold <- function(data, lim = 100, pixx = 100, pixy = 100, sam = 0.9, bins = 100, grad = TRUE) {
  
  start <- Sys.time()
  
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
    
    imag <- matrix(0, nrow = pixx, ncol = pixy)
    
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
    
    ts[,1] <- ((((ts[,1] - ts[,1][which.min(ts[,2])]) / (per[i])) + 0.25) %% 1)
    
    ts <- ts[order(ts[,1]),]
    
    fold <- ts[,1] - 1
    
    folded <- c(fold, ts[,1])
    
    foldts <- c(ts[,2], ts[,2])
    
    pos <- weighted.mean(ts[,2], (1 / ts[,3])) + amplitude
    
    neg <- weighted.mean(ts[,2], (1 / ts[,3])) - amplitude
    
    magerr <- c(ts[,3], ts[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    fts <- fts[which(fts[,2] <= pos & fts[,2] >= neg),]
    
    if (bins < 2){
      
      print("Number of bins must be 2 or more, setting number of bins to 2.")
      
      bins <- 2
      
    }
    
    splitphase <- seq(from = -1, to = 1, length.out = bins)
    
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
    
    splitx <- seq(from = -1, to = 1, length.out = pixx+1)
    
    splity <- seq(from = neg, to = pos, length.out = pixy+1)
    
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
    
    imagm <- imag - 0.5
    
    imag <- imag[, rev(seq_len(ncol(imag)))]
    
    if (grad == FALSE){
      
      imag[which(imag >= sam)] <- 1
      
      imag[which(imag < sam)] <- 0
      
      imagm <- imag - 0.5
      
    }
    else if (grad == TRUE){
      
      imag[which(imag < sam)] <- 0
      
      imagm <- imag - 0.5
      
    }
    
    png(paste("F:/Clean/Light Curves Set/", data$Class[i], "-", data$Name[i], ".png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    imagm <- t(imagm)
    
    write.table(imagm, file = paste("F:/Clean/Matrix Set/", data$Class[i], "-", data$Name[i], ".csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    imagm <- as.vector(imagm)
    
    write.table(imagm, file = paste("F:/Clean/Vector Set/", data$Class[i], "-", data$Name[i], ".csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}