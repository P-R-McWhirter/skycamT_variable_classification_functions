trainfold <- function(data, lim = 100, pixsize = 100, sam = 0.9, grad = TRUE) {
  
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
    
    imag <- matrix(0, nrow = pixsize, ncol = pixsize)
    
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
    
    splitx <- seq(from = -1, to = 1, length.out = pixsize+1)
    
    splity <- seq(from = neg, to = pos, length.out = pixsize+1)
    
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
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], ".png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    imagm <- t(imagm)
    
    write.table(imagm, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], ".csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    imagm <- as.vector(imagm)
    
    write.table(imagm, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], ".csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}




trainfold.old <- function(data, lim = 100, pixsize = 100, sam = 0.9, grad = TRUE) {
  
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
    
    imag3 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    imag4 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    imag5 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    imag6 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
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
    
    pos3 <- weighted.mean(ts[,2], (1 / ts[,3])) + 3.0 * sd(ts[,2])
    
    neg3 <- weighted.mean(ts[,2], (1 / ts[,3])) - 3.0 * sd(ts[,2])
    
    pos4 <- weighted.mean(ts[,2], (1 / ts[,3])) + 4.0 * sd(ts[,2])
    
    neg4 <- weighted.mean(ts[,2], (1 / ts[,3])) - 4.0 * sd(ts[,2])
    
    pos5 <- weighted.mean(ts[,2], (1 / ts[,3])) + 5.0 * sd(ts[,2])
    
    neg5 <- weighted.mean(ts[,2], (1 / ts[,3])) - 5.0 * sd(ts[,2])
    
    pos6 <- weighted.mean(ts[,2], (1 / ts[,3])) + 6.0 * sd(ts[,2])
    
    neg6 <- weighted.mean(ts[,2], (1 / ts[,3])) - 6.0 * sd(ts[,2])
    
    magerr <- c(ts[,3], ts[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    fts3 <- fts[which(fts[,2] <= pos3 & fts[,2] >= neg3),]
    
    fts4 <- fts[which(fts[,2] <= pos4 & fts[,2] >= neg4),]
    
    fts5 <- fts[which(fts[,2] <= pos5 & fts[,2] >= neg5),]
    
    fts6 <- fts[which(fts[,2] <= pos6 & fts[,2] >= neg6),]
    
    splitx <- seq(from = -1, to = 1, length.out = pixsize+1)
    
    splity3 <- seq(from = neg3, to = pos3, length.out = pixsize+1)
    
    splity4 <- seq(from = neg4, to = pos4, length.out = pixsize+1)
    
    splity5 <- seq(from = neg5, to = pos5, length.out = pixsize+1)
    
    splity6 <- seq(from = neg6, to = pos6, length.out = pixsize+1)
    
    for (j in 1:(length(splitx)-1)){
      
      imgcol3 <- matrix(fts3[which(fts3[,1] > splitx[j] & fts3[,1] <= splitx[j+1]),], ncol = 3)
      
      imgcol4 <- matrix(fts4[which(fts4[,1] > splitx[j] & fts4[,1] <= splitx[j+1]),], ncol = 3)
      
      imgcol5 <- matrix(fts5[which(fts5[,1] > splitx[j] & fts5[,1] <= splitx[j+1]),], ncol = 3)
      
      imgcol6 <- matrix(fts6[which(fts6[,1] > splitx[j] & fts6[,1] <= splitx[j+1]),], ncol = 3)
      
      if (as.numeric(dim(imgcol3)[1]) != 0){
        
        for (k in 1:(length(splity3)-1)){
          
          imgrow3 <- matrix(imgcol3[which(imgcol3[,2] > splity3[k] & imgcol3[,2] <= splity3[k+1]),], ncol = 3)
          
          imag3[j,k] <- dim(imgrow3)[1]
          
        }
        
      }
      else{
        
        imag3[j,] <- 0
        
      }
      
      if (as.numeric(dim(imgcol4)[1]) != 0){
        
        for (k in 1:(length(splity4)-1)){
          
          imgrow4 <- matrix(imgcol4[which(imgcol4[,2] > splity4[k] & imgcol4[,2] <= splity4[k+1]),], ncol = 3)
          
          imag4[j,k] <- dim(imgrow4)[1]
          
        }
        
      }
      else{
        
        imag4[j,] <- 0
        
      }
      
      if (as.numeric(dim(imgcol5)[1]) != 0){
        
        for (k in 1:(length(splity5)-1)){
          
          imgrow5 <- matrix(imgcol5[which(imgcol5[,2] > splity5[k] & imgcol5[,2] <= splity5[k+1]),], ncol = 3)
          
          imag5[j,k] <- dim(imgrow5)[1]
          
        }
        
      }
      else{
        
        imag5[j,] <- 0
        
      }
      
      if (as.numeric(dim(imgcol6)[1]) != 0){
        
        for (k in 1:(length(splity6)-1)){
          
          imgrow6 <- matrix(imgcol6[which(imgcol6[,2] > splity6[k] & imgcol6[,2] <= splity6[k+1]),], ncol = 3)
          
          imag6[j,k] <- dim(imgrow6)[1]
          
        }
        
      }
      else{
        
        imag6[j,] <- 0
        
      }
      
    }    
    
    imgcol3 <- matrix(fts3[which(fts3[,1] > splitx[1] & fts3[,1] <= splitx[2]),], ncol = 3)
    
    imgrow3 <- matrix(imgcol3[which(imgcol3[,2] > splity3[1] & imgcol3[,2] <= splity3[2]),], ncol = 3)
    
    imgcol4 <- matrix(fts4[which(fts4[,1] > splitx[1] & fts4[,1] <= splitx[2]),], ncol = 3)
    
    imgrow4 <- matrix(imgcol4[which(imgcol4[,2] > splity4[1] & imgcol4[,2] <= splity4[2]),], ncol = 3)
    
    imgcol5 <- matrix(fts5[which(fts5[,1] > splitx[1] & fts5[,1] <= splitx[2]),], ncol = 3)
    
    imgrow5 <- matrix(imgcol5[which(imgcol5[,2] > splity5[1] & imgcol5[,2] <= splity5[2]),], ncol = 3)
    
    imgcol6 <- matrix(fts6[which(fts6[,1] > splitx[1] & fts6[,1] <= splitx[2]),], ncol = 3)
    
    imgrow6 <- matrix(imgcol6[which(imgcol6[,2] > splity6[1] & imgcol6[,2] <= splity6[2]),], ncol = 3)
    
    imag3[1,1] <- as.numeric(dim(imgrow3)[1])
    
    imag4[1,1] <- as.numeric(dim(imgrow4)[1])
    
    imag5[1,1] <- as.numeric(dim(imgrow5)[1])
    
    imag6[1,1] <- as.numeric(dim(imgrow6)[1])
    
    is.na(imag3) <- 0
    
    is.na(imag4) <- 0
    
    is.na(imag5) <- 0
    
    is.na(imag6) <- 0
    
    imag3 <- imag3 / max(imag3)
    
    imag3m <- imag3 - 0.5
    
    imag3 <- imag3[, rev(seq_len(ncol(imag3)))]
    
    imag4 <- imag4 / max(imag4)
    
    imag4m <- imag4 - 0.5
    
    imag4 <- imag4[, rev(seq_len(ncol(imag4)))]
    
    imag5 <- imag5 / max(imag5)
    
    imag5m <- imag5 - 0.5
    
    imag5 <- imag5[, rev(seq_len(ncol(imag5)))]
    
    imag6 <- imag6 / max(imag6)
    
    imag6m <- imag6 - 0.5
    
    imag6 <- imag6[, rev(seq_len(ncol(imag6)))]
    
    if (grad == FALSE){
      
      imag3[which(imag3 >= sam)] <- 1
      
      imag3[which(imag3 < sam)] <- 0
      
      imag4[which(imag4 >= sam)] <- 1
      
      imag4[which(imag4 < sam)] <- 0
      
      imag5[which(imag5 >= sam)] <- 1
      
      imag5[which(imag5 < sam)] <- 0
      
      imag6[which(imag6 >= sam)] <- 1
      
      imag6[which(imag6 < sam)] <- 0
      
      imag3m <- imag3 - 0.5
      
      imag4m <- imag4 - 0.5
      
      imag5m <- imag5 - 0.5
      
      imag6m <- imag6 - 0.5
      
    }
    else if (grad == TRUE){
      
      imag3[which(imag3 < sam)] <- 0
      
      imag4[which(imag4 < sam)] <- 0
      
      imag5[which(imag5 < sam)] <- 0
      
      imag6[which(imag6 < sam)] <- 0
      
      imag3m <- imag3 - 0.5
      
      imag4m <- imag4 - 0.5
      
      imag5m <- imag5 - 0.5
      
      imag6m <- imag6 - 0.5
      
    }
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-3.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag3, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-4.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag4, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-5.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag5, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-6.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag6, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    imag3m <- t(imag3m)
    
    imag4m <- t(imag4m)
    
    imag5m <- t(imag5m)
    
    imag6m <- t(imag6m)
    
    write.table(imag3m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-3.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag4m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-4.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag5m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-5.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag6m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-6.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    imag3m <- as.vector(imag3m)
    
    imag4m <- as.vector(imag4m)
    
    imag5m <- as.vector(imag5m)
    
    imag6m <- as.vector(imag6m)
    
    write.table(imag3m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-3.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag4m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-4.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag5m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-5.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag6m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-6.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}





trainraw <- function(data, lim = 100, xpixsize = 100, ypixsize = 100, sam = 0.9, grad = TRUE) {
  
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
    
    imag <- matrix(0, nrow = xpixsize, ncol = ypixsize)
    
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
    
    splitx <- seq(from = 54891, to = 55987, length.out = xpixsize+1)
    
    splity <- seq(from = -1, to = 12, length.out = ypixsize+1)
    
    for (j in 1:(length(splitx)-1)){
      
      imgcol <- matrix(ts[which(ts[,1] > splitx[j] & ts[,1] <= splitx[j+1]),], ncol = 3)
      
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
    
    imgcol <- matrix(ts[which(ts[,1] > splitx[1] & ts[,1] <= splitx[2]),], ncol = 3)
    
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
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], ".png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    imagm <- t(imagm)
    
    write.table(imagm, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], ".csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    imagm <- as.vector(imagm)
    
    write.table(imagm, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], ".csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}




trainfold.old <- function(data, lim = 100, pixsize = 100, sam = 0.9, grad = TRUE) {
  
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
    
    imag3 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    imag4 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    imag5 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
    imag6 <- matrix(0, nrow = pixsize, ncol = pixsize)
    
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
    
    pos3 <- weighted.mean(ts[,2], (1 / ts[,3])) + 3.0 * sd(ts[,2])
    
    neg3 <- weighted.mean(ts[,2], (1 / ts[,3])) - 3.0 * sd(ts[,2])
    
    pos4 <- weighted.mean(ts[,2], (1 / ts[,3])) + 4.0 * sd(ts[,2])
    
    neg4 <- weighted.mean(ts[,2], (1 / ts[,3])) - 4.0 * sd(ts[,2])
    
    pos5 <- weighted.mean(ts[,2], (1 / ts[,3])) + 5.0 * sd(ts[,2])
    
    neg5 <- weighted.mean(ts[,2], (1 / ts[,3])) - 5.0 * sd(ts[,2])
    
    pos6 <- weighted.mean(ts[,2], (1 / ts[,3])) + 6.0 * sd(ts[,2])
    
    neg6 <- weighted.mean(ts[,2], (1 / ts[,3])) - 6.0 * sd(ts[,2])
    
    magerr <- c(ts[,3], ts[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    fts3 <- fts[which(fts[,2] <= pos3 & fts[,2] >= neg3),]
    
    fts4 <- fts[which(fts[,2] <= pos4 & fts[,2] >= neg4),]
    
    fts5 <- fts[which(fts[,2] <= pos5 & fts[,2] >= neg5),]
    
    fts6 <- fts[which(fts[,2] <= pos6 & fts[,2] >= neg6),]
    
    splitx <- seq(from = -1, to = 1, length.out = pixsize+1)
    
    splity3 <- seq(from = neg3, to = pos3, length.out = pixsize+1)
    
    splity4 <- seq(from = neg4, to = pos4, length.out = pixsize+1)
    
    splity5 <- seq(from = neg5, to = pos5, length.out = pixsize+1)
    
    splity6 <- seq(from = neg6, to = pos6, length.out = pixsize+1)
    
    for (j in 1:(length(splitx)-1)){
      
      imgcol3 <- matrix(fts3[which(fts3[,1] > splitx[j] & fts3[,1] <= splitx[j+1]),], ncol = 3)
      
      imgcol4 <- matrix(fts4[which(fts4[,1] > splitx[j] & fts4[,1] <= splitx[j+1]),], ncol = 3)
      
      imgcol5 <- matrix(fts5[which(fts5[,1] > splitx[j] & fts5[,1] <= splitx[j+1]),], ncol = 3)
      
      imgcol6 <- matrix(fts6[which(fts6[,1] > splitx[j] & fts6[,1] <= splitx[j+1]),], ncol = 3)
      
      if (as.numeric(dim(imgcol3)[1]) != 0){
        
        for (k in 1:(length(splity3)-1)){
          
          imgrow3 <- matrix(imgcol3[which(imgcol3[,2] > splity3[k] & imgcol3[,2] <= splity3[k+1]),], ncol = 3)
          
          imag3[j,k] <- dim(imgrow3)[1]
          
        }
        
      }
      else{
        
        imag3[j,] <- 0
        
      }
      
      if (as.numeric(dim(imgcol4)[1]) != 0){
        
        for (k in 1:(length(splity4)-1)){
          
          imgrow4 <- matrix(imgcol4[which(imgcol4[,2] > splity4[k] & imgcol4[,2] <= splity4[k+1]),], ncol = 3)
          
          imag4[j,k] <- dim(imgrow4)[1]
          
        }
        
      }
      else{
        
        imag4[j,] <- 0
        
      }
      
      if (as.numeric(dim(imgcol5)[1]) != 0){
        
        for (k in 1:(length(splity5)-1)){
          
          imgrow5 <- matrix(imgcol5[which(imgcol5[,2] > splity5[k] & imgcol5[,2] <= splity5[k+1]),], ncol = 3)
          
          imag5[j,k] <- dim(imgrow5)[1]
          
        }
        
      }
      else{
        
        imag5[j,] <- 0
        
      }
      
      if (as.numeric(dim(imgcol6)[1]) != 0){
        
        for (k in 1:(length(splity6)-1)){
          
          imgrow6 <- matrix(imgcol6[which(imgcol6[,2] > splity6[k] & imgcol6[,2] <= splity6[k+1]),], ncol = 3)
          
          imag6[j,k] <- dim(imgrow6)[1]
          
        }
        
      }
      else{
        
        imag6[j,] <- 0
        
      }
      
    }    
    
    imgcol3 <- matrix(fts3[which(fts3[,1] > splitx[1] & fts3[,1] <= splitx[2]),], ncol = 3)
    
    imgrow3 <- matrix(imgcol3[which(imgcol3[,2] > splity3[1] & imgcol3[,2] <= splity3[2]),], ncol = 3)
    
    imgcol4 <- matrix(fts4[which(fts4[,1] > splitx[1] & fts4[,1] <= splitx[2]),], ncol = 3)
    
    imgrow4 <- matrix(imgcol4[which(imgcol4[,2] > splity4[1] & imgcol4[,2] <= splity4[2]),], ncol = 3)
    
    imgcol5 <- matrix(fts5[which(fts5[,1] > splitx[1] & fts5[,1] <= splitx[2]),], ncol = 3)
    
    imgrow5 <- matrix(imgcol5[which(imgcol5[,2] > splity5[1] & imgcol5[,2] <= splity5[2]),], ncol = 3)
    
    imgcol6 <- matrix(fts6[which(fts6[,1] > splitx[1] & fts6[,1] <= splitx[2]),], ncol = 3)
    
    imgrow6 <- matrix(imgcol6[which(imgcol6[,2] > splity6[1] & imgcol6[,2] <= splity6[2]),], ncol = 3)
    
    imag3[1,1] <- as.numeric(dim(imgrow3)[1])
    
    imag4[1,1] <- as.numeric(dim(imgrow4)[1])
    
    imag5[1,1] <- as.numeric(dim(imgrow5)[1])
    
    imag6[1,1] <- as.numeric(dim(imgrow6)[1])
    
    is.na(imag3) <- 0
    
    is.na(imag4) <- 0
    
    is.na(imag5) <- 0
    
    is.na(imag6) <- 0
    
    imag3 <- imag3 / max(imag3)
    
    imag3m <- imag3 - 0.5
    
    imag3 <- imag3[, rev(seq_len(ncol(imag3)))]
    
    imag4 <- imag4 / max(imag4)
    
    imag4m <- imag4 - 0.5
    
    imag4 <- imag4[, rev(seq_len(ncol(imag4)))]
    
    imag5 <- imag5 / max(imag5)
    
    imag5m <- imag5 - 0.5
    
    imag5 <- imag5[, rev(seq_len(ncol(imag5)))]
    
    imag6 <- imag6 / max(imag6)
    
    imag6m <- imag6 - 0.5
    
    imag6 <- imag6[, rev(seq_len(ncol(imag6)))]
    
    if (grad == FALSE){
      
      imag3[which(imag3 >= sam)] <- 1
      
      imag3[which(imag3 < sam)] <- 0
      
      imag4[which(imag4 >= sam)] <- 1
      
      imag4[which(imag4 < sam)] <- 0
      
      imag5[which(imag5 >= sam)] <- 1
      
      imag5[which(imag5 < sam)] <- 0
      
      imag6[which(imag6 >= sam)] <- 1
      
      imag6[which(imag6 < sam)] <- 0
      
      imag3m <- imag3 - 0.5
      
      imag4m <- imag4 - 0.5
      
      imag5m <- imag5 - 0.5
      
      imag6m <- imag6 - 0.5
      
    }
    else if (grad == TRUE){
      
      imag3[which(imag3 < sam)] <- 0
      
      imag4[which(imag4 < sam)] <- 0
      
      imag5[which(imag5 < sam)] <- 0
      
      imag6[which(imag6 < sam)] <- 0
      
      imag3m <- imag3 - 0.5
      
      imag4m <- imag4 - 0.5
      
      imag5m <- imag5 - 0.5
      
      imag6m <- imag6 - 0.5
      
    }
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-3.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag3, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-4.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag4, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-5.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag5, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    png(paste("F:/Documents/PhD/Light Curves Set/", data$Type[i], "-", data$Name[i], "-6.png", sep=""))
    
    par(mar = rep(0, 4))
    
    image(imag6, axes = FALSE, col = grey(seq(0, 1, length = 256)))
    
    dev.off()
    
    imag3m <- t(imag3m)
    
    imag4m <- t(imag4m)
    
    imag5m <- t(imag5m)
    
    imag6m <- t(imag6m)
    
    write.table(imag3m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-3.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag4m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-4.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag5m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-5.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag6m, file = paste("F:/Documents/PhD/Matrix Set/", data$Type[i], "-", data$Name[i], "-6.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    imag3m <- as.vector(imag3m)
    
    imag4m <- as.vector(imag4m)
    
    imag5m <- as.vector(imag5m)
    
    imag6m <- as.vector(imag6m)
    
    write.table(imag3m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-3.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag4m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-4.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag5m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-5.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(imag6m, file = paste("F:/Documents/PhD/Vector Set/", data$Type[i], "-", data$Name[i], "-6.csv", sep=""), sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}