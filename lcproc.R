lcproc <- function(data, lim = 100, trials = 10000, radii = 0.1, maxp = 20000, bins = 100, norm = TRUE, len = 20, seedno = 10) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  if (bins < 2){
    
    print("Number of bins must be 2 or more, setting number of bins to 2.")
    
    bins <- 2
    
  }
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  tper <- seq(from = 0.05, to = 1000, length.out = trials)
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = len*trials)
  
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
    
    allp <- length(ts[,1])
    
    print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
    
    if (maxp < length(ts[,1])){
      
      print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    cfeats <- matrix(0, nrow = trials, ncol = len)
    
    for (z in 1:trials){
      
      if (z == 1 | z == trials | z %% floor(trials/10) == 0){
        
        print(paste("Executing trial number: '", z, "' of value: '", tper[z], "'.", sep = ""))
        
      }
      
      fts <- ts
      
      fts[,1] <- ((((fts[,1] - fts[,1][which.min(fts[,2])]) / (tper[z])) + 0.25) %% 1)
      
      fts <- fts[order(fts[,1]),]
      
      pos <- weighted.mean(ts[,2], (1 / ts[,3])) + amplitude
      
      neg <- weighted.mean(ts[,2], (1 / ts[,3])) - amplitude
      
      fts <- fts[which(fts[,2] <= pos & fts[,2] >= neg),]
      
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
      
      fts <- fts[order(fts[,1]),]
      
      ftsst <- c(0, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
      
      ftsend <- c(1, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
      
      ftsin <- rbind(ftsst, fts, ftsend)
      
      lo <- loess(ftsin[,2]~ftsin[,1], data = as.data.frame(ftsin[,1:2]), weights = 1/ftsin[,3], span = 0.5, degree = 2, model = TRUE)
      
      # Put New models here
      
      #negfold <- cbind((fts[,1]-1), fts[,2], fts[,3])
      
      #fulfts <- rbind(negfold, fts)
      
      #plot(fulfts[,1], fulfts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0", sep=""), xlab = paste("Phase at a period of ", tper[z], " days", sep=""),
      #      ylab = "Apparent Magnitude", pch=19, ylim = c(max(fulfts[,2]), min(fulfts[,2])), xlim = c(-1, 1))
      
      #predlo <- predict(lo, se = TRUE)
      
      #lines((ftsin[,1]-1), predict(lo), col='red', lwd=2)
      
      #lines(ftsin[,1], predict(lo), col='red', lwd=2)
      
      #lines(ftsin[,1], predlo$fit-1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      #lines(ftsin[,1], predlo$fit+1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      #lines((ftsin[,1]-1), predlo$fit-1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      #lines((ftsin[,1]-1), predlo$fit+1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      feats <- predict(lo, seq(0, 1, length.out = len)) - weighted.mean(fts[,2], (1 / fts[,3]))
      
      if (norm == TRUE){
        
        maxabsfeat <- which.max(abs(feats))
        
        feats <- c(feats[maxabsfeat:len], feats[1:maxabsfeat-1])
        
      }
      
      cfeats[z,] <- feats
      
    }
    
    dim(cfeats) <- c(1, len*trials)
    
    cfeats <- as.vector(cfeats)
    
    fullfeats[i,] <- cfeats
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", paste("V", 1:len, sep = ""))
  
  fullfeats
  
}