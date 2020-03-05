lcsplineper <- function(data, lim = 100, detrend = FALSE, spur = FALSE, radii = 0.1, maxp = 100000, pero = 1.0, bins = 100, norm = TRUE, len = 200, seedno = 10) {
  
  start <- Sys.time()
  
  set.seed(seedno)
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fullfeats <- matrix(0, nrow = length(data$Name), ncol = len)
  
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
    
    if (detrend == TRUE){
      
      fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts), weights = 1/ts[,3])
      
      coeff <- coefficients(fit)
      
      ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
      
    }
    
    allp <- length(ts[,1])
    
    print(paste("Number of datapoints for '", data$Name[i], "' is ", allp, ".", sep = ""))
    
    if (maxp < length(ts[,1])){
      
      print(paste("Reducing the number of datapoints from ", allp, " to ", maxp, ".", sep = ""))
      
      ts <- ts[sample(nrow(ts)),][1:maxp,]
      
    }
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    if (spur == TRUE){
      
      per[i] <- pero
      
    }
    
    pos <- weighted.mean(ts[,2], (1 / ts[,3])) + amplitude
    
    neg <- weighted.mean(ts[,2], (1 / ts[,3])) - amplitude
    
    ts <- ts[which(ts[,2] <= pos & ts[,2] >= neg),]
    
    goodness <- rep(0, 10000)
    
    leng <- goodness
    
    for (z in 1:10000){
      
      fts <- ts
      
      fts[,1] <- ((((fts[,1] - fts[,1][which.max(abs(fts[,2]))]) / (seq(0.1, 10, length.out = 10000)[z])) + 0.75) %% 1) - 0.5
      
      #fts <- fts[order(fts[,1]),]
      
      if (bins < 2){
        
        print("Number of bins must be 2 or more, setting number of bins to 2.")
        
        bins <- 2
        
      }
      
      splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
      
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
      
      #fts <- fts[order(fts[,1]),]
      
      fts[,2] <- (fts[,2] - min(fts[,2])) / abs(max(fts[,2]) - min(fts[,2]))
      
      leng[z] <- length(fts[,1])
      
      ftsst <- c(-0.5, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
      
      ftsend <- c(0.5, (fts[1,2] + fts[length(fts[,2]),2])/2, (fts[1,3] + fts[length(fts[,3]),3])/2)
      
      ftsin <- rbind(ftsst, fts, ftsend)
      
      lo <- loess(ftsin[,2]~ftsin[,1], data = as.data.frame(ftsin[,1:2]), weights = 1/ftsin[,3], span = 0.5, degree = 2, model = TRUE)
      
      # Put New models here
      
      #plot(fts[,1], fts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = paste("Phase at a period of ", seq(0.1, 10, length.out = 10000)[z], " days", sep=""),
      #     ylab = "Apparent Magnitude", pch=19, ylim = c(max(fts[,2]), min(fts[,2])), xlim = c(-0.5, 0.5))
      
      predlo <- predict(lo, se = TRUE)
      
      goodness[z] <- mean(abs((predlo$fit+1.96*predlo$se.fit) - (predlo$fit-1.96*predlo$se.fit)))
      
      #lines(ftsin[,1], predict(lo), col='red', lwd=2)
      
      #lines(ftsin[,1], predlo$fit-1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      #lines(ftsin[,1], predlo$fit+1.96*predlo$se.fit, lty="dashed", col="blue", lwd=1)
      
      feats <- predict(lo, seq(-0.5, 0.5, length.out = len)) - weighted.mean(fts[,2], (1 / fts[,3]))
      
      if (norm == TRUE){
        
        maxabsfeat <- which.max(abs(feats))
        
        feats <- c(feats[maxabsfeat:len], feats[1:maxabsfeat-1])
        
      }
      
      fullfeats[i,] <- feats
      
    }
    
  }
  
  print(seq(0.1, 10, length.out = 10000)[order(goodness, decreasing = T)][1:10])
  
  print("Operation Complete.")
  
  plot(seq(0.1, 10, length.out = 10000), goodness, type = "l")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullfeats <- as.data.frame(fullfeats)
  
  fullfeats <- cbind(data$Name, data$Type, data$Period, fullfeats)
  
  rownames(fullfeats) <- 1:length(data$Name)
  
  colnames(fullfeats) <- c("Name", "Type", "Period", paste("V", 1:len, sep = ""))
  
  fullfeats
  
}