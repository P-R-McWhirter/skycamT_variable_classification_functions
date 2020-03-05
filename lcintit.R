lcintit <- function(data, lim = 100, bins = 360000, fold = FALSE) {
  
  library(ash)
  
  start <- Sys.time()
  
  print("Plotting all objects...")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  fin <- matrix(0, nrow = length(data$Name), ncol = bins)
  
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
    
    ts[,2] <- ts[,2] - weighted.mean(ts[,2], (1 / ts[,3]))
    
    amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    pos <- amplitude
    
    neg <- -1.0*amplitude
    
    tsstart <- c(54884, weighted.mean(ts[,2], (1 / ts[,3])), 0)
    
    tsend <- c(56079, weighted.mean(ts[,2], (1 / ts[,3])), 0)
    
    ts <- rbind(tsstart, ts, tsend)
    
    memg <- as.numeric(tapply(ts[,2], cut(ts[,1], bins), mean))
    
    plot(ts[,1], ts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Raw Light Curve", sep=""), xlab = "Modified Julian Date",
         ylab = "Apparent Magnitude", pch=10, col = "red", ylim = c(pos, neg))
    
    points(seq(min(ts[,1]), max(ts[,1]), length.out = bins), memg, pch=19)
    
    if (!is.na(data$Period[1]) & fold == TRUE){
      
      meph <- (seq(min(ts[,1]), max(ts[,1]), length.out = bins) * (1/data$Period[i])) %% 1
      
      ph <- (ts[,1] * (1/data$Period[i])) %% 1
      
      plot(ph, ts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Folded Light Curve", sep=""), xlab = "Phase",
           ylab = "Apparent Magnitude", pch=10, col = "red", ylim = c(pos, neg))
      
      points(meph, memg, pch=19)
      
    }
    
    fin[i,] <- memg

  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fin <- as.data.frame(fin)
  
  fin
  
}




lclsptest <- function(data){
  
  test <- lcintit(data, bins = 1584000, fold = T)
  
  set <- as.data.frame((cbind(seq(54884, 56079, length.out = 1584000), t(test))))
  
  set <- set[complete.cases(set),]
  
  setlsp <- lsp(set, from = 0.001, to = 10, ofac = 5)
  
  setlsp
  
}





lctrendtrain <- function(data, bins = 1584000){
  
  lcnum <- rep(0, bins)
  
  lctrendmod <- rep(0, bins)
  
  for (i in 1:nrow(data)){
    
    print(i)
    
    lccur <- lcintit(data[i,], bins = bins)
    
    curloc <- which(!is.na(lccur))
    
    lctrendmod[curloc] <- as.numeric(as.vector(lctrendmod[curloc])) + as.numeric(as.vector(lccur[curloc]))
    
    lcnum[curloc] <- as.numeric(as.vector(lcnum[curloc])) + 1
    
  }
  
  lctrendmod <- as.numeric(as.vector(lctrendmod))
  
  lcnum <- as.numeric(as.vector(lcnum))
  
  zeroloc <- which(lcnum != 0)
  
  lctrendmod[zeroloc] <- as.numeric(as.vector(lctrendmod))[zeroloc] / as.numeric(as.vector(lcnum))[zeroloc]
  
  plot(seq(54884,56079, length.out = bins)[zeroloc], lctrendmod[zeroloc], main = "Trend Light Curve", xlab = "Modified Julian Date",
       ylab = "Apparent Magnitude", pch=19, ylim = c(max(lctrendmod), min(lctrendmod)))
  
  sp.out <- list(lctrendmod = lctrendmod, lcnum = lcnum)
  
  sp.out
  
}