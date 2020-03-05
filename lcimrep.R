lcimrep <- function(data, lim = 100, lctrenddata, sigk = 4, errs = FALSE, quiet = FALSE) {
  
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
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[i], "'", sep=""))
    
    keep <- sigclip(info$Rcat, sig = sigk, tol = 0.000001)
    
    info <- info[keep,]
    
    ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    ts <- unique(ts)
    
    preav <- weighted.mean(ts[,2], (1 / ts[,3]))
    
    ts[,2] <- ts[,2] - preav
    
    tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    tsin <- rbind(tsin, ts)
    
    tsin <- tsin[order(tsin[,1]),]
    
    seastrend <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/365.24217) + cos(2*pi*tsin[,1]/365.24217) +
                      sin(4*pi*tsin[,1]/365.24217) + cos(4*pi*tsin[,1]/365.24217), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    seascoe <- seastrend$coefficients
    
    seascoe[is.na(seascoe)] <- 0
    
    ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    prelen <- length(ts[,2])
    
    ts <- unique(ts)
    
    #ts[which(ts[,1] < 55200),2] <- ts[which(ts[,1] < 55200),2] - 0.1
    
    plot(ts[,1:2], pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
    lines(tsin[,1], seastrend$fitted.values+preav, col = "red")
    
    seascoe <- lcsysco(as.numeric(objects$RA[1]), as.numeric(objects$DEClin[1]), lctrenddata, trcut = 50, model = "harm", quiet = TRUE)
    
    seafit <- (seascoe[1] + seascoe[2]*tsin[,1] + seascoe[3]*sin(2*pi*tsin[,1]/365.24217) + seascoe[4]*cos(2*pi*tsin[,1]/365.24217) +
                 seascoe[5]*sin(4*pi*tsin[,1]/365.24217) + seascoe[6]*cos(4*pi*tsin[,1]/365.24217))
    
    firstamp <- max(ts[,2]) - min(ts[,2])
    
    trerang <- max(seafit) - min(seafit)
    
    rangcut <- min(c(trerang, 0.2))
    
    lines(tsin[,1], seafit+preav, col = "blue")
    
    trendcor <- cor(seastrend$fitted.values, seafit)
    
    print(trendcor)
    
    if (trendcor >= 0 & trerang <= firstamp){
      
      ts[,2] <- ts[,2] - (seascoe[1] + seascoe[2]*ts[,1] + seascoe[3]*sin(2*pi*ts[,1]/365.24217) + seascoe[4]*cos(2*pi*ts[,1]/365.24217) +
                            seascoe[5]*sin(4*pi*ts[,1]/365.24217) + seascoe[6]*cos(4*pi*ts[,1]/365.24217))
      
      trendrem <- TRUE
      
    }
    else{
      
      trendrem <- FALSE
      
    }
    
    lcrpart <- try(lcrpartfit(ts, cut = rangcut), TRUE)
    
    if (class(lcrpart) != "try-error"){
      
      ts <- lcrpart$ts
      
    }
    
    #keep <- sigclip(info$Rcat, sig = sigk, tol = 0.000001)
    
    #info <- info[keep,]
    
    if (quiet == FALSE){
      
      print(paste("Object with reference id '", objects$usnoref[i],
                  "' successfully read in. Object has ", objects$entries[i], " observation(s).", sep=""))
      
    }
    
    ts <- unique(ts)
    
    ts <- ts[order(ts[,1]),]
    
    postlen <- length(ts[,2])
    
    difflen <- abs(prelen - postlen)
    
    keeplen <- abs(objects$entries[i] - prelen)
    
    if (quiet == FALSE){
      
      print(paste("There were ", keeplen, " sigma-clipped outliers found within the data.", sep = ""))
      
      print(paste("There were ", difflen, " duplicate observations found within the remaining data.", sep = ""))
      
      print(paste("These have been removed leaving ", postlen, " unique observations.", sep=""))
      
    }
    
    t <- ts[,1]
    
    timdif <- max(t) - min(t)
    
    plot(ts[,1], ts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Raw Light Curve", sep=""), xlab = "Modified Julian Date",
         ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
    
    magerr <- ts[,3]
    
    if (errs == TRUE){
      
      segments(ts[,1], ts[,2]-magerr, ts[,1], ts[,2]+magerr)
      
      epsilon = timdif/500
      
      segments(ts[,1]-epsilon, ts[,2]-magerr, ts[,1]+epsilon, ts[,2]-magerr)
      
      segments(ts[,1]-epsilon, ts[,2]+magerr, ts[,1]+epsilon, ts[,2]+magerr)
      
    }
    
    dm <- matrix(0, nrow = length(ts[,1]), ncol = length(ts[,1]))
    
    dt <- matrix(0, nrow = length(ts[,1]), ncol = length(ts[,1]))
    
    for (j in 1:length(ts[,1])){
      
      for (k in 1:length(ts[,1])){
        
        if (j != k & j < k){
          
          dm[j,k] <- ts[k,2] - ts[j,2]
          
          dt[j,k] <- ts[k,1] - ts[j,1]
          
        }
        
      }
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  mbin <- c(-5,-3,-2.5,-2,-1.5,-1,-0.5,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,0.5,1,1.5,2,2.5,3,5)
  
  tbin <- log10(c(0, (1/145),(2/145),(3/145),(4/145),(1/25),(2/25),(3/25),1.5,2.5,3.5,4.5,5.5,7,10,20,30,60,90,120,240,600,960,2000))
  
  print(length(mbin))
  
  print(length(tbin))
  
  dm <- c(dm)
  
  dt <- c(dt)
  
  imag <- matrix(0, nrow = length(mbin)-1, ncol = length(tbin)-1)
  
  imag <- matrix(0, nrow = 24, ncol = 24)
  
  del <- cbind(dt, dm)
  
  del <- del[which(del[,1] > 0),]
  
  del[,1] <- log10(del[,1])
  
  print(head(del[order(del[,2], decreasing = T),]))
  
  ab = matrix(c(-2, -2, 3.6, 2), 2, 2)
  
  datbin <- bin2(del, ab = ab, nbin = c(24, 24))
  
  imag <- datbin$nc / (nrow(del)*(nrow(del)-1))
  
  #for (j in 1:length(tbin)-1){
    
  #  for (k in 1:length(mbin)-1){
      
  #    imag[j,k] <- nrow(del[which(del[,1] > tbin[j] & del[,1] <= tbin[j+1] & del[,2] > mbin[k] & del[,2] <= mbin[k+1]),])
      
  #    print(nrow(del[which(del[,1] > tbin[j] & del[,1] <= tbin[j+1] & del[,2] > mbin[k] & del[,2] <= mbin[k+1]),]))
      
  #  }
    
  #}
  
  dm <- dm[dm != 0]
  
  dt <- dt[dt != 0]
  
  #hist(dm, breaks = 50)
  
  #hist(dt, breaks = 50)
  
  #hist(dt[dt < 365], breaks = 50)
  
  #hist(dt[dt < 10], breaks = 50)
  
  #hist(dt[dt < 1], breaks = 50)
  
  #hist(dt[dt < 5/24/60], breaks = 50)
  
  out <- list(dm = dm, dt = dt, imag = imag)
  
  imag <- imag * 255
  
  imag <- imag[, rev(seq_len(ncol(imag)))]
  
  image(imag, axes = FALSE, col = my.colors(256))
  
  out
  
}




lcimrepcom <- function(data, lctrenddata, sigk = 4){
  
  out <- matrix(0, nrow = nrow(data), ncol = 576)
  
  for (i in 1:nrow(data)){
    
    print(paste("Processing Light Curve ", i, " out of ", nrow(data), ".", sep = ""))
    
    out[i,] <- c(lcimrep(data[i,], lctrenddata = lctrenddata, sigk = sigk)$imag)
    
    gc()
    
  }
  
  out <- as.data.frame(out)
  
  out <- cbind(data, out)
  
  out
  
}