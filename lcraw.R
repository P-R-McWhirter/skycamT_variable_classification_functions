lcraw <- function(data, lim = 100, errs = FALSE) {
  
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
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    timdif <- max(ts[,1]) - min(ts[,1])
    
    plot(ts[,1], ts[,2], main = paste("", data$Name[i], " (", data$Type[i], ") Raw Light Curve", sep=""), xlab = "Modified Julian Date",
         ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])), col = rgb(0, 0, 0, 0.25))
    
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




lcrawfit <- function(data, lim = 100, lctrenddata, sigk = 4, errs = FALSE, quiet = FALSE) {
  
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
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}