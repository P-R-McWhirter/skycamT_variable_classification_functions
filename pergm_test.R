pergm_test <- function(data, lim = 100, sigk = 4, ofac = 5, nt = 150, b = 3, type = "LSP"){
  
  print("Calculating the Periodogram for all objects...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  #data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  out <- rep(0, length(per))
  
  runtime <- rep(0, length(per))
  
  entryset <- rep(0, length(per))
  
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
    
    seascoe <- lcsysco(objects$RA, objects$DEClin, lctrenddata, trcut = 50, model = "harm", quiet = TRUE)
    
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
    
    #ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    ts <- unique(ts)
    
    ts <- ts[order(ts[,1]),]
    
    postlen <- length(ts[,2])
    
    difflen <- abs(prelen - postlen)
    
    keeplen <- abs(objects$entries[1] - prelen)
      
    print(paste("There were ", keeplen, " sigma-clipped outliers found within the data.", sep = ""))
      
    print(paste("There were ", difflen, " duplicate observations found within the remaining data.", sep = ""))
      
    print(paste("These have been removed leaving ", postlen, " unique observations.", sep=""))
    
    entryset[i] <- postlen
    
    t <- ts[,1]
    
    permax <- max(t) - min(t)
    
    gap <- rep.int(0, times = length(2:length(t)))
    
    for (k in 2:length(t)){
      
      gap[k-1] <- (t[k] - t[k-1])
      
    }
    
    gapn <- ((gap * 24 * 60) %% 1) * 60
    
    for (k in 1:length(gapn)){
      
      if (gapn[k] > 30){
        
        gapn[k] <- gapn[k] - 60
        
      }
      
    }
    
    permin <- 1 / 20
    
    freqmin <- 1 / permax
    
    freqmax <- 20
    
    if (type == "LSP"){
      
      start <- as.numeric(Sys.time())
      
      result <- lsp(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ofac, plot = FALSE)
      
      runtime[i] <- as.numeric(Sys.time()) - start
      
      out[i] <- as.numeric(result$peak.at[1])
      
    }
    else if (type == "SLLK"){
      
      start <- as.numeric(Sys.time())
      
      result <- SLLK(ts, from = freqmin, to = freqmax, ofac = ofac, nt = nt, b = b, plot = FALSE)
      
      runtime[i] <- as.numeric(Sys.time()) - start
      
      out[i] <- as.numeric(result$peak.at[1])
      
    }
    else if (type == "VRP"){
      
      start <- as.numeric(Sys.time())
      
      result <- VRP2(ts, from = freqmin, to = freqmax, ofac = ofac, nt = nt, b = b, plot = FALSE)
      
      runtime[i] <- as.numeric(Sys.time()) - start
      
      out[i] <- as.numeric(result$peak.at[1])
      
    }
    else if (type == "CKP"){
      
      start <- as.numeric(Sys.time())
      
      IP <- InPot(ts)
      
      result <- CKP2(ts, IP = IP, from = freqmin, to = freqmax, ofac = ofac, nt = nt, b = b, plot = FALSE)
      
      runtime[i] <- as.numeric(Sys.time()) - start
      
      out[i] <- as.numeric(result$peak.at[1])
      
    }
    else if (type == "CE"){
      
      start <- as.numeric(Sys.time())
      
      result <- CE(ts, from = freqmin, to = freqmax, ofac = ofac, plot = T)
      
      runtime[i] <- as.numeric(Sys.time()) - start
      
      out[i] <- as.numeric(result$peak.at[1])
      
    }
    else{
      
      stop("Incorrect type.")
      
    }
    
  }
  
  out <- 1/out
  
  data <- cbind(data, entryset, out, runtime)
  
  colnames(data) <- c("Name", "AUID", "Type", "RA", "DEC", "Input.Period", "Magnitude", "Observations", "Vuong.Period", "Runtime")
  
  data
  
}







bands_test <- function(data, lim = 100, sigk = 4, ofac = 10, nt = 150, b = 3){
  
  print("Calculating the Periodogram for all objects...")
  
  print("Objects without an AAVSO period have been removed.")
  
  print(paste("Using oversampling factor of ", ofac, ".", sep = ""))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  #data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  out <- rep(0, length(per))
  
  out2 <- rep(0, length(per))
  
  runtime <- rep(0, length(per))
  
  entryset <- rep(0, length(per))
  
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
    
    seascoe <- lcsysco(objects$RA, objects$DEClin, lctrenddata, trcut = 50, model = "harm", quiet = TRUE)
    
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
    
    #ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
    
    ts <- unique(ts)
    
    ts <- ts[order(ts[,1]),]
    
    postlen <- length(ts[,2])
    
    difflen <- abs(prelen - postlen)
    
    keeplen <- abs(objects$entries[1] - prelen)
    
    print(paste("There were ", keeplen, " sigma-clipped outliers found within the data.", sep = ""))
    
    print(paste("There were ", difflen, " duplicate observations found within the remaining data.", sep = ""))
    
    print(paste("These have been removed leaving ", postlen, " unique observations.", sep=""))
    
    entryset[i] <- postlen
    
    t <- ts[,1]
    
    permax <- max(t) - min(t)
    
    gap <- rep.int(0, times = length(2:length(t)))
    
    for (k in 2:length(t)){
      
      gap[k-1] <- (t[k] - t[k-1])
      
    }
    
    gapn <- ((gap * 24 * 60) %% 1) * 60
    
    for (k in 1:length(gapn)){
      
      if (gapn[k] > 30){
        
        gapn[k] <- gapn[k] - 60
        
      }
      
    }
    
    permin <- 1 / 20
    
    freqmin <- 1 / permax
    
    freqmax <- 20
    
    start <- as.numeric(Sys.time())
    
    freqs <- bands2(ts, ofac = 5, from = freqmin, to = freqmax, nt = nt, b = b)
    
    out[i] <- 1/freqs[which.min(abs(per[i] - (1/freqs)))]
    
    out2[i] <- 1/freqs[which.min(abs(0.5*per[i] - (1/freqs)))]
    
    runtime[i] <- as.numeric(Sys.time()) - start
    
  }
  
  data <- cbind(data, entryset, out, out2, runtime)
  
  colnames(data) <- c("Name", "AUID", "Type", "RA", "DEC", "Input.Period", "Magnitude", "Observations", "Vuong.Period", "Vuong.Period2", "Runtime")
  
  data
  
}