lcfeatgen <- function(ra, declin, radii = 0.01, lim = 2, seed = 20, run = 5000, alldb = FALSE, geo = FALSE, quiet = FALSE) {
  
  closeAllConnections()
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  
  start <- Sys.time()
  
  if (quiet == FALSE){
    
    print("Executing Software...")
    
  }
  
  set.seed(seed)
  
  ovsm <- 5
  
  eps <- 1e-3
  
  if (lim < 2){
    
    print("Limit too low, changing to '2' and continuing...")
    
    lim = 2
    
  }
  
  radra <- radii / abs(cos(declin*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  if (alldb == TRUE){
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries FROM objdat WHERE entries > '", lim, "'", sep=""))
    
  }
  
  else {
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                       declin+radii, "' AND entries > '", lim, "'", sep=""))
    
  }
  
  if (length(objects$entries) < 1){
    
    close(channel)
    
    if (alldb == TRUE){
      
      stop(paste("No objects detected across the whole database with greater than ", lim, " entries. Function will now exit.", sep=""))
      
    }
    
    else{
      
      stop(paste("No objects detected within search area of ", radii, " degree(s) with greater than ", lim, " entries. Function will now exit.", sep=""))
      
    }
    
  }
  
  else{
    
    minmag <- as.numeric(sqlQuery(channel, "SELECT max(r2mag) FROM catdat"))
    
    maxmag <- as.numeric(sqlQuery(channel, "SELECT min(r2mag) FROM catdat"))
    
    brcol <- rep(0, length(objects$usnoref))
    
    rmag <- rep(maxmag, length(objects$usnoref))
    
    namecol <- rep("#00FF00", length(objects$usnoref))
    
    for (i in 1:length(objects$usnoref)){
      
      colandmag <- sqlQuery(channel, paste("SELECT BRcolour, r2mag FROM catdat WHERE usnoref = '", objects$usnoref[i], "'", sep=""))
      
      brcol[i] <- colandmag$BRcolour[1]
      
      rmag[i] <- colandmag$r2mag[1]
      
      if (is.na(rmag[i])){
        
        rmag[i] <- 12
        
      }
      
      if (is.na(brcol[i])){
        
        brcol[i] <- 1000
        
      }
      
      if (brcol[i] <= -0.45){
        
        namecol[i] <- "#0000FF"
        
      }
      else if (brcol[i] > -0.45 & brcol[i] <= -0.2){
        
        namecol[i] <- "#00FFFF"
        
      }
      else if (brcol[i] > -0.2 & brcol[i] <= 0.2){
        
        namecol[i] <- "#FFFFFF"
        
      }
      else if (brcol[i] > 0.2 & brcol[i] <= 0.64){
        
        namecol[i] <- "#FFFF80"
        
      }
      else if (brcol[i] > 0.64 & brcol[i] <= 1.06){
        
        namecol[i] <- "#FFFF00"
        
      }
      else if (brcol[i] > 1.06 & brcol[i] <= 1.42){
        
        namecol[i] <- "#FF8000"
        
      }
      else if (brcol[i] > 1.42 & brcol[i] <= 100){
        
        namecol[i] <- "#FF0000"
        
      }
      else{
        
        namecol[i] <- "#00FF00"
        
      }
      
    }
    
    if (quiet == FALSE){
      
      if (alldb == TRUE){
        
        print(paste("Reading in all ", length(objects$usnoref), " object(s) with greater than ", lim, " observations.", sep=""))
        
      }
      
      else{
        
        print(paste("A total of ", length(objects$usnoref),
                    " object(s) identified within the search radius of ", radii, " degree(s) with greater than ", lim, " observations.", sep=""))
        
      }
      
      if (alldb == FALSE){
        
        par(bg = "black", col = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white")
        
        plot(objects$RA, objects$DEClin, main = paste("Coordinates RA = ", ra, " deg and DEC = ", declin, " deg, radius = ", radii, " deg.", sep=""), xlab = "Right Ascension (degrees)",
             ylab = "Declination (degrees)", pch = 19, cex = ((2.5 * ((rmag - minmag) / (maxmag - minmag)) + 0.5) * max((log10(1 / (2 * radii)) + 1), 0.1)), col = namecol, xlim = c(ra-radra, ra+radra), ylim = c(declin-radii, declin+radii))
        
        if (length(objects$usnoref) < 100){
          
          text(objects$RA, objects$DEClin, labels = objects$usnoref, cex = 0.6, pos = 4, offset = (0.5 * max((log10(1 / (2 * radii)) + 1), 0.01)), col = "#808080")
          
        }
        
        grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
             lwd = par("lwd"), equilogs = TRUE)
        
        par(bg = "white", col = "black", col.axis = "black", col.lab = "black", col.main = "black", col.sub = "black")
        
      }
      
    }
    
    featvec <- matrix(0, length(objects$entries), 16)
    
    if (geo == TRUE){
      
      if (quiet == FALSE){
        
        print("Computing geometrically closest object to the reference coordinates...")
        
      }
      
      geodist <- sqrt((declin - objects$DEClin)^2.0 + ((ra - objects$RA) / abs(cos(declin*pi/180)))^2.0)
      
      mindist <- which.min(geodist)
      
      objects <- objects[mindist,]
      
    }
    
    for (i in 1:length(objects$entries)){
      
      if (quiet == FALSE){
        
        if (geo == FALSE){
          
          print(paste("Reading next object from database...   {", i, "/", length(objects$entries), "}", sep=""))
          
        }
        else{
          
          print("Reading this object from database...")
          
        }
        
      }
      
      info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                      objects$usnoref[i], "'", sep=""))
      
      if (quiet == FALSE){
        
        print(paste("Object with reference id '", objects$usnoref[i],
                    "' successfully read in. Object has ", objects$entries[i], " observation(s).", sep=""))
        
      }
      
      ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
      
      prelen <- length(ts[,2])
      
      ts <- unique(ts)
      
      ts <- ts[order(ts[,1]),]
      
      postlen <- length(ts[,2])
      
      difflen <- abs(prelen - postlen)
      
      if (quiet == FALSE){
        
        print(paste("There were ", difflen, " duplicate observations found within the data.", sep = ""))
        
        print(paste("These have been removed leaving ", postlen, " unique observations.", sep=""))
        
      }
      
      if (run > 0 & run < as.numeric(objects$entries[i])){
        
        if (quiet == FALSE){
          
          print(paste("Number of data points greater than limit. Setting number of data points to ", run, ".", sep=""))
          
        }
        
        mints <- t(as.matrix(ts[1,]))
        
        maxts <- t(as.matrix(ts[length(ts[,1]),]))
        
        ts <- ts[2:(length(ts[,1]) - 1),]
        
        ts <- ts[order(ts[,3]),]
        
        ts <- ts[1:run,]
        
        ts <- rbind(ts, mints)
        
        ts <- rbind(ts, maxts)
        
        ts <- ts[order(ts[,1]),]
        
      }
      
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
      
      gap2 <- gap[gap < 2*sd(gap)]
      
      if (quiet == FALSE){
        
        plot(1:length(gap), gap, pch=19, type="n", main = paste("The Gaps between observations for ", objects$usnoref[i], "", sep=""),
             xlab = "Interval Number", ylab = "Interval Length (days)")
        
        lines(1:length(gap), gap, type="l")
        
        plot(1:length(gapn), gapn, pch=19, type="n", main = paste("The deviation of uneven observations from the sampling rate (every minute) for ", objects$usnoref[i], "", sep=""),
             xlab = "Interval Number", ylab = "Observational Uneven deviations (seconds)")
        
        lines(1:length(gapn), gapn, type="l")
        
        hist(ts[,1], breaks = as.integer(permax / 10), main = paste("Histogram of the distribution of time instants for ", objects$usnoref[i], "", sep=""),
             xlab = "Modified Julian Date (Days)")
        
      }
      
      dt <- mean(gap2)
        
      nyq <- 20
      
      permin <- 1 / nyq
      
      freqmin <- 1 / permax
      
      freqmax <- nyq
      
      if (quiet == FALSE){
        
        print(paste("Searching over a period range of ", permin,
                    " days to ", ceiling(permax), " days, the time-span of observations.", sep=""))
        
        print(paste("This is a frequency range of ", freqmin, " days^-1 to ", freqmax, " days^-1.", sep=""))
        
      }
      
      df <- freqmin * (1 / ovsm)
      
      if (quiet == FALSE){
        
        print(paste("The frequency step is ", df,
                    " days^-1 as defined by the oversampling rate = ", ovsm, ".", sep=""))
        
      }
      
      freqlen <- length(seq(freqmin, freqmax, by = df))
      
      if (quiet == FALSE){
        
        print(paste("The frequency grid is ", freqlen, " candidates in size.", sep=""))
        
      }
      
      
      features <- rep(0, 16)
      
      meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
      
      stddev <- sd(ts[,2])
      
      ts <- ts[which(ts[,2] < (meanmag + 3*stddev)),]
      
      ts <- ts[which(ts[,2] > (meanmag - 3*stddev)),]
      
      tsts <- ts
      
      print("Calculating Non-Periodic Features...")
      
      features[8] <- skewness(ts[,2])
      
      if (length(ts[,2]) > 5000){
        
        mags <- sample(ts[,2], 5000)
        
      }
      else{
        
        mags <- ts[,2]
        
      }
      
      features[14] <- as.numeric(shapiro.test(mags)[1])
      
      n <- length(ts[,2])
      
      numf <- length(ts[which(ts[,2] > meanmag),2])
      
      numb <- length(ts[which(ts[,2] < meanmag),2])
      
      sig2f <- sum((ts[which(ts[,2] > meanmag),2] - meanmag)^2.0) / numf
      
      sig2b <- sum((ts[which(ts[,2] < meanmag),2] - meanmag)^2.0) / numb
      
      features[12] <- sig2f / sig2b
      
      w <- 1 / ((ts[(2:n),1] - ts[(1:n-1),1]) ^ 2.0)
      
      wmean <- mean(w)
      
      S1 <- sum(w * (ts[2:n,2] - ts[1:n-1,2]) ^ 2.0)
      
      S2 <- sum(w)
      
      quan <- quantile(ts[,2])
      
      features[11] <- as.numeric(quan[4]) - as.numeric(quan[2])
      
      features[9] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((ts[,2] - meanmag) / stddev) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
      
      amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
      
      medianmag <- median(ts[,2])
      
      sortedmags <- sort(ts[,2])
      
      sigmap <- (sqrt(n / (n - 1)) * (ts[,2] - meanmag) / (1 / ts[,3]))
      
      features[10] <- (1 / sqrt(n) * sum(abs(sigmap)) / sqrt(sum(sigmap ^ 2)))
          
      sper <- spurper(ts, from = freqmin, to = freqmax, ofac = ovsm, name = objects$usnoref[i], plot = !quiet)
          
      spurs <- sper$scanned[order(sper$power, decreasing = TRUE)]
          
      spurs <- spurs[which(sper$power > (max(sper$power)/3))]
          
      spurs <- c((1 / 29.53), spurs)
            
      print("Fitting Lomb-Scargle Periodogram...")
            
      ans <- lsp(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
      pwr <- ans$power
          
      if (quiet == FALSE){
            
        plot(ans$scanned, pwr, pch=19, type="n", main = paste("'LSP' Periodogram for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
              ylab = "LSP Statistic")
            
        lines(ans$scanned, pwr, type="l")
            
      }
          
      f <- ans$scanned[order(pwr, decreasing = T)]
          
      for (zz in 1:length(spurs)){
            
        f <- f[which(!(abs((1 / f) - (1 / spurs[zz])) < eps*(1 / spurs[zz])))]
            
        f <- f[which(!(((1 / f) > (1 / spurs[zz])) & (abs(((1 / f) / (1 / spurs[zz])) - floor(((1 / f) / (1 / spurs[zz])))) < eps)))]
            
        f <- f[which(!(((1 / f) < (1 / spurs[zz])) & (abs(((1 / spurs[zz]) / (1 / f)) - floor(((1 / spurs[zz]) / (1 / f)))) < eps)))]
            
      }
          
      fo <- f[1]
          
      fine <- fineper(ts, f, "LSP")
          
      f <- fine$freq
          
      p <- 1 / f
          
      features[1] <- p
          
      effm <- 2 * fine$len / ovsm
          
      level <- -log(1 - (1 - ans$alpha)^(1/effm))
          
      exPN <- exp(-1.0 * fine$power)
          
      pv <- effm * exPN
          
      if(is.na(pv)){
            
        pv <- 1
            
      }
      else{
            
        if (pv > 0.01){
              
          pv <- 1 - (1 - exPN)^effm
              
        }
            
      }
      
    }
    
    print("Displaying #1 period detected:")
        
    print(p)
        
    if (quiet == FALSE){
          
      print("Displaying #1 period p-values:")
          
      print(pv)
          
    }
        
    magerr <- ts[,3]
        
    tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(meanmag, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
        
    tsin <- rbind(tsin, ts)
        
    tsin <- tsin[order(tsin[,1]),]
        
    if (quiet == FALSE){
          
      print("Calculating best 1st to 3rd harmonic fit for the #1 period...")
          
    }
        
    SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                  sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                  cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]) +
                  + sin(10*pi*f*tsin[,1]) + cos(10*pi*f*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
    co <- SSTlm$coefficients
        
    co[is.na(co)] <- 0
        
    ts <- tsts
        
    if (quiet == FALSE){
          
      print("Plotting Light Curve folded onto #1 period...")
          
    }
        
    pos <- weighted.mean(ts[,2], (1 / ts[,3])) + 2 * sd(ts[,2])
        
    neg <- weighted.mean(ts[,2], (1 / ts[,3])) - 2 * sd(ts[,2])
        
    folded <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (features[1])) + 1.25) %% 2) - 1
        
    folded2 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (features[1])) + 0.25) %% 1)
        
    fts <- matrix(c(folded2, ts[,2]), nrow = length(folded2))
        
    fts <- fts[order(fts[,1]),]
        
    date_diff <- fts[2:length(fts[,1]),1] - fts[1:(length(fts[,1])-1),1]
        
    mag_diff <- fts[2:length(fts[,1]),2] - fts[1:(length(fts[,1])-1),2]
        
    date_diff <- date_diff[which(mag_diff != 0)]
        
    mag_diff <- mag_diff[which(mag_diff != 0)]
        
    slope <- date_diff / mag_diff
        
    features[15] <- quantile(slope, 0.1)
        
    features[16] <- quantile(slope, 0.9)
        
    pcs <- (1 / (n * stddev)) * cumsum(fts[,2] - meanmag)
        
    features[3] <- max(pcs) - min(pcs)
        
    varfold <- (sd(fts[,2])) ^ 2.0
        
    features[2] <-  (1.0 / ((n - 1) * varfold) * sum((fts[2:n,2] - fts[1:(n-1),2]) ^ 2.0))
        
    folded3 <- folded2 - 1
        
    folded4 <- c(folded3, folded2)
        
    foldts <- c(ts[,2], ts[,2])
        
    if (quiet == FALSE){
          
      plot(folded, ts[,2], main = paste("", objects$usnoref[i], " Folded Light Curve spread from -1.0 to 1.0 for first period", sep=""), xlab = paste("Phase at the first period of ", p, " days", sep=""),
            ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
            
      magerr <- ts[,3]
            
      segments(folded, ts[,2]-magerr, folded, ts[,2]+magerr)
            
      epsilon = 0.01
            
      segments(folded-epsilon, ts[,2]-magerr, folded+epsilon, ts[,2]-magerr)
            
      segments(folded-epsilon, ts[,2]+magerr, folded+epsilon, ts[,2]+magerr)
        
      plot(folded4, foldts, main = paste("", objects$usnoref[i], " Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0 for first period", sep=""), xlab = paste("Phase at the first period of ", p, " days", sep=""),
            ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
            
      magerr <- c(ts[,3], ts[,3])
            
      segments(folded4, foldts-magerr, folded4, foldts+magerr)
            
      epsilon = 0.01
            
      segments(folded4-epsilon, foldts-magerr, folded4+epsilon, foldts-magerr)
            
      segments(folded4-epsilon, foldts+magerr, folded4+epsilon, foldts+magerr)
          
      print("Plotting Raw Light Curve with fitted model...")
          
      plot(ts[,1], ts[,2], main = paste("", objects$usnoref[i], " Raw Light Curve with fitted model", sep=""), xlab = "Modified Julian Date (Days)",
            ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
      lines(sort(tsin[,1]), SSTlm$fitted.values, col = 2)
            
      magerr <- ts[,3]
            
      segments(ts[,1], ts[,2]-magerr, ts[,1], ts[,2]+magerr)
            
      epsilon2 = 4
            
      segments(ts[,1]-epsilon2, ts[,2]-magerr, ts[,1]+epsilon2, ts[,2]-magerr)
            
      segments(ts[,1]-epsilon2, ts[,2]+magerr, ts[,1]+epsilon2, ts[,2]+magerr)
          
    }
        
    if (quiet == FALSE){
          
      print(paste("Regression Flag = ", features[2], "", sep=""))
          
    }
        
    R1 <- (co[3]^2.0 + co[4]^2.0)^0.5
        
    R2 <- (co[5]^2.0 + co[6]^2.0)^0.5
        
    R3 <- (co[7]^2.0 + co[8]^2.0)^0.5
        
    features[13] <- R1
        
    features[4] <- R2/R1
        
    features[5] <- R3/R1
        
    PH1 <- atan2(-co[3], co[4])
        
    PH2 <- atan2(-co[5], co[6]) - (2.0 * PH1)
        
    PH3 <- atan2(-co[7], co[8]) - (3.0 * PH1)
        
    features[6] <- PH2 - PH1
        
    features[7] <- PH3 - PH1
        
    featvec[i,] <- features
        
  }
      
  featvec <- as.data.frame(featvec)
      
  featvec <- cbind(objects$usnoref, featvec)
      
  colnames(featvec, do.NULL = FALSE)
  
  close(channel)
    
  print(Sys.time() - start)
    
  gc()
    
  if (geo == TRUE){
      
    featvec <- featvec[1,]
      
  }
  
  colnames(featvec) <- c("USnoref", "LS.Period", "Psi-eta", "Psi-cs",
                        "R21", "R31", "PH21", "PH31", "Skewness", "Kurtosis", "Stetson K", "Q31", "A", "H1",
                        "W", "mp10", "mp90")
    
  featvec
  
}