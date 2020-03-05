lcan <- function(ra, declin, radii = 0.01, lim = 2, ovsm = 2, run = 5000, type = "LSP", rang = 2, ploterr = FALSE, alldb = FALSE, usenyq = TRUE, geo = FALSE, nonper = FALSE, quiet = FALSE) {
  
  closeAllConnections()
  
  #Sys.setenv(SPARK_HOME = "C:/Apache/spark-1.6.1")
  
  #.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  #library(SparkR)
  
  #sparkR.stop()
  
  start <- Sys.time()
  
  if (quiet == FALSE){
    
    print("Executing Software...")
    
  }
  
  #sc <- sparkR.init(master = "spark://150.204.48.130:7077")
  
  #sqlContext <- sparkRSQL.init(sc)
  
  set.seed(20)
  eps <- 1e-3
  
  nt = 2500
  b = 3
  
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
    
    featvec <- matrix(0, length(objects$entries), 62)
    
    if (geo == TRUE){
      
      if (quiet == FALSE){
        
        print("Computing geometrically closest object to the reference coordinates...")
        
      }
      
      geodist <- sqrt((declin - objects$DEClin)^2.0 + (ra - objects$RA)^2.0)
      
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
      
      #ts <- macho_6_6545_5[,1:3]
      
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
      
      if (usenyq == TRUE){
        
        nyq <- 0.5 * (1 / dt)
        
      }
      
      else{
        
        nyq <- 20
        
      }
      
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
      
      if (type == "LSPL"){
        
        if (permin > 10){
          
          freqlen <- length(seq(permin, permax, by = 100*df))
          
        }
        else if (permax < 10){
          
          freqlen <- length(seq(freqmin/10, freqmax/10, by = df))
          
        }
        else{
          
          freqlen <- length(seq(0.1, freqmax/10, by = df)) + length(seq(1, permax, by = 100*df))
          
        }
        
      }
      else{
        
        freqlen <- length(seq(freqmin, freqmax, by = df))
        
      }
      
      if (quiet == FALSE){
        
        print(paste("The frequency grid is ", freqlen, " candidates in size.", sep=""))
        
      }
      
      #tn <- seq(from = min(info$MJD), to = max(info$MJD), length.out = 3110)
      
      #t <- info$MJD
      
      #y <- (0.2*sin(2*pi*t/123.3) + 0.3*cos(2*pi*t/123.3) +
      #        0.09*sin(4*pi*t/123.3) + 0.01*cos(4*pi*t/123.3) + 0.03*sin(6*pi*t/123.3) +
      #        0.06*cos(6*pi*t/123.3) + 0.05*sin(8*pi*t/123.3) + 0.02*cos(8*pi*t/123.3) + 0.4*sin(2*pi*t/17.2) + 0.1*cos(2*pi*t/17.2) +
      #        0.085*sin(4*pi*t/17.2) + 0.096*cos(4*pi*t/17.2) + 0.03*sin(6*pi*t/17.2) +
      #        0.02*cos(6*pi*t/17.2) + 0.05*sin(8*pi*t/17.2) + 0.028*cos(8*pi*t/17.2) +
      #        0.35*sin(2*pi*t/42.8) + 0.25*cos(2*pi*t/42.8) + 0.04*sin(4*pi*t/42.8) + 0.032*cos(4*pi*t/42.8) + 0.024*sin(6*pi*t/42.8) +
      #        0.092*cos(6*pi*t/42.8) + 0.05*sin(8*pi*t/42.8) + 0.07*cos(8*pi*t/42.8)) +  rnorm(length(t), mean = 0, sd = 0.5) + 3
      
      #yn <- (0.2*sin(2*pi*tn/123.3) + 0.3*cos(2*pi*tn/123.3) +
      #         0.09*sin(4*pi*tn/123.3) + 0.01*cos(4*pi*tn/123.3) + 0.03*sin(6*pi*tn/123.3) +
      #         0.06*cos(6*pi*tn/123.3) + 0.05*sin(8*pi*tn/123.3) + 0.02*cos(8*pi*tn/123.3) + 0.4*sin(2*pi*tn/17.2) + 0.1*cos(2*pi*tn/17.2) +
      #         0.085*sin(4*pi*tn/17.2) + 0.096*cos(4*pi*tn/17.2) + 0.03*sin(6*pi*tn/17.2) +
      #         0.02*cos(6*pi*tn/17.2) + 0.05*sin(8*pi*tn/17.2) + 0.028*cos(8*pi*tn/17.2) +
      #         0.35*sin(2*pi*tn/42.8) + 0.25*cos(2*pi*tn/42.8) + 0.04*sin(4*pi*tn/42.8) + 0.032*cos(4*pi*tn/42.8) + 0.024*sin(6*pi*tn/42.8) +
      #         0.092*cos(6*pi*tn/42.8) + 0.05*sin(8*pi*tn/42.8) + 0.07*cos(8*pi*tn/42.8)) + 3
      
      #ts <- matrix(c(t, y), nrow = length(t))
      
      #print("Removing Linear Trends...")
      
      features <- rep(0, 62)
      
      #fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts), weights = 1/info$Rcaterr)
      
      #coeff <- coefficients(fit)
      
      features[1] <- objects$entries[i]
      
      tsts <- ts
      
      #ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
      
      print("Calculating Non-Periodic Features...")
      
      IP <- InPot(ts)
      
      stddev <- sd(ts[,2])
      
      features[41] <- stddev
      
      features[39] <- skewness(ts[,2])
      
      n <- length(ts[,2])
      
      meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
      
      w <- 1 / ((ts[(2:n),1] - ts[(1:n-1),1]) ^ 2.0)
      
      wmean <- mean(w)
      
      S1 <- sum(w * (ts[2:n,2] - ts[1:n-1,2]) ^ 2.0)
      
      S2 <- sum(w)
      
      features[34] <- meanmag
      
      features[35] <- stddev / meanmag
      
      features[36] <- (wmean * ((ts[(n-1),1] - ts[1,1]) ^ 2.0) * S1 / ((stddev ^ 2.0) * S2 * (n ^ 2.0)))
      
      quan <- quantile(ts[,2])
      
      features[38] <- as.numeric(quan[4]) - as.numeric(quan[2])
      
      features[40] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((ts[,2] - meanmag) / stddev) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
      
      cs <- (1 / (n * stddev)) * cumsum(ts[,2] - meanmag)
      
      features[43] <- max(cs) - min(cs)
      
      features[45] <- (length(ts[,2][ts[,2] > (meanmag + stddev)]) + length(ts[,2][ts[,2] < (meanmag - stddev)])) / length(ts[,2])
      
      sigmap <- (sqrt(n / (n - 1)) * (ts[,2] - meanmag) / (1 / ts[,3]))
      
      features[46] <- (1 / sqrt(n) * sum(abs(sigmap)) / sqrt(sum(sigmap ^ 2)))
      
      features[48] <- 0
      
      for (k in 2:n){
        
        slope <- (ts[k,2] - ts[(k-1),2]) / (ts[k,1] - ts[(k-1),1])
        
        if ((abs(slope) > abs(features[48])) & !is.na(slope)){
          
          features[48] <- slope
          
        }
        
      }
      
      features[48] <- abs(features[48])
      
      #amplitude <- (max(ts[,2]) - min(ts[,2])) / 2
      
      amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
      
      features[49] <- amplitude
      
      features[50] <- median(abs(ts[,2] - median(ts[,2])))
      
      medianmag <- median(ts[,2])
      
      features[51] <- 1 - ((length(ts[,2][ts[,2] > (medianmag + (amplitude/10))]) + length(ts[,2][ts[,2] < (medianmag - (amplitude/10))])) / length(ts[,2]))
      
      if (n >= 30){
        
        datalast <- ts[,2][length(ts[,2]) - (29:0)]
        
        features[52] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / 30)
        
      }
      else
      {
        
        datalast <- ts[,2][length(ts[,2]) - ((n-1):0)]
        
        features[52] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / n)
        
      }
      
      sortedmags <- sort(ts[,2])
      
      F_60_index <- ceiling(0.60 * n)
      
      F_40_index <- ceiling(0.40 * n)
      
      F_5_index <- ceiling(0.05 * n)
      
      F_95_index <- ceiling(0.95 * n)
      
      F_40_60 <- sortedmags[F_60_index] - sortedmags[F_40_index]
      
      F_5_95 <- sortedmags[F_95_index] - sortedmags[F_5_index]
      
      features[53] <- F_40_60 / F_5_95
      
      F_325_index <- ceiling(0.325 * n)
      
      F_675_index <- ceiling(0.675 * n)
      
      F_325_675 <- sortedmags[F_675_index] - sortedmags[F_325_index]
      
      features[54] <- F_325_675 / F_5_95
      
      F_25_index <- ceiling(0.25 * n)
      
      F_75_index <- ceiling(0.75 * n)
      
      F_25_75 <- sortedmags[F_75_index] - sortedmags[F_25_index]
      
      features[55] <- F_25_75 / F_5_95
      
      F_175_index <- ceiling(0.175 * n)
      
      F_825_index <- ceiling(0.825 * n)
      
      F_175_825 <- sortedmags[F_825_index] - sortedmags[F_175_index]
      
      features[56] <- F_175_825 / F_5_95
      
      F_10_index <- ceiling(0.10 * n)
      
      F_90_index <- ceiling(0.90 * n)
      
      F_10_90 <- sortedmags[F_90_index] - sortedmags[F_10_index]
      
      features[57] <- F_10_90 / F_5_95
      
      distmed <- abs(ts[,2] - medianmag)
      
      maxdist <- max(distmed)
      
      features[58] <- maxdist / medianmag
      
      features[59] <- F_5_95 / medianmag
      
      adnt <- ad.test(ts[,2])
      
      features[60] <- (1 / (1.0 + exp(-10.0 * (as.numeric(adnt$statistic) - 0.3))))
      
      lagval <- 100
      
      autocorlen <- 0
      
      while (autocorlen == 0){
        
        AC <- acf(ts[,2], lag.max = lagval, type = "correlation", plot = FALSE)
        
        racf <- as.numeric(AC$acf)
        
        autocorlen <- match(max(racf[racf < exp(-1)]), racf)
        
        features[61] <- autocorlen
        
        if (is.na(autocorlen)){
          
          autocorlen <- 0
          
        }
        
        lagval <- lagval + 100
        
      }
      
      sac <- sacf(ts[,1:2])
      
      features[62] <- sac$Length
      
      sacl <- sac$SAC
      
      nac <- length(sacl)
      
      sigmapac <- sqrt(nac / (nac - 1)) * ((sacl - mean(sacl)) / sd(sacl))
      
      features[47] <- (1 / sqrt(nac) * sum(abs(sigmapac)) / sqrt(sum(sigmapac ^ 2)))
      
      features[42] <- -1.0 * log(IP)
      
      if (nonper == TRUE){
        
        featvec[i,] <- features
        
      }
      
      if (nonper == FALSE){
        
        sper <- spurper(ts, from = freqmin, to = freqmax, ofac = ovsm, name = objects$usnoref[i], plot = !quiet)
        
        spurs <- sper$scanned[order(sper$power, decreasing = TRUE)]
        
        spurs <- spurs[which(sper$power > (max(sper$power)/3))]
        
        spurs <- c((1 / 29.53), spurs)
        
        if (type == "LSP"){
          
          print("Fitting Lomb-Scargle Periodogram...   [1/3]")
          
          ans <- lsp(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "BGLS"){
          
          print("Fitting Bayesian Generalised Lomb-Scargle Periodogram...   [1/3]")
          
          ans <- lsp(ts[,1], ts[,2], rep(1, length(ts[,1])), phigh = 1/freqmin, plow = 1/freqmax, ofac = 1)
          
        }
        else if (type == "LSPL"){
          
          print("Fitting Hybrid-Log Lomb-Scargle Periodogram...   [1/3]")
          
          ans <- lspl(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "GLS"){
          
          print("Fitting Generalised Lomb-Scargle Periodogram...   [1/3]")
          
          ans <- GLS(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
          
        }
        else if (type == "VRP"){
          
          print("Fitting Variance Ratio Periodogram...   [1/3]")
          
          ans <- VRP(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
          
        }
        else if (type == "SLLK"){
          
          print("Fitting String Length Lafler-Kinman Periodogram...   [1/3]")
          
          ans <- SLLK(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
          
        }
        else if (type == "BKR"){
          
          print("Fitting Blum-Kiefer-Rosenblatt Periodogram...   [1/3]")
          
          ans <- BKR(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
          
        }
        else if (type == "CKP"){
          
          print("Fitting Correntropy Kernelised Periodogram...   [1/3]")
          
          ans <- CKP(ts, IP, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
          
        }
        else if (type == "IFA"){
          
          print("Fitting Image Feature Analysis Periodogram...   [1/3]")
          
          ans <- IFA(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, pixsize = 28)
          
        }
        else if (type == "CNN"){
          
          print("Fitting using a Convolutional Neural Network...   [1/3]")
          
          ans <- CNN(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, pixsize = 28, amp = amplitude, sam = 0.2, grad = FALSE)
          
        }
        else if (type == "FSF"){
          
          print("Fitting using folded shape features...   [1/3]")
          
          ans <- FSF(ts, filt = spurs, melc = meanmag, amplc = amplitude, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b)
          
        }
        else if (type == "GEN"){
          
          print("Fitting using a genetic algorithm...   [1/3]")
          
          ans <- GEN(ts, filt = spurs, melc = meanmag, amplc = amplitude, ofac = ovsm, from = freqmin, to = freqmax, pop = 1000, pairups = 200)
          
        }
        else
        {
          
          close(channel)
          
          stop(paste("'", type, "' is an invalid Periodogram. Program will now exit.", sep=""))
          
        }
        
        pwr <- ans$power
        
        if (quiet == FALSE){
          
          plot(ans$scanned, pwr, pch=19, type="n", main = paste("'", type, "' Periodogram for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
               ylab = paste("'", type, "' Statistic", sep=""))
          
          lines(ans$scanned, pwr, type="l")
          
        }
        
        # Possibly add low period penalty term here log(f/fn) where fn = 0.5<1/delT>
        
        f <- ans$scanned[order(pwr, decreasing = T)]
        
        for (zz in 1:length(spurs)){
          
          f <- f[which(!(abs((1 / f) - (1 / spurs[zz])) < eps*(1 / spurs[zz])))]
          
          f <- f[which(!(((1 / f) > (1 / spurs[zz])) & (abs(((1 / f) / (1 / spurs[zz])) - floor(((1 / f) / (1 / spurs[zz])))) < eps)))]
          
          f <- f[which(!(((1 / f) < (1 / spurs[zz])) & (abs(((1 / spurs[zz]) / (1 / f)) - floor(((1 / spurs[zz]) / (1 / f)))) < eps)))]
          
        }
        
        fo <- f[1]
        
        fine <- fineper(ts, f, type)
        
        f <- fine$freq
        
        p <- 1 / f
        
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
          
          print("Calculating best 1st to 4th harmonic fit for the #1 period...")
          
        }
        
        SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                      sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                      cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
        coeffper1 <- SSTlm$coefficients
        
        coeffper1[is.na(coeffper1)] <- 0
        
        prevar <- var(ts[,2])
        
        if (quiet == FALSE){
          
          print("Subtracting harmonic fit of period #1 from light curve (prewhitening)...")
          
        }
        
        ts[,2] <- ts[,2] - (coeffper1[2]*ts[,1] + coeffper1[3]*sin(2*pi*f*ts[,1]) + coeffper1[4]*cos(2*pi*f*ts[,1]) +
                              coeffper1[5]*sin(4*pi*f*ts[,1]) + coeffper1[6]*cos(4*pi*f*ts[,1]) + coeffper1[7]*sin(6*pi*f*ts[,1]) +
                              coeffper1[8]*cos(6*pi*f*ts[,1]) + coeffper1[9]*sin(8*pi*f*ts[,1]) + coeffper1[10]*cos(8*pi*f*ts[,1]))
        
        postvar <- var(ts[,2])
        
        varratio <- postvar / prevar
        
        features[33] <- varratio
        
        if (quiet == FALSE){
          
          print(paste("Variance Ratio before and after harmonic subtraction is ", varratio, ".", sep=""))
          
        }
        
        features[4] <- f
        
        features[5] <- pv
        
        if (quiet == FALSE){
          
          print("Calculating #2 and #3 periods independent from period #1...")
          
        }
        
        f23 <- rep.int(0, 2)
        
        pv23 <- rep.int(0, 2)
        
        for (z in 2:3){
          
          if (quiet == FALSE){
            
            print(paste("Calculating the #", z, " period:", sep=""))
            
          }
          
          if (type == "LSP"){
            
            print(paste("Fitting Lomb-Scargle Periodogram...   [", z, "/3]", sep=""))
            
            ans <- lsp(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
            
          }
          else if (type == "LSPL"){
            
            print(paste("Fitting Hybrid-Log Lomb-Scargle Periodogram...   [", z, "/3]", sep=""))
            
            ans <- lspl(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
            
          }
          else if (type == "GLS"){
            
            print(paste("Fitting Generalised Lomb-Scargle Periodogram...   [", z, "/3]", sep=""))
            
            ans <- GLS(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
            
          }
          else if (type == "VRP"){
            
            print(paste("Fitting Variance Ratio Periodogram...   [", z, "/3]", sep=""))
            
            ans <- VRP(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
            
          }
          else if (type == "SLLK"){
            
            print(paste("Fitting String Length Lafler-Kinman Periodogram...   [", z, "/3]", sep=""))
            
            ans <- SLLK(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
            
          }
          else if (type == "BKR"){
            
            print(paste("Fitting Blum-Kiefer-Rosenblatt Periodogram...   [", z, "/3]", sep=""))
            
            ans <- BKR(ts, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
            
          }
          else if (type == "CKP"){
            
            print(paste("Fitting Correntropy Kernelised Periodogram...   [", z, "/3]", sep=""))
            
            ans <- CKP(ts, IP, filt = spurs, ofac = ovsm, from = freqmin, to = freqmax, nt = nt, b = b, plot = FALSE)
            
          }
          
          pwr <- ans$power
          
          if (quiet == FALSE){
            
            plot(ans$scanned, pwr, pch=19, type="n", main = paste("Prewhitening #", z-1, " '", type, "' Periodogram for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
                 ylab = paste("'", type, "' Statistic", sep=""))
            
            lines(ans$scanned, pwr, type="l")
            
          }
          
          f <- ans$scanned[order(pwr, decreasing = T)]
          
          for (zz in 1:length(spurs)){
            
            f <- f[which(!(abs((1 / f) - (1 / spurs[zz])) < eps*(1 / spurs[zz])))]
            
            f <- f[which(!(((1 / f) > (1 / spurs[zz])) & (abs(((1 / f) / (1 / spurs[zz])) - floor(((1 / f) / (1 / spurs[zz])))) < eps)))]
            
            f <- f[which(!(((1 / f) < (1 / spurs[zz])) & (abs(((1 / spurs[zz]) / (1 / f)) - floor(((1 / spurs[zz]) / (1 / f)))) < eps)))]
            
          }
          
          fo <- f[1]
          
          fine <- fineper(ts, f, type)
          
          f <- fine$freq
          
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
          
          p <- 1 / f
          
          f23[z-1] <- f
          
          pv23[z-1] <- pv
          
          if (quiet == FALSE){
            
            print(paste("Displaying #", z, " period detected:", sep=""))
            
            print(p)
            
            print(paste("Displaying #", z, " period p-values:", sep=""))
            
            print(pv)
            
            print(paste("Calculating best 1st to 4th harmonic fit for the #", z, " period...", sep=""))
            
          }
          
          tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
          
          tsin <- rbind(tsin, ts)
          
          tsin <- tsin[order(tsin[,1]),]
          
          SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                        sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                        cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
          
          coeffper2 <- SSTlm$coefficients
          
          coeffper2[is.na(coeffper2)] <- 0
          
          if (z == 2){
            
            if (quiet == FALSE){
              
              print("Subtracting harmonic fit of period #2 from light curve (prewhitening #2)...")
              
            }
            
            ts[,2] <- ts[,2] - (coeffper2[2]*ts[,1] + coeffper2[3]*sin(2*pi*f*ts[,1]) + coeffper2[4]*cos(2*pi*f*ts[,1]) +
                                  coeffper2[5]*sin(4*pi*f*ts[,1]) + coeffper2[6]*cos(4*pi*f*ts[,1]) + coeffper2[7]*sin(6*pi*f*ts[,1]) +
                                  coeffper2[8]*cos(6*pi*f*ts[,1]) + coeffper2[9]*sin(8*pi*f*ts[,1]) + coeffper2[10]*cos(8*pi*f*ts[,1]))
            
          }
          
        }
        
        f <- features[4]
        
        p <- 1 / f
        
        if (quiet == FALSE){
          
          print("Calculating harmonic best-fit to the original light curve...")
          
        }
        
        ts <- tsts
        
        tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
        
        tsin <- rbind(tsin, ts)
        
        tsin <- tsin[order(tsin[,1]),]
        
        SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                      sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                      cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]) + sin(2*pi*f23[1]*tsin[,1]) + cos(2*pi*f23[1]*tsin[,1]) +
                      sin(4*pi*f23[1]*tsin[,1]) + cos(4*pi*f23[1]*tsin[,1]) + sin(6*pi*f23[1]*tsin[,1]) +
                      cos(6*pi*f23[1]*tsin[,1]) + sin(8*pi*f23[1]*tsin[,1]) + cos(8*pi*f23[1]*tsin[,1]) +
                      sin(2*pi*f23[2]*tsin[,1]) + cos(2*pi*f23[2]*tsin[,1]) + sin(4*pi*f23[2]*tsin[,1]) + cos(4*pi*f23[2]*tsin[,1]) + sin(6*pi*f23[2]*tsin[,1]) +
                      cos(6*pi*f23[2]*tsin[,1]) + sin(8*pi*f23[2]*tsin[,1]) + cos(8*pi*f23[2]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
        #print(summary(SSTlm))
        
        co <- SSTlm$coefficients
        
        co[is.na(co)] <- 0
        
        if (quiet == FALSE){
          
          print("Plotting Light Curve folded onto #1 period...")
          
        }
        
        pos <- weighted.mean(ts[,2], (1 / ts[,3])) + rang * sd(ts[,2])
        
        neg <- weighted.mean(ts[,2], (1 / ts[,3])) - rang * sd(ts[,2])
        
        folded <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (1 / features[4])) + 1.25) %% 2) - 1
        
        folded2 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (1 / features[4])) + 0.25) %% 1)
        
        folded_2 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (1 / f23[1])) + 1.25) %% 2) - 1
        
        folded2_2 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (1 / f23[1])) + 0.25) %% 1)
        
        folded_3 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (1 / f23[2])) + 1.25) %% 2) - 1
        
        folded2_3 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / (1 / f23[2])) + 0.25) %% 1)
        
        fts <- matrix(c(folded2, ts[,2]), nrow = length(folded2))
        
        fts <- fts[order(fts[,1]),]
        
        pcs <- (1 / (n * stddev)) * cumsum(fts[,2] - meanmag)
        
        features[44] <- max(pcs) - min(pcs)
        
        varfold <- (sd(fts[,2])) ^ 2.0
        
        features[37] <-  (1.0 / ((n - 1) * varfold) * sum((fts[2:n,2] - fts[1:(n-1),2]) ^ 2.0))
        
        folded3 <- folded2 - 1
        
        folded4 <- c(folded3, folded2)
        
        foldts <- c(ts[,2], ts[,2])
        
        folded3_2 <- folded2_2 - 1
        
        folded4_2 <- c(folded3_2, folded2_2)
        
        folded3_3 <- folded2_3 - 1
        
        folded4_3 <- c(folded3_3, folded2_3)
        
        if (quiet == FALSE){
          
          plot(folded, ts[,2], main = paste("", objects$usnoref[i], " Folded Light Curve spread from -1.0 to 1.0 for first period", sep=""), xlab = paste("Phase at the first period of ", (1 / features[4]), " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          if (ploterr == TRUE){
            
            magerr <- ts[,3]
            
            segments(folded, ts[,2]-magerr, folded, ts[,2]+magerr)
            
            epsilon = 0.01
            
            segments(folded-epsilon, ts[,2]-magerr, folded+epsilon, ts[,2]-magerr)
            
            segments(folded-epsilon, ts[,2]+magerr, folded+epsilon, ts[,2]+magerr)
            
          }
          
          plot(folded4, foldts, main = paste("", objects$usnoref[i], " Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0 for first period", sep=""), xlab = paste("Phase at the first period of ", (1 / features[4]), " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          if (ploterr == TRUE){
            
            magerr <- c(ts[,3], ts[,3])
            
            segments(folded4, foldts-magerr, folded4, foldts+magerr)
            
            epsilon = 0.01
            
            segments(folded4-epsilon, foldts-magerr, folded4+epsilon, foldts-magerr)
            
            segments(folded4-epsilon, foldts+magerr, folded4+epsilon, foldts+magerr)
            
          }
          
          plot(folded_2, ts[,2], main = paste("", objects$usnoref[i], " Folded Light Curve spread from -1.0 to 1.0 for second period", sep=""), xlab = paste("Phase at the second period of ", (1 / f23[1]), " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          if (ploterr == TRUE){
            
            magerr <- ts[,3]
            
            segments(folded_2, ts[,2]-magerr, folded_2, ts[,2]+magerr)
            
            epsilon = 0.01
            
            segments(folded_2-epsilon, ts[,2]-magerr, folded_2+epsilon, ts[,2]-magerr)
            
            segments(folded_2-epsilon, ts[,2]+magerr, folded_2+epsilon, ts[,2]+magerr)
            
          }
          
          plot(folded4_2, foldts, main = paste("", objects$usnoref[i], " Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0 for second period", sep=""), xlab = paste("Phase at the second period of ", (1 / f23[1]), " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          if (ploterr == TRUE){
            
            magerr <- c(ts[,3], ts[,3])
            
            segments(folded4_2, foldts-magerr, folded4_2, foldts+magerr)
            
            epsilon = 0.01
            
            segments(folded4_2-epsilon, foldts-magerr, folded4_2+epsilon, foldts-magerr)
            
            segments(folded4_2-epsilon, foldts+magerr, folded4_2+epsilon, foldts+magerr)
            
          }
          
          plot(folded_3, ts[,2], main = paste("", objects$usnoref[i], " Folded Light Curve spread from -1.0 to 1.0 for third period", sep=""), xlab = paste("Phase at the third period of ", (1 / f23[2]), " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          if (ploterr == TRUE){
            
            magerr <- ts[,3]
            
            segments(folded_3, ts[,2]-magerr, folded_3, ts[,2]+magerr)
            
            epsilon = 0.01
            
            segments(folded_3-epsilon, ts[,2]-magerr, folded_3+epsilon, ts[,2]-magerr)
            
            segments(folded_3-epsilon, ts[,2]+magerr, folded_3+epsilon, ts[,2]+magerr)
            
          }
          
          plot(folded4_3, foldts, main = paste("", objects$usnoref[i], " Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0 for third period", sep=""), xlab = paste("Phase at the third period of ", (1 / f23[2]), " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          if (ploterr == TRUE){
            
            magerr <- c(ts[,3], ts[,3])
            
            segments(folded4_3, foldts-magerr, folded4_3, foldts+magerr)
            
            epsilon = 0.01
            
            segments(folded4_3-epsilon, foldts-magerr, folded4_3+epsilon, foldts-magerr)
            
            segments(folded4_3-epsilon, foldts+magerr, folded4_3+epsilon, foldts+magerr)
            
          }
          
          print("Plotting Raw Light Curve with fitted model...")
          
          plot(ts[,1], ts[,2], main = paste("", objects$usnoref[i], " Raw Light Curve with fitted model", sep=""), xlab = "Modified Julian Date (Days)",
               ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
          
          lines(sort(tsin[,1]), SSTlm$fitted.values, col = 2)
          
          if (ploterr == TRUE){
            
            magerr <- ts[,3]
            
            segments(ts[,1], ts[,2]-magerr, ts[,1], ts[,2]+magerr)
            
            epsilon2 = 4
            
            segments(ts[,1]-epsilon2, ts[,2]-magerr, ts[,1]+epsilon2, ts[,2]-magerr)
            
            segments(ts[,1]-epsilon2, ts[,2]+magerr, ts[,1]+epsilon2, ts[,2]+magerr)
            
          }
          
          print("Plotting Residuals from fitted model...")
          
          res <- ts[,2] - (co[1] + co[2]*ts[,1] + co[3]*sin(2*pi*features[4]*ts[,1]) + co[4]*cos(2*pi*features[4]*ts[,1]) +
                             co[5]*sin(4*pi*features[4]*ts[,1]) + co[6]*cos(4*pi*features[4]*ts[,1]) + co[7]*sin(6*pi*features[4]*ts[,1]) +
                             co[8]*cos(6*pi*features[4]*ts[,1]) + co[9]*sin(8*pi*features[4]*ts[,1]) + co[10]*cos(8*pi*features[4]*ts[,1]) + 
                             co[3]*sin(2*pi*f23[1]*ts[,1]) + co[4]*cos(2*pi*f23[1]*ts[,1]) +
                             co[5]*sin(4*pi*f23[1]*ts[,1]) + co[6]*cos(4*pi*f23[1]*ts[,1]) + co[7]*sin(6*pi*f23[1]*ts[,1]) +
                             co[8]*cos(6*pi*f23[1]*ts[,1]) + co[9]*sin(8*pi*f23[1]*ts[,1]) + co[10]*cos(8*pi*f23[1]*ts[,1]) +
                             co[3]*sin(2*pi*f23[2]*ts[,1]) + co[4]*cos(2*pi*f23[2]*ts[,1]) +
                             co[5]*sin(4*pi*f23[2]*ts[,1]) + co[6]*cos(4*pi*f23[2]*ts[,1]) + co[7]*sin(6*pi*f23[2]*ts[,1]) +
                             co[8]*cos(6*pi*f23[2]*ts[,1]) + co[9]*sin(8*pi*f23[2]*ts[,1]) + co[10]*cos(8*pi*f23[2]*ts[,1]))
          
          plot(ts[,1], res, main = paste("", objects$usnoref[i], " Residuals from fitted model", sep=""), xlab = "Modified Julian Date (Days)",
               ylab = "Residual Magnitude", pch=19, ylim = c(max(res), min(res)))
          
          radnt <- ad.test(res)
          
          resnorm <- as.numeric(radnt$statistic)
          
          print(paste("Normality Test on Residuals (0.3 ~ normal, larger is non-normal): ", resnorm, sep=""))
          
        }
        
        features[2] <- 0
        
        if (any(SSTlm$fitted.values > pos) | any(SSTlm$fitted.values < neg)){
          
          features[2] <- 1
          
        }
        
        if (quiet == FALSE){
          
          print(paste("Regression Flag = ", features[2], "", sep=""))
          
        }
        
        features[3] <- co[2]
        
        features[6] <- f23[1]
        
        features[7] <- pv23[1]
        
        features[8] <- f23[2]
        
        features[9] <- pv23[2]
        
        features[10] <- (co[3]^2.0 + co[4]^2.0)^0.5
        
        features[11] <- (co[5]^2.0 + co[6]^2.0)^0.5
        
        features[12] <- (co[7]^2.0 + co[8]^2.0)^0.5
        
        features[13] <- (co[9]^2.0 + co[10]^2.0)^0.5
        
        features[14] <- (co[11]^2.0 + co[12]^2.0)^0.5
        
        features[15] <- (co[13]^2.0 + co[14]^2.0)^0.5
        
        features[16] <- (co[15]^2.0 + co[16]^2.0)^0.5
        
        features[17] <- (co[17]^2.0 + co[18]^2.0)^0.5
        
        features[18] <- (co[19]^2.0 + co[20]^2.0)^0.5
        
        features[19] <- (co[21]^2.0 + co[22]^2.0)^0.5
        
        features[20] <- (co[23]^2.0 + co[24]^2.0)^0.5
        
        features[21] <- (co[25]^2.0 + co[26]^2.0)^0.5
        
        PH11 <- atan2(co[4], co[3])
        
        ratfreq2 <- features[6] / features[4]
        
        ratfreq3 <- features[8] / features[4]
        
        PH12 <- atan2(co[6], co[5]) - (2.0 * PH11)
        
        PH13 <- atan2(co[8], co[7]) - (3.0 * PH11)
        
        PH14 <- atan2(co[10], co[9]) - (4.0 * PH11)
        
        PH21 <- atan2(co[12], co[11]) - (ratfreq2 * PH11)
        
        PH22 <- atan2(co[14], co[13]) - (2.0 * ratfreq2 * PH11)
        
        PH23 <- atan2(co[16], co[15]) - (3.0 * ratfreq2 * PH11)
        
        PH24 <- atan2(co[18], co[17]) - (4.0 * ratfreq2 * PH11)
        
        PH31 <- atan2(co[20], co[19]) - (ratfreq3 * PH11)
        
        PH32 <- atan2(co[22], co[21]) - (2.0 * ratfreq3 * PH11)
        
        PH33 <- atan2(co[24], co[23]) - (3.0 * ratfreq3 * PH11)
        
        PH34 <- atan2(co[26], co[25]) - (4.0 * ratfreq3 * PH11)
        
        features[22] <- atan2(sin(PH12), cos(PH12))
        
        features[23] <- atan2(sin(PH13), cos(PH13))
        
        features[24] <- atan2(sin(PH14), cos(PH14))
        
        features[25] <- atan2(sin(PH21), cos(PH21))
        
        features[26] <- atan2(sin(PH22), cos(PH22))
        
        features[27] <- atan2(sin(PH23), cos(PH23))
        
        features[28] <- atan2(sin(PH24), cos(PH24))
        
        features[29] <- atan2(sin(PH31), cos(PH31))
        
        features[30] <- atan2(sin(PH32), cos(PH32))
        
        features[31] <- atan2(sin(PH33), cos(PH33))
        
        features[32] <- atan2(sin(PH34), cos(PH34))
        
        featvec[i,] <- features
        
      }
      
    }
    
    if (nonper == FALSE){
      
      featvec <- as.data.frame(featvec)
      
      featvec <- cbind(objects$usnoref, featvec)
      
      colnames(featvec, do.NULL = FALSE)
      
      colnames(featvec) <- c("USnoref", "# Observations", "Regression Flag", "Slope of Linear Trend", "#1 Freq", "#1 p-value", "#2 Freq",
                             "#2 p-value", "#3 Freq", "#3 p-value", "A11", "A12", "A13", "A14", "A21", "A22", "A23",
                             "A24", "A31", "A32", "A33", "A34", "PH12", "PH13", "PH14", "PH21", "PH22", "PH23", "PH24",
                             "PH31", "PH32", "PH33", "PH34", "Variance Ratio", "Mean Magnitude", "Mean Variance", "Variability Index",
                             "Folded Variability Index", "Q31", "Skewness", "Small Kurtosis", "StD", "Renyi's Quadratic Entropy", "Rcs", "Psi-cs", "Beyond 1 StD",
                             "Stetson Kurtosis Measure", "Stetson-K AC", "Maximum Slope", "Amplitude", "Median Absolute Deviation",
                             "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20", "Flux Percentile Ratio mid35",
                             "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                             "Percent Amplitude", "Percent Difference Flux Percentile", "Anderson-Darling Statistic",
                             "Autocorrelation Length", "Slotted Autocorrelation Length")
      
    }
    
    else{
      
      featvec <- as.data.frame(featvec)
      
      featvec <- cbind(objects$usnoref, featvec)
      
      featvec <- featvec[,c(1:2, 35:37, 39:44, 46:63)]
      
      colnames(featvec, do.NULL = FALSE)
      
      colnames(featvec) <- c("USnoref", "# Observations", "Mean Magnitude", "Mean Variance", "Variability Index",
                             "Q31", "Skewness", "Small Kurtosis", "StD", "Renyi's Quadratic Entropy", "Rcs", "Beyond 1 StD",
                             "Stetson Kurtosis Measure", "Stetson-K AC", "Maximum Slope", "Amplitude", "Median Absolute Deviation",
                             "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20", "Flux Percentile Ratio mid35",
                             "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                             "Percent Amplitude", "Percent Difference Flux Percentile", "Anderson-Darling Statistic",
                             "Autocorrelation Length", "Slotted Autocorrelation Length")
      
    }
    
    close(channel)
    
    #sparkR.stop()
    
    print(Sys.time() - start)
    
    gc()
    
    if (geo == TRUE){
      
      featvec <- featvec[1,]
      
    }
    
    featvec
    
  }
  
}