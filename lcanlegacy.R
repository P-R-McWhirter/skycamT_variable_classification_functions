lcanlegacy <- function(ra, declin, radii = 0.01, lim = 2, ovsm = 2, run = 5000, type = "LSP", rang = 2, rand = FALSE, ploterr = FALSE, alldb = FALSE, filtper = FALSE, usenyq = TRUE) {
  
  #Sys.setenv(SPARK_HOME = "C:/Apache/spark-1.6.1")
  
  #.libPaths(c(file.path(Sys.getenv("SPARK_HOME"), "R", "lib"), .libPaths()))
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  #library(SparkR)
  
  #sparkR.stop()
  
  start <- Sys.time()
  
  #sc <- sparkR.init(master = "spark://150.204.48.130:7077")
  
  #sqlContext <- sparkRSQL.init(sc)
  
  set.seed(20)
  
  if (lim < 2){
    
    print("Limit too low, changing to '2' and continuing...")
    
    lim = 2
    
  }
  
  radra <- radii / abs(cos(declin))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  if (alldb == TRUE){
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries FROM objdat WHERE entries > '", lim, "'", sep=""))
    
  }
  
  else {
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries FROM objdat WHERE RA BETWEEN '", 
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
    
    if (alldb == TRUE){
      
      print(paste("Reading in all ", length(objects$usnoref), " object(s) with greater than ", lim, " observations.", sep=""))
      
    }
    
    else{
      
      print(paste("A total of ", length(objects$usnoref),
                  " object(s) identified within the search radius of ", radii, " degree(s) with greater than ", lim, " observations.", sep=""))
      
    }
    
    featvec <- matrix(0, length(objects$entries), 63)
    
    for (i in 1:length(objects$entries)){
      
      print("Reading next object from database...")
      
      info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                      objects$usnoref[i], "'", sep=""))
      
      print(paste("Object with reference id '", objects$usnoref[i],
                  "' successfully read in. Object has ", objects$entries[i], " observation(s).", sep=""))
      
      ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
      
      ts <- ts[order(ts[,1]),]
      
      if (run > 0 & run < as.numeric(objects$entries[i])){
        
        print(paste("Number of data points greater than limit. Setting number of data points to ", run, ".", sep=""))
        
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
      
      gapsd <- sd(gap)
      
      gap2 <- gap[gap < 2*gapsd]
      
      plot(1:length(gap), gap, pch=19, type="n", main = paste("The Gaps between observations for ", objects$usnoref[i], "", sep=""),
           xlab = "Interval Number", ylab = "Interval Length (days)")
      
      lines(1:length(gap), gap, type="l")
      
      plot(1:length(gapn), gapn, pch=19, type="n", main = paste("The deviation of uneven observations from the sampling rate (every minute) for ", objects$usnoref[i], "", sep=""),
           xlab = "Interval Number", ylab = "Observational Uneven deviations (seconds)")
      
      lines(1:length(gapn), gapn, type="l")
      
      hist(ts[,1], breaks = as.integer(permax / 10), main = paste("Histogram of the distribution of datapoints in time for ", objects$usnoref[i], "", sep=""),
           xlab = "Modified Julian Date (Days)")
      
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
      
      print(paste("Searching over a period range of ", permin,
                  " days to ", ceiling(permax), " days, the time-span of observations.", sep=""))
      
      print(paste("This is a frequency range of ", freqmin, " days^-1 to ", freqmax, " days^-1.", sep=""))
      
      df <- freqmin * (1 / ovsm)
      
      print(paste("The frequency step is ", df,
                  " days^-1 as defined by the oversampling rate = ", ovsm, ".", sep=""))
      
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
      
      print(paste("The frequency grid is ", freqlen, " candidates in size.", sep=""))
      
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
      
      features <- rep(0, 63)
      
      #fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts), weights = 1/info$Rcaterr)
      
      #coeff <- coefficients(fit)
      
      features[1] <- objects$entries[i]
      
      tsts <- ts
      
      #ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
      
      if (type == "LSP"){
        
        print("Fitting Lomb-Scargle Periodogram...")
        
        ans <- lsp(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
        
      }
      else if (type == "LSPL"){
        
        print("Fitting Hybrid-Log Lomb-Scargle Periodogram...")
        
        ans <- lspl(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
        
      }
      else if (type == "SLLK"){
        
        print("Fitting String Length Lafler-Kinman Periodogram...")
        
        ans <- SLLK(ts[,1:2], from = permin, to = permax, ofac = ovsm, plot = FALSE)
        
      }
      else if (type == "BKR"){
        
        print("Fitting Blum-Kiefer-Rosenblatt Periodogram...")
        
        ans <- BKR(ts[,1:2], from = permin, to = permax, ofac = ovsm, plot = FALSE)
        
      }
      else
      {
        
        close(channel)
        
        stop(paste("'", type, "' is an invalid Periodogram. Program will now exit.", sep=""))
        
      }
      
      if (rand == TRUE){
        
        print("Using random dataset to remove sampling periodicities...")
        
        tsran <- matrix(c(ts[,1], rnorm(length(ts[,2]), mean = weighted.mean(ts[,2], (1/ts[,3])), sd = sd(ts[,2]))), nrow = length(ts[,2]))
        
        if (type == "LSP"){
          
          ansran <- lsp(tsran, from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "LSPL"){
          
          ansran <- lspl(tsran, from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "SLLK"){
          
          ansran <- SLLK(tsran, from = permin, to = permax, ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "BKR"){
          
          ansran <- BKR(tsran, from = permin, to = permax, ofac = ovsm, plot = FALSE)
          
        }
        
        pwr <- ans$power
        
        pwrran <- ansran$power[order(pwr, decreasing = T)]
        
        pwr <- pwr - pwrran
        
      }
      else
      {
        
        print("Not using random dataset. Variable 'rand' set to FALSE.")
        
        pwr <- ans$power
        
      }
      
      if (rand == TRUE){
        
        plot(ans$scanned, pwr, pch=19, type="n", main = paste("'", type, "' Periodogram with randomised data Power subtracted for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
             ylab = "Normalised Power")
        
        lines(ans$scanned, pwr, type="l")
        
      }
      
      else{
        
        plot(ans$scanned, pwr, pch=19, type="n", main = paste("'", type, "' Periodogram for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
             ylab = "Normalised Power")
        
        lines(ans$scanned, pwr, type="l")
        
      }
      
      # Possibly add low period penalty term here log(f/fn) where fn = 0.5<1/delT>
      
      print("Displaying #1 period detected:")
      
      f <- ans$scanned[order(pwr, decreasing = T)]
      
      if (filtper == TRUE){
        
        f <- f[findInterval(f, c(0.95238,1.05263)) != 1L & findInterval(f, c(1.81818,2.22222)) != 1L & findInterval(f, c(0.03378,0.03401)) != 1L]
        
      }
      
      f <- f[1]
      
      p <- 1 / f
      
      effm <- 2 * length(ans$scanned) / ovsm
      
      level <- -log(1 - (1 - ans$alpha)^(1/effm))
      
      exPN <- exp(-1.0 * pwr[ans$scanned == f])
      
      pv <- effm * exPN
      
      if (pv > 0.01){
        
        pv <- 1 - (1 - exPN)^effm
        
      }
      
      print(p)
      
      print("Displaying #1 period p-values:")
      
      print(pv)
      
      stddev <- sd(ts[,2])
      
      features[40] <- stddev
      
      features[38] <- skewness(ts[,2])
      
      n <- length(ts[,2])
      
      meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
      
      w <- rep(0, (n-1))
      
      wn <- 0
      
      for (k in 1:(n-1)){
        
        w[k] <- 1 / ((ts[(k+1),1] - ts[k,1]) ^ 2)
        
        if (ts[(k+1),1] == ts[k,1]){
          
          w[k] <- 0
          
        }
        
        wn <- wn + (w[k] * ((ts[(k+1),2] - ts[k,2]) ^ 2))
        
      }
      
      features[33] <- meanmag
      
      features[34] <- stddev / meanmag
      
      features[35] <- mean(w) * ((ts[(n-1),1] - ts[1,1]) ^ 2) * (wn / ((stddev ^ 2) * sum(w) * n^2))
      
      quan <- quantile(ts[,2])
      
      features[37] <- as.numeric(quan[4]) - as.numeric(quan[2])
      
      features[39] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((ts[,2] - meanmag) / stddev) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
      
      print(paste("Standard Deviation of the Apparent Magnitude: ", stddev, sep=""))
      
      cs <- (1 / (n * stddev)) * cumsum(ts[,2] - meanmag)
      
      features[41] <- max(cs) - min(cs)
      
      features[43] <- (length(ts[,2][ts[,2] > (meanmag + stddev)]) + length(ts[,2][ts[,2] < (meanmag - stddev)])) / length(ts[,2])
      
      sigmap <- (sqrt(n / (n - 1)) * (ts[,2] - meanmag) / (1 / ts[,3]))
      
      features[44] <- (1 / sqrt(n) * sum(abs(sigmap)) / sqrt(sum(sigmap ^ 2)))
      
      features[46] <- 0
      
      for (k in 2:n){
        
        slope <- (ts[k,2] - ts[(k-1),2]) / (ts[k,1] - ts[(k-1),1])
        
        if ((abs(slope) > abs(features[46])) & !is.na(slope)){
          
          features[46] <- slope
          
        }
        
      }
      
      #amplitude <- (max(ts[,2]) - min(ts[,2])) / 2
      
      amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
      
      features[47] <- amplitude
      
      features[48] <- median(abs(ts[,2] - median(ts[,2])))
      
      medianmag <- median(ts[,2])
      
      features[49] <- 1 - ((length(ts[,2][ts[,2] > (medianmag + (amplitude/10))]) + length(ts[,2][ts[,2] < (medianmag - (amplitude/10))])) / length(ts[,2]))
      
      if (n >= 30){
        
        datalast <- ts[,2][length(ts[,2]) - (29:0)]
        
        features[50] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / 30)
        
      }
      else
      {
        
        datalast <- ts[,2][length(ts[,2]) - ((n-1):0)]
        
        features[50] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / n)
        
      }
      
      sortedmags <- sort(ts[,2])
      
      F_60_index <- ceiling(0.60 * n)
      
      F_40_index <- ceiling(0.40 * n)
      
      F_5_index <- ceiling(0.05 * n)
      
      F_95_index <- ceiling(0.95 * n)
      
      F_40_60 <- sortedmags[F_60_index] - sortedmags[F_40_index]
      
      F_5_95 <- sortedmags[F_95_index] - sortedmags[F_5_index]
      
      features[51] <- F_40_60 / F_5_95
      
      F_325_index <- ceiling(0.325 * n)
      
      F_675_index <- ceiling(0.675 * n)
      
      F_325_675 <- sortedmags[F_675_index] - sortedmags[F_325_index]
      
      features[52] <- F_325_675 / F_5_95
      
      F_25_index <- ceiling(0.25 * n)
      
      F_75_index <- ceiling(0.75 * n)
      
      F_25_75 <- sortedmags[F_75_index] - sortedmags[F_25_index]
      
      features[53] <- F_25_75 / F_5_95
      
      F_175_index <- ceiling(0.175 * n)
      
      F_825_index <- ceiling(0.825 * n)
      
      F_175_825 <- sortedmags[F_825_index] - sortedmags[F_175_index]
      
      features[54] <- F_175_825 / F_5_95
      
      F_10_index <- ceiling(0.10 * n)
      
      F_90_index <- ceiling(0.90 * n)
      
      F_10_90 <- sortedmags[F_90_index] - sortedmags[F_10_index]
      
      features[55] <- F_10_90 / F_5_95
      
      distmed <- abs(ts[,2] - medianmag)
      
      maxdist <- max(distmed)
      
      features[56] <- maxdist / medianmag
      
      features[57] <- F_5_95 / medianmag
      
      adnt <- ad.test(ts[,2])
      
      features[58] <- (1 / (1.0 + exp(-10.0 * (as.numeric(adnt$statistic) - 0.3))))
      
      lagval <- 100
      
      autocorlen <- 0
      
      while (autocorlen == 0){
        
        AC <- acf(ts[,2], lag.max = lagval, type = "correlation", plot = FALSE)
        
        racf <- as.numeric(AC$acf)
        
        autocorlen <- match(max(racf[racf < exp(-1)]), racf)
        
        features[59] <- autocorlen
        
        if (is.na(autocorlen)){
          
          autocorlen <- 0
          
        }
        
        lagval <- lagval + 100
        
      }
      
      nac <- length(racf)
      
      sigmapac <- sqrt(nac / (nac - 1)) * ((racf - mean(racf)) / sd(racf))
      
      features[45] <- (1 / sqrt(n) * sum(abs(sigmapac)) / sqrt(sum(sigmapac ^ 2)))
      
      features[60] <- 0 #sacf(ts)
      
      features[61] <- 0
      
      features[62] <- 0
      
      magerr <- ts[,3]
      
      tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(meanmag, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
      
      tsin <- rbind(tsin, ts)
      
      tsin <- tsin[order(tsin[,1]),]
      
      print("Calculating best 1st to 4th harmonic fit for the #1 period...")
      
      SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                    sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                    cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      coeffper1 <- SSTlm$coefficients
      
      coeffper1[is.na(coeffper1)] <- 0
      
      prevar <- var(ts[,2])
      
      print("Subtracting harmonic fit of period #1 from light curve (prewhitening)...")
      
      ts[,2] <- ts[,2] - (SSTlm$coefficients[2]*ts[,1] + SSTlm$coefficients[3]*sin(2*pi*f*ts[,1]) + SSTlm$coefficients[4]*cos(2*pi*f*ts[,1]) +
                            SSTlm$coefficients[5]*sin(4*pi*f*ts[,1]) + SSTlm$coefficients[6]*cos(4*pi*f*ts[,1]) + SSTlm$coefficients[7]*sin(6*pi*f*ts[,1]) +
                            SSTlm$coefficients[8]*cos(6*pi*f*ts[,1]) + SSTlm$coefficients[9]*sin(8*pi*f*ts[,1]) + SSTlm$coefficients[10]*cos(8*pi*f*ts[,1]))
      
      postvar <- var(ts[,2])
      
      varratio <- postvar / prevar
      
      features[32] <- varratio
      
      print(paste("Variance Ratio before and after harmonic subtraction is ", varratio, ".", sep=""))
      
      features[3] <- f
      
      features[4] <- pv
      
      print("Calculating #2 and #3 periods independent from period #1...")
      
      f23 <- rep.int(0, 2)
      
      pv23 <- rep.int(0, 2)
      
      for (z in 2:3){
        
        print(paste("Calculating the #", z, " period:", sep=""))
        
        if (type == "LSP"){
          
          print("Fitting Lomb-Scargle Periodogram...")
          
          ans <- lsp(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "LSPL"){
          
          print("Fitting Hybrid-Log Lomb-Scargle Periodogram...")
          
          ans <- lspl(ts[,1:2], from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "SLLK"){
          
          print("Fitting String Length Lafler-Kinman Periodogram...")
          
          ans <- SLLK(ts[,1:2], from = permin, to = permax, ofac = ovsm, plot = FALSE)
          
        }
        else if (type == "BKR"){
          
          print("Fitting Blum-Kiefer-Rosenblatt Periodogram...")
          
          ans <- BKR(ts[,1:2], from = permin, to = permax, ofac = ovsm, plot = FALSE)
          
        }
        
        if (rand == TRUE){
          
          print("Using random dataset to remove sampling periodicities...")
          
          tsran <- matrix(c(ts[,1], rnorm(length(ts[,2]), mean = weighted.mean(ts[,2], (1/ts[,3])), sd = sd(ts[,2]))), nrow = length(ts[,2]))
          
          if (type == "LSP"){
            
            ansran <- lsp(tsran, from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
            
          }
          else if (type == "LSPL"){
            
            ansran <- lspl(tsran, from = freqmin, to = freqmax, type = "frequency", ofac = ovsm, plot = FALSE)
            
          }
          else if (type == "SLLK"){
            
            ansran <- SLLK(tsran, from = permin, to = permax, ofac = ovsm, plot = FALSE)
            
          }
          else if (type == "BKR"){
            
            ansran <- BKR(tsran, from = permin, to = permax, ofac = ovsm, plot = FALSE)
            
          }
          
          pwr <- ans$power
          
          pwrran <- ansran$power[order(pwr, decreasing = T)]
          
          pwr <- pwr - pwrran
          
        }
        else
        {
          
          print("Not using random dataset. Variable 'rand' set to FALSE.")
          
          pwr <- ans$power
          
        }
        
        if (rand == TRUE){
          
          plot(ans$scanned, pwr, pch=19, type="n", main = paste("Prewhitening #", z-1, " '", type, "' Periodogram with randomised data Power subtracted for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
               ylab = "Normalised Power")
          
          lines(ans$scanned, pwr, type="l")
          
        }
        
        else{
          
          plot(ans$scanned, pwr, pch=19, type="n", main = paste("Prewhitening #", z-1, " '", type, "' Periodogram for ", objects$usnoref[i], "", sep=""), xlab = "Frequency (cycles/day)",
               ylab = "Normalised Power")
          
          lines(ans$scanned, pwr, type="l")
          
        }
        
        print(paste("Displaying #", z, " period detected:", sep=""))
        
        f <- ans$scanned[order(pwr, decreasing = T)]
        
        if (filtper == TRUE){
          
          f <- f[findInterval(f, c(0.95238,1.05263)) != 1L & findInterval(f, c(1.81818,2.22222)) != 1L & findInterval(f, c(0.03378,0.03401)) != 1L]
          
        }
        
        f <- f[1]
        
        effm <- 2 * length(ans$scanned) / ovsm
        
        level <- -log(1 - (1 - ans$alpha)^(1/effm))
        
        exPN <- exp(-1.0 * pwr[ans$scanned == f])
        
        pv <- effm * exPN
        
        if (pv > 0.01){
          
          pv <- 1 - (1 - exPN)^effm
          
        }
        
        p <- 1 / f
        
        f23[z-1] <- f
        
        pv23[z-1] <- pv
        
        print(p)
        
        print(paste("Displaying #", z, " period p-values:", sep=""))
        
        print(pv)
        
        print(paste("Calculating best 1st to 4th harmonic fit for the #", z, " period...", sep=""))
        
        tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
        
        tsin <- rbind(tsin, ts)
        
        tsin <- tsin[order(tsin[,1]),]
        
        SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                      sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                      cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
        
        if (z == 2){
          
          print("Subtracting harmonic fit of period #2 from light curve (prewhitening #2)...")
          
          ts[,2] <- ts[,2] - (SSTlm$coefficients[2]*ts[,1] + SSTlm$coefficients[3]*sin(2*pi*f*ts[,1]) + SSTlm$coefficients[4]*cos(2*pi*f*ts[,1]) +
                                SSTlm$coefficients[5]*sin(4*pi*f*ts[,1]) + SSTlm$coefficients[6]*cos(4*pi*f*ts[,1]) + SSTlm$coefficients[7]*sin(6*pi*f*ts[,1]) +
                                SSTlm$coefficients[8]*cos(6*pi*f*ts[,1]) + SSTlm$coefficients[9]*sin(8*pi*f*ts[,1]) + SSTlm$coefficients[10]*cos(8*pi*f*ts[,1]))
          
        }
        
      }
      
      f <- features[3]
      
      p <- 1 / f
      
      print("Calculating harmonic best-fit to the original light curve...")
      
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
      
      print("Plotting Light Curve folded onto #1 period...")
      
      pos <- weighted.mean(ts[,2], (1 / ts[,3])) + rang * sd(ts[,2])
      
      neg <- weighted.mean(ts[,2], (1 / ts[,3])) - rang * sd(ts[,2])
      
      folded <- ((((ts[,1] - ts[,1][which.max(ts[,1])]) / (1 / features[3])) + 1.25) %% 2) - 1
      
      folded2 <- ((((ts[,1] - ts[,1][which.max(ts[,1])]) / (1 / features[3])) + 0.25) %% 1)
      
      fts <- matrix(c(folded2, ts[,2]), nrow = length(folded2))
      
      fts <- fts[order(fts[,1]),]
      
      pcs <- (1 / (n * stddev)) * cumsum(fts[,2] - meanmag)
      
      features[42] <- max(pcs) - min(pcs)
      
      w <- rep(0, (length(fts[,2])-1))
      
      wn <- 0
      
      for (k in 1:(n-1)){
        
        w[k] <- 1 / ((fts[(k+1),1] - fts[k,1]) ^ 2)
        
        if (fts[(k+1),1] == fts[k,1]){
          
          w[k] <- 0
          
        }
        
        wn <- wn + (w[k] * ((fts[(k+1),2] - fts[k,2]) ^ 2))
        
      }
      
      features[36] <- mean(w) * ((fts[((length(fts[,2]))-1),1] - fts[1,1]) ^ 2) * (wn / ((sd(fts[,2]) ^ 2) * sum(w) * (length(fts[,2]))^2))
      
      folded3 <- folded2 - 1
      
      folded4 <- c(folded3, folded2)
      
      foldts <- c(ts[,2], ts[,2])
      
      plot(folded, ts[,2], main = paste("", objects$usnoref[i], " Folded Light Curve spread from -1.0 to 1.0", sep=""), xlab = paste("Phase at a period of ", (1 / features[3]), " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
      
      if (ploterr == TRUE){
        
        magerr <- ts[,3]
        
        segments(folded, ts[,2]-magerr, folded, ts[,2]+magerr)
        
        epsilon = 0.01
        
        segments(folded-epsilon, ts[,2]-magerr, folded+epsilon, ts[,2]-magerr)
        
        segments(folded-epsilon, ts[,2]+magerr, folded+epsilon, ts[,2]+magerr)
        
      }
      
      plot(folded4, foldts, main = paste("", objects$usnoref[i], " Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0", sep=""), xlab = paste("Phase at a period of ", (1 / features[3]), " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
      
      if (ploterr == TRUE){
        
        magerr <- c(ts[,3], ts[,3])
        
        segments(folded4, foldts-magerr, folded4, foldts+magerr)
        
        epsilon = 0.01
        
        segments(folded4-epsilon, foldts-magerr, folded4+epsilon, foldts-magerr)
        
        segments(folded4-epsilon, foldts+magerr, folded4+epsilon, foldts+magerr)
        
      }
      
      print("Plotting Raw Light Curve with fitted model...")
      
      plot(ts[,1], ts[,2], main = paste("", objects$usnoref[i], " Raw Light Curve with fitted model", sep=""), xlab = "Modified Julian Date (Days)",
           ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
      
      lines(sort(tsin[,1]), SSTlm$fitted.values, col = 2)
      
      features[63] <- 0
      
      if (any(SSTlm$fitted.values > pos) | any(SSTlm$fitted.values < neg)){
        
        features[63] <- 1
        
      }
      
      #lines(tn, yn, col = 4)
      
      if (ploterr == TRUE){
        
        magerr <- ts[,3]
        
        segments(ts[,1], ts[,2]-magerr, ts[,1], ts[,2]+magerr)
        
        epsilon2 = 4
        
        segments(ts[,1]-epsilon2, ts[,2]-magerr, ts[,1]+epsilon2, ts[,2]-magerr)
        
        segments(ts[,1]-epsilon2, ts[,2]+magerr, ts[,1]+epsilon2, ts[,2]+magerr)
        
      }
      
      print(paste("Regression Flag = ", features[63], "", sep=""))
      
      features[2] <- co[2]
      
      features[5] <- f23[1]
      
      features[6] <- pv23[1]
      
      features[7] <- f23[2]
      
      features[8] <- pv23[2]
      
      features[9] <- (co[3]^2.0 + co[4]^2.0)^0.5
      
      features[10] <- (co[5]^2.0 + co[6]^2.0)^0.5
      
      features[11] <- (co[7]^2.0 + co[8]^2.0)^0.5
      
      features[12] <- (co[9]^2.0 + co[10]^2.0)^0.5
      
      features[13] <- (co[11]^2.0 + co[12]^2.0)^0.5
      
      features[14] <- (co[13]^2.0 + co[14]^2.0)^0.5
      
      features[15] <- (co[15]^2.0 + co[16]^2.0)^0.5
      
      features[16] <- (co[17]^2.0 + co[18]^2.0)^0.5
      
      features[17] <- (co[19]^2.0 + co[20]^2.0)^0.5
      
      features[18] <- (co[21]^2.0 + co[22]^2.0)^0.5
      
      features[19] <- (co[23]^2.0 + co[24]^2.0)^0.5
      
      features[20] <- (co[25]^2.0 + co[26]^2.0)^0.5
      
      PH11 <- atan2(co[4], co[3])
      
      ratfreq2 <- features[5] / features[3]
      
      ratfreq3 <- features[7] / features[3]
      
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
      
      features[21] <- atan2(sin(PH12), cos(PH12))
      
      features[22] <- atan2(sin(PH13), cos(PH13))
      
      features[23] <- atan2(sin(PH14), cos(PH14))
      
      features[24] <- atan2(sin(PH21), cos(PH21))
      
      features[25] <- atan2(sin(PH22), cos(PH22))
      
      features[26] <- atan2(sin(PH23), cos(PH23))
      
      features[27] <- atan2(sin(PH24), cos(PH24))
      
      features[28] <- atan2(sin(PH31), cos(PH31))
      
      features[29] <- atan2(sin(PH32), cos(PH32))
      
      features[30] <- atan2(sin(PH33), cos(PH33))
      
      features[31] <- atan2(sin(PH34), cos(PH34))
      
      featvec[i,] <- features
      
    }
    
  }
  
  featvec <- as.data.frame(featvec)
  
  featvec <- cbind(objects$usnoref, featvec)
  
  colnames(featvec, do.NULL = FALSE)
  
  colnames(featvec) <- c("USnoref", "# Observations", "Slope of Linear Trend", "#1 Freq", "#1 p-value", "#2 Freq",
                         "#2 p-value", "#3 Freq", "#3 p-value", "A11", "A12", "A13", "A14", "A21", "A22", "A23",
                         "A24", "A31", "A32", "A33", "A34", "PH12", "PH13", "PH14", "PH21", "PH22", "PH23", "PH24",
                         "PH31", "PH32", "PH33", "PH34", "Variance Ratio", "Mean Magnitude", "Mean Variance", "Variability Index",
                         "Folded Variability Index", "Q31", "Skewness", "Small Kurtosis", "StD", "Rcs", "Psi-cs", "Beyond 1 StD",
                         "Stetson Kurtosis Measure", "Stetson-K AC", "Maximum Slope", "Amplitude", "Median Absolute Deviation",
                         "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20", "Flux Percentile Ratio mid35",
                         "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                         "Percent Amplitude", "Percent Difference Flux Percentile", "Anderson-Darling Statistic",
                         "Autocorrelation Length", "Slotted Autocorrelation Length", "QSO", "Non-QSO", "Regression Flag")
  
  close(channel)
  
  #sparkR.stop()
  
  print(Sys.time() - start)
  
  featvec
  
}