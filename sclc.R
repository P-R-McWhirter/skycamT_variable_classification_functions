sclc <- function(ra, declin, radii = 0.01, lim = 2, pca, spca, sigk = 2.5, lctrenddata, seedno = 20, redo = 5, ploterr = FALSE, alldb = FALSE, geo = FALSE, nonper = FALSE, quiet = FALSE) {
  
  closeAllConnections()
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  
  start <- Sys.time()
  
  if (quiet == FALSE){
    
    print("Executing Software...")
    
  }
  
  set.seed(seedno)
  eps <- 1e-3
  clipcut <- 0.2
  ovsm <- 5
  rang <- 2
  gaufil <- c(1, 100)
  
  if (lim < 2){
    
    print("Limit too low, changing to '2' and continuing...")
    
    lim = 2
    
  }
  
  radra <- radii / abs(cos(declin*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  if (alldb == TRUE){
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries FROM objdat WHERE entries >= '", lim, "'", sep=""))
    
  }
  
  else {
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                       declin+radii, "' AND entries >= '", lim, "'", sep=""))
    
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
    
    if (nonper == TRUE){
      
      featvec <- matrix(0, length(objects$entries), 38)
      
    }
    else{
      
      featvec <- matrix(0, length(objects$entries), 112)
      
    }
    
    if (geo == TRUE){
      
      if (quiet == FALSE){
        
        print("Computing geometrically closest object to the reference coordinates...")
        
      }
      
      geodist <- sqrt((declin - objects$DEClin)^2.0 + (ra - objects$RA)^2.0)
      
      mindist <- which.min(geodist)
      
      objects <- objects[mindist,]
      
    }
    else{
      
      geodist <- sqrt((declin - objects$DEClin)^2.0 + (ra - objects$RA)^2.0)
      
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
      
      if (quiet == FALSE){
        
        plot(1:length(gap), gap, pch=19, type="n", main = paste("The Gaps between observations for ", objects$usnoref[i], "", sep=""),
             xlab = "Interval Number", ylab = "Interval Length (days)")
        
        lines(1:length(gap), gap, type="l")
        
        plot(1:length(gapn), gapn, pch=19, type="n", main = paste("The deviation of uneven observations from the sampling rate (every minute) for ", objects$usnoref[i], "", sep=""),
             xlab = "Interval Number", ylab = "Observational Uneven deviations (seconds)")
        
        lines(1:length(gapn), gapn, type="l")
        
        hist(ts[,1], breaks = max(c(as.integer(permax / 10), 2)), main = paste("Histogram of the distribution of time instants for ", objects$usnoref[i], "", sep=""),
             xlab = "Modified Julian Date (Days)")
        
      }
      
      permin <- 1 / 20
      
      freqmin <- 1 / permax
      
      freqmax <- 20
      
      if (quiet == FALSE){
        
        print(paste("Searching over a period range of ", permin,
                    " days to ", ceiling(permax), " days, the time-span of observations.", sep=""))
        
        print(paste("This is a frequency range of ", freqmin, " days^-1 to ", freqmax, " days^-1.", sep=""))
        
      }
      
      features <- rep(0, 112)
      
      n <- length(ts[,2])
        
      varm <- try(lsp(ts[,1:2], from = 2/3, to = 2, ofac = 1, plot = F), TRUE)
      
      if (class(varm) == "try-error"){
        
        varpeak <- 1
        
        varval <- 0
        
      }
      else{
      
        varpeak <- varm$peak.at[2]
      
        varval <- -log10(varm$p.value)
       
        if (varval == Inf){
          
          varval <- 324
          
        }
        
      }
      
      features[1] <- n
      
      features[2] <- permax
      
      features[3] <- as.numeric(as.vector(geodist[i]))*60*60
      
      if (geo == TRUE){
        
        features[3] <- as.numeric(as.vector(geodist[mindist]))*60*60
        
      }
      
      tsts <- ts
      
      print("Calculating Non-Periodic Features...")
      
      meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
      
      features[4] <- meanmag
      
      stddev <- sd(ts[,2])
      
      features[5] <- stddev / meanmag
      
      features[6] <- median(ts[,2])
      
      features[7] <- sum(abs((ts[,2] - features[6])/ts[,3])) / n
      
      features[12] <- stddev
      
      features[13] <- skewness(ts[,2])
      
      features[8] <- (sum(diff(ts[,2])^2)/(n-1))/stddev^2
      
      features[9] <- sum((ts[,2]-meanmag)^2.0)/sum(diff(ts[,2])^2.0)
      
      w <- 1 / diff(ts[,1])^2.0
      
      wmean <- mean(w)
      
      S1 <- sum(w * (diff(ts[,2])^2.0))
      
      S2 <- sum(w)
      
      features[10] <- (wmean * ((ts[(n-1),1] - ts[1,1]) ^ 2.0) * S1 / ((stddev ^ 2.0) * S2 * (n ^ 2.0)))
      
      quan <- quantile(ts[,2])
      
      features[11] <- as.numeric(quan[4]) - as.numeric(quan[2])
      
      features[14] <- ((n-1)/((n-2)*(n-3)))*((n+1)*(((sum((ts[,2] - meanmag)^4.0)/n)/(sum(((ts[,2] - meanmag)^2.0)/n)^2.0))-3)+6)
      
      features[15] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((ts[,2] - meanmag) / stddev) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
      
      cs <- (1 / (n * stddev)) * cumsum(ts[,2] - meanmag)
      
      features[16] <- max(cs) - min(cs)
      
      constar <- 3
      
      if (n < 3){
        
        features[17] <- 0
        
      }
      else{
        
        concount <- 0
        
        for (k in 1:(n - constar)){
          
          conflag <- 0
          
          for (z in 1:3){
            
            if (ts[k+z,2] > meanmag + 2*stddev | ts[k+z,2] < meanmag - 2*stddev){
              
              conflag <- 1
              
            }
            else{
              
              conflag <- 0
              break
              
            }
            
          }
          
          if (conflag == 1){
            
            concount <- concount + 1
            
          }
          
        }
        
        features[17] <- concount * (1.0 / (n - constar))
        
      }
      
      features[18] <- (length(ts[,2][ts[,2] > (meanmag + stddev)]) + length(ts[,2][ts[,2] < (meanmag - stddev)])) / length(ts[,2])
      
      deln <- (ts[,2] - meanmag) / (1 / ts[,3])
      
      features[19] <- 0
      
      features[20] <- 0
      
      for (k in 1:(n-1)){
        
        features[19] <- features[19] + (deln[k]*deln[k+1])
        
        features[20] <- features[20] + ((sign(deln[k]*deln[k+1]))*sqrt(abs(deln[k]*deln[k+1])))
        
      }
      
      features[19] <- sqrt(1 / ((n - 1)*(n - 2))) * features[19]
      
      sigmap <- (sqrt(n / (n - 1)) * (ts[,2] - meanmag) / (1 / ts[,3]))
      
      features[21] <- (1 / sqrt(n) * sum(abs(sigmap)) / sqrt(sum(sigmap ^ 2)))
      
      features[22] <- (features[20]*features[21])/0.798
      
      features[23] <- 0
      
      for (k in 2:n){
        
        slope <- (ts[k,2] - ts[(k-1),2]) / (ts[k,1] - ts[(k-1),1])
        
        if ((abs(slope) > abs(features[23])) & !is.na(slope)){
          
          features[23] <- slope
          
        }
        
      }
      
      features[23] <- abs(features[23])
      
      amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
      
      features[24] <- amplitude
      
      features[25] <- median(abs(ts[,2] - median(ts[,2])))
      
      medianmag <- median(ts[,2])
      
      features[26] <- 1 - ((length(ts[,2][ts[,2] > (medianmag + (amplitude/10))]) + length(ts[,2][ts[,2] < (medianmag - (amplitude/10))])) / length(ts[,2]))
      
      if (n >= 30){
        
        datalast <- ts[,2][length(ts[,2]) - (29:0)]
        
        features[27] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / 30)
        
      }
      else
      {
        
        datalast <- ts[,2][length(ts[,2]) - ((n-1):0)]
        
        features[27] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / n)
        
      }
      
      sortedmags <- sort(ts[,2])
      
      F_60_index <- ceiling(0.60 * n)
      
      F_40_index <- ceiling(0.40 * n)
      
      F_5_index <- ceiling(0.05 * n)
      
      F_95_index <- ceiling(0.95 * n)
      
      F_40_60 <- sortedmags[F_60_index] - sortedmags[F_40_index]
      
      F_5_95 <- sortedmags[F_95_index] - sortedmags[F_5_index]
      
      features[28] <- F_40_60 / F_5_95
      
      F_325_index <- ceiling(0.325 * n)
      
      F_675_index <- ceiling(0.675 * n)
      
      F_325_675 <- sortedmags[F_675_index] - sortedmags[F_325_index]
      
      features[29] <- F_325_675 / F_5_95
      
      F_25_index <- ceiling(0.25 * n)
      
      F_75_index <- ceiling(0.75 * n)
      
      F_25_75 <- sortedmags[F_75_index] - sortedmags[F_25_index]
      
      features[30] <- F_25_75 / F_5_95
      
      F_175_index <- ceiling(0.175 * n)
      
      F_825_index <- ceiling(0.825 * n)
      
      F_175_825 <- sortedmags[F_825_index] - sortedmags[F_175_index]
      
      features[31] <- F_175_825 / F_5_95
      
      F_10_index <- ceiling(0.10 * n)
      
      F_90_index <- ceiling(0.90 * n)
      
      F_10_90 <- sortedmags[F_90_index] - sortedmags[F_10_index]
      
      features[32] <- F_10_90 / F_5_95
      
      distmed <- abs(ts[,2] - medianmag)
      
      maxdist <- max(distmed)
      
      features[33] <- maxdist / medianmag
      
      features[34] <- F_5_95 / medianmag
      
      lagval <- 100
      
      autocorlen <- 0
      
      while (autocorlen == 0){
        
        AC <- acf(ts[,2], lag.max = lagval, type = "correlation", plot = FALSE)
        
        racf <- as.numeric(AC$acf)
        
        autocorlen <- match(max(racf[racf < exp(-1)]), racf)
        
        features[35] <- autocorlen
        
        if (is.na(autocorlen)){
          
          autocorlen <- 0
          
        }
        
        lagval <- lagval + 100
        
      }
      
      features[36] <- brcol[i]
      
      if (geo == TRUE){
        
        features[36] <- brcol[mindist]
        
      }
      
      features[37] <- varpeak
      
      features[38] <- varval
      
      if (nonper == TRUE){
        
        features <- features[1:38]
        
        featvec[i,] <- features
        
      }
      else if (nonper == FALSE){
        
        print("Fitting using a genetic algorithm...")
        
        noimod <- 1
        
        jit <- 0.6*amplitude
        
        target <- as.data.frame(cbind(as.character(objects$usnoref[i]), " ", "NEW", objects$RA, objects$DEClin, 1.0, " "))
        
        colnames(target) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude")
        
        target <- target[i,]
        
        ans <- try(LCPE(target, pop = 200, ps = 50, nogen = 100, mut = seq(0.8, 0.0, length.out = 100), seed = seedno, redo = redo, quiet = T, fit = "bgls", jit = jit, noimod = noimod, sigk = sigk, seascoe = seascoe, trendrem = trendrem, rangcut = rangcut, radii = 0.00001, lim = lim, tenper = F), TRUE)
        
        if (class(ans) == "try-error"){
          
          print("Genetic Algorithm failed. Likely need to add more white noise jitter...")
          
          jit <- 0.9*amplitude
          
          ans <- try(LCPE(target, pop = 200, ps = 50, nogen = 100, mut = seq(0.8, 0.0, length.out = 100), seed = seedno, redo = redo, quiet = T, fit = "bgls", jit = jit, noimod = noimod, sigk = sigk, seascoe = seascoe, trendrem = trendrem, rangcut = rangcut, radii = 0.00001, lim = lim, tenper = F), TRUE)
          
        }
        
        if (class(ans) == "try-error"){
          
          print("Genetic Algorithm failed. Likely need to add more white noise jitter...")
          
          jit <- 0.9*amplitude
          
          ans <- try(LCPE(target, pop = 200, ps = 50, nogen = 100, mut = seq(0.8, 0.0, length.out = 100), seed = seedno, redo = redo, quiet = T, fit = "bgls", jit = jit, noimod = noimod, sigk = sigk, seascoe = seascoe, trendrem = trendrem, rangcut = rangcut, radii = 0.00001, lim = lim, tenper = F), TRUE)
          
        }
        
        if (class(ans) == "try-error"){
          
          print("Genetic Algorithm failed. Likely need to add even more white noise jitter...")
          
          jit <- 1.2*amplitude
          
          ans <- try(LCPE(target, pop = 200, ps = 50, nogen = 100, mut = seq(0.8, 0.0, length.out = 100), seed = seedno, redo = redo, quiet = T, fit = "bgls", jit = jit, noimod = noimod, sigk = sigk, seascoe = seascoe, trendrem = trendrem, rangcut = rangcut, radii = 0.00001, lim = lim, tenper = F), TRUE)
          
        }
        
        if (class(ans) == "try-error"){
          
          print("Genetic Algorithm failed. Going to 3 mag noisy light curve estimation...")
          
          jit <- 3.0
          
          ans <- LCPE(target, pop = 200, ps = 50, nogen = 100, mut = seq(0.8, 0.0, length.out = 100), seed = seedno, redo = redo, quiet = T, fit = "bgls", jit = jit, noimod = noimod, sigk = sigk, seascoe = seascoe, trendrem = trendrem, rangcut = rangcut, radii = 0.00001, lim = lim, tenper = F)
          
        }
        
        f <- 1 / as.numeric(as.vector(ans$Vuong.Period))
        
        p <- 1 / f
        
        fs <- as.numeric(as.vector(ans$Fit.Stat))
        
        pv <- as.numeric(as.vector(ans$Vuong.Confidence))
        
        pvd <- as.numeric(as.vector(ans$Vuong.Day.Confidence))
        
        ns <- as.numeric(as.vector(ans$Num.Spurious))
        
        nearspur <- as.numeric(as.vector(ans$Nearest.Spurious))
        
        lowper <- 0.99 * p
        
        highper <- 1.01 * p
        
        tuneper <- vartune(ts, from = 1/highper, to = 1/lowper, ofac = 10, plot = F)
        
        tuneper2 <- vartune(ts, from = 0.5/highper, to = 0.5/lowper, ofac = 10, plot = F)
        
        #if (p < 100){
        
        #if (tuneper$peak > tuneper2$peak){
        
        f <- as.numeric(tuneper$peak.at[1])
        
        p <- as.numeric(tuneper$peak.at[2])
        
        #}
        #else{
        
        #f <- as.numeric(tuneper2$peak.at[1])
        
        #p <- as.numeric(tuneper2$peak.at[2])
        
        #}
        
        #}
        
        f2 <- as.numeric(tuneper2$peak.at[1])
        
        p2 <- as.numeric(tuneper2$peak.at[2])
        
        print(paste("Variance Ratio of P: ", tuneper$peak, ". Variance Ratio of 2P: ", tuneper2$peak, ".", sep=""))
        
        perratio <- as.numeric(tuneper$peak) / as.numeric(tuneper2$peak)
        
        print("Displaying #1 period detected:")
        
        print(p)
        
        perspurlvl <- gaufil[1]*exp(-1.0 * ((f - nearspur)^2.0 / 2.0*(gaufil[2]^2.0)))
        
        perspurlvl <- -log10(perspurlvl)
        
        if (perspurlvl == Inf){
          
          perspurlvl <- 324
          
        }
        
        if (quiet == FALSE){
          
          print("Displaying #1 period Vuong Confidence:")
          
          print(pv)
          
        }
        
        magerr <- ts[,3]
        
        tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(meanmag, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
        
        tsin <- rbind(tsin, ts)
        
        tsin <- tsin[order(tsin[,1]),]
        
        if (quiet == FALSE){
          
          print("Calculating harmonic best-fit to the light curve...")
          
        }
        
        ts <- tsts
        
        tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
        
        tsin <- rbind(tsin, ts)
        
        tsin <- tsin[order(tsin[,1]),]
        
        co <- try(harm4fit(tsin, f), TRUE)
        
        if (class(co) == "try-error"){
          
          SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f*tsin[,1]) + cos(2*pi*f*tsin[,1]) +
                        sin(4*pi*f*tsin[,1]) + cos(4*pi*f*tsin[,1]) + sin(6*pi*f*tsin[,1]) +
                        cos(6*pi*f*tsin[,1]) + sin(8*pi*f*tsin[,1]) + cos(8*pi*f*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
          
          co <- SSTlm$coefficients
          
          co[is.na(co)] <- 0
          
        }
        
        prevar <- var(ts[,2])
        
        ts[,2] <- ts[,2] - (co[2]*ts[,1] + co[3]*sin(2*pi*f*ts[,1]) + co[4]*cos(2*pi*f*ts[,1]) +
                              co[5]*sin(4*pi*f*ts[,1]) + co[6]*cos(4*pi*f*ts[,1]) + co[7]*sin(6*pi*f*ts[,1]) +
                              co[8]*cos(6*pi*f*ts[,1]) + co[9]*sin(8*pi*f*ts[,1]) + co[10]*cos(8*pi*f*ts[,1]))
        
        postvar <- var(ts[,2])
        
        varratio <- postvar / prevar
        
        features[62] <- varratio
        
        ts <- tsts
        
        if (quiet == FALSE){
          
          print("Plotting Light Curve folded onto #1 period...")
          
        }
        
        folded <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / p) + 1.25) %% 2) - 1
        
        folded2 <- ((((ts[,1] - ts[which.min(ts[,2]),1]) / p) + 0.25) %% 1)
        
        fts <- matrix(c(folded2, ts[,2]), nrow = length(folded2))
        
        fts <- fts[order(fts[,1]),]
        
        date_diff <- fts[2:length(fts[,1]),1] - fts[1:(length(fts[,1])-1),1]
        
        mag_diff <- fts[2:length(fts[,1]),2] - fts[1:(length(fts[,1])-1),2]
        
        date_diff <- date_diff[which(mag_diff != 0)]
        
        mag_diff <- mag_diff[which(mag_diff != 0)]
        
        slope <- date_diff / mag_diff
        
        features[63] <- quantile(slope, 0.1)
        
        features[64] <- quantile(slope, 0.9)
        
        pcs <- (1 / (n * stddev)) * cumsum(fts[,2] - meanmag)
        
        features[65] <- max(pcs) - min(pcs)
        
        varfold <- (sd(fts[,2])) ^ 2.0
        
        features[66] <-  (1.0 / ((n - 1) * varfold) * sum((diff(fts[,2])) ^ 2.0))
        
        folded3 <- folded2 - 1
        
        folded4 <- c(folded3, folded2)
        
        foldts <- c(ts[,2], ts[,2])
        
        if (quiet == FALSE){
          
          plot(folded, ts[,2], main = paste("", objects$usnoref[i], " Folded Light Curve spread from -1.0 to 1.0", sep=""), xlab = paste("Phase at the first period of ", p, " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, xlim = c(-1, 1), ylim = c(max(ts[,2]), min(ts[,2])), col = rgb(0, 0, 0, 0.25))
          
          if (ploterr == TRUE){
            
            magerr <- ts[,3]
            
            segments(folded, ts[,2]-magerr, folded, ts[,2]+magerr)
            
            epsilon = 0.01
            
            segments(folded-epsilon, ts[,2]-magerr, folded+epsilon, ts[,2]-magerr)
            
            segments(folded-epsilon, ts[,2]+magerr, folded+epsilon, ts[,2]+magerr)
            
          }
          
          plot(folded4, foldts, main = paste("", objects$usnoref[i], " Folded Light Curve duplicated from 0.0 - 1.0 to -1.0 to 1.0", sep=""), xlab = paste("Phase at the first period of ", p, " days", sep=""),
               ylab = "Apparent Magnitude", pch=19, xlim = c(-1, 1), ylim = c(max(foldts), min(foldts)), col = rgb(0, 0, 0, 0.25))
          
          if (ploterr == TRUE){
            
            magerr <- c(ts[,3], ts[,3])
            
            segments(folded4, foldts-magerr, folded4, foldts+magerr)
            
            epsilon = 0.01
            
            segments(folded4-epsilon, foldts-magerr, folded4+epsilon, foldts-magerr)
            
            segments(folded4-epsilon, foldts+magerr, folded4+epsilon, foldts+magerr)
            
          }
          
          print("Plotting Raw Light Curve with fitted model...")
          
          plot(ts[,1], ts[,2], main = paste("", objects$usnoref[i], " Raw Light Curve with fitted model", sep=""), xlab = "Modified Julian Date (Days)",
               ylab = "Apparent Magnitude", pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
          
          tspt <- sort(tsin[,1])
          
          tspr <- co[1] + co[2]*tspt + co[3]*sin(2*pi*f*tspt) + co[4]*cos(2*pi*f*tspt) +
            co[5]*sin(4*pi*f*tspt) + co[6]*cos(4*pi*f*tspt) + co[7]*sin(6*pi*f*tspt) +
            co[8]*cos(6*pi*f*tspt) + co[9]*sin(8*pi*f*tspt) + co[10]*cos(8*pi*f*tspt)
          
          lines(tspt, tspr, col = 2)
          
          #lines(sort(tsin[,1]), SSTlm$fitted.values, col = 2)
          
          if (ploterr == TRUE){
            
            magerr <- ts[,3]
            
            segments(ts[,1], ts[,2]-magerr, ts[,1], ts[,2]+magerr)
            
            epsilon2 = 4
            
            segments(ts[,1]-epsilon2, ts[,2]-magerr, ts[,1]+epsilon2, ts[,2]-magerr)
            
            segments(ts[,1]-epsilon2, ts[,2]+magerr, ts[,1]+epsilon2, ts[,2]+magerr)
            
          }
          
        }
        
        res <- ts[,2] - (co[1] + co[2]*ts[,1] + co[3]*sin(2*pi*f*ts[,1]) + co[4]*cos(2*pi*f*ts[,1]) +
                           co[5]*sin(4*pi*f*ts[,1]) + co[6]*cos(4*pi*f*ts[,1]) + co[7]*sin(6*pi*f*ts[,1]) +
                           co[8]*cos(6*pi*f*ts[,1]) + co[9]*sin(8*pi*f*ts[,1]) + co[10]*cos(8*pi*f*ts[,1]))
        
        keepres <- sigclip(res, sig = sigk, tol = 0.000001)
        
        res <- res[keepres]
        
        tres <- ts[keepres,1]
        
        if (quiet == FALSE){
          
          print("Plotting Residuals from fitted model...")
          
          plot(tres, res, main = paste("", objects$usnoref[i], " Residuals from fitted model", sep=""), xlab = "Modified Julian Date (Days)",
               ylab = "Residual Magnitude", pch=19, ylim = c(max(res), min(res)))
          
        }
        
        radnt <- ad.test(res)
        
        resnorm <- as.numeric(radnt$statistic)
        
        print(paste("Normality Test on Residuals (0.3 ~ normal, larger is non-normal): ", resnorm, sep=""))
        
        features[39] <- co[2]
        
        features[40] <- p
        
        features[41] <- p2
        
        features[42] <- pv
        
        features[43] <- pvd
        
        features[44] <- log10(fs)
        
        features[45] <- ns
        
        features[46] <- perspurlvl
        
        features[47] <- jit/amplitude
        
        R1 <- (co[3]^2.0 + co[4]^2.0)^0.5
        
        R2 <- (co[5]^2.0 + co[6]^2.0)^0.5
        
        R3 <- (co[7]^2.0 + co[8]^2.0)^0.5
        
        R4 <- (co[9]^2.0 + co[10]^2.0)^0.5
        
        features[48] <- R1
        
        features[49] <- R2
        
        features[50] <- R3
        
        features[51] <- R4
        
        features[52] <- R2/R1
        
        features[53] <- R3/R1
        
        features[54] <- R4/R1
        
        PH1 <- atan2(-co[3], co[4])
        
        PH2 <- atan2(-co[5], co[6]) - (2.0 * PH1)
        
        PH3 <- atan2(-co[7], co[8]) - (3.0 * PH1)
        
        PH4 <- atan2(-co[9], co[10]) - (4.0 * PH1)
        
        features[55] <- PH1
        
        features[56] <- PH2
        
        features[57] <- PH3
        
        features[58] <- PH4
        
        features[59] <- PH2 - PH1
        
        features[60] <- PH3 - PH1
        
        features[61] <- PH4 - PH1
        
        polyfit <- try(lcgbinpolyfit(target, spur = T, pero = p, delta = 1, seedno = seedno), TRUE)
        
        if (class(polyfit) == "try-error" | is.na(polyfit[38])){
          
          print("Polyfit failed to converge... Retrying #1...")
          
          polyfit <- try(lcgbinpolyfit(target, spur = T, pero = p, delta = 1, seedno = seedno+1), TRUE)
          
        }
        
        if (class(polyfit) == "try-error" | is.na(polyfit[38])){
          
          print("Polyfit failed to converge... Retrying #2...")
          
          polyfit <- try(lcgbinpolyfit(target, spur = T, pero = p, delta = 0.1, seedno = seedno+2), TRUE)
          
        }
        
        if (class(polyfit) == "try-error" | is.na(polyfit[38])){
          
          print("Polyfit failed to converge... Retrying #3...")
          
          polyfit <- try(lcgbinpolyfit(target, spur = T, pero = p, delta = 0.1, seedno = seedno+3), TRUE)
          
        }
        
        if (class(polyfit) == "try-error" | is.na(polyfit[38])){
          
          print("Polyfit failed to converge... Using a spline algorithm...")
          
          splines <- lcsplinefit(target, spur = T, pero = p)
          
          features[67] <- as.numeric(as.vector(splines$Binned.Ratio))
          
          features[68] <- as.numeric(as.vector(splines$Goodness.of.Fit))
          
          if (features[68] == Inf | is.nan(features[68])){
            
            features[68] <- 10
            
          }
          
          features[73] <- as.numeric(as.vector(splines$Amplitude))
          
          splines <- splines[7:106]
          
          splines <- as.numeric(as.vector(splines))
          
          sdmag <- sd(splines)
          
          features[69] <- sdmag
          
          n <- length(splines)
          
          meanmag <- mean(splines)
          
          features[70] <- skewness(splines)
          
          features[71] <- ((n-1)/((n-2)*(n-3)))*((n+1)*(((sum((splines - meanmag)^4.0)/n)/sum(((splines - meanmag)^2.0)/n)^2.0)-3)+6)
          
          features[72] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((splines - meanmag) / sdmag) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
          
          features[74] <- (length(splines[splines > (meanmag + sdmag)]) + length(splines[splines < (meanmag - sdmag)])) / length(splines)
          
          ics <- (1 / (n * sdmag)) * cumsum(splines - meanmag)
          
          features[75] <- max(ics) - min(ics)
          
          splines <- splines[2:length(splines)]
          
          splinepca <- t(as.matrix(splines)) %*% as.matrix(spca$rotation[,1:10])
          
          features[76:85] <- splinepca
          
        }
        else{
          
          features[67] <- as.numeric(as.vector(polyfit$Act.Bins)) / as.numeric(as.vector(polyfit$No.Bins))
          
          features[68] <- as.numeric(as.vector(polyfit$Chi.Squared))
          
          features[73] <- as.numeric(as.vector(polyfit$Int.Amplitude))
          
          normag <- polyfit$Norm.Magnitude
          
          polyfit <- polyfit[c(38:136)]
          
          polyfit <- as.numeric(as.vector(polyfit))
          
          n <- length(polyfit)
          
          meanmag <- mean(polyfit)
          
          sdmag <- sd(polyfit)
          
          features[69] <- sdmag
          
          features[70] <- skewness(polyfit)
          
          features[71] <- ((n-1)/((n-2)*(n-3)))*((n+1)*(((sum((polyfit - meanmag)^4.0)/n)/sum(((polyfit - meanmag)^2.0)/n)^2.0)-3)+6)
          
          features[72] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((polyfit - meanmag) / sdmag) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
          
          features[74] <- (length(polyfit[polyfit > (meanmag + sdmag)]) + length(polyfit[polyfit < (meanmag - sdmag)])) / length(polyfit)
          
          ics <- (1 / (n * sdmag)) * cumsum(polyfit - meanmag)
          
          features[75] <- max(ics) - min(ics)
          
          #polyfit <- polyfit[c(1:49, 51:length(polyfit))]
          
          polyfitpca <- t(as.matrix(polyfit)) %*% as.matrix(pca$rotation[,1:10])
          
          pcarecon <- as.vector(as.matrix(pca$rotation[,1:10]) %*% t(as.matrix(polyfitpca)))
          
          lines(seq(-0.49, 0.49, by = 0.01), ((pcarecon*features[73]*2.0) + normag), lwd = 2, col = "blue")
          
          features[76:85] <- polyfitpca
          
          ## MAKE ECLIPSE TRACKING FEATURES USING THE MINIMUM (MAX) OF EACH OF THE PHASE REGIONS @ 2P.
          
          ## ALSO PERHAPS PREVENTING MINIMA IN ADJACENT BINS TO FIX THE KNOTS AT MINIMA PROBLEM.
          
        }
        
        polyfit2 <- try(lcgbinpolyfit(target, spur = T, pero = p2, delta = 1, seedno = seedno), TRUE)
        
        if (class(polyfit2) == "try-error" | is.na(polyfit2[38])){
          
          print("Polyfit failed to converge... Retrying #1...")
          
          polyfit2 <- try(lcgbinpolyfit(target, spur = T, pero = p2, delta = 1, seedno = seedno+1), TRUE)
          
        }
        
        if (class(polyfit2) == "try-error" | is.na(polyfit2[38])){
          
          print("Polyfit failed to converge... Retrying #2...")
          
          polyfit2 <- try(lcgbinpolyfit(target, spur = T, pero = p2, delta = 0.1, seedno = seedno+2), TRUE)
          
        }
        
        if (class(polyfit2) == "try-error" | is.na(polyfit2[38])){
          
          print("Polyfit failed to converge... Retrying #3...")
          
          polyfit2 <- try(lcgbinpolyfit(target, spur = T, pero = p2, delta = 0.1, seedno = seedno+3), TRUE)
          
        }
        
        if (class(polyfit2) == "try-error" | is.na(polyfit2[38])){
          
          print("Polyfit failed to converge... Using a spline algorithm...")
          
          splines2 <- lcsplinefit(target, spur = T, pero = p2)
          
          features[86] <- as.numeric(as.vector(splines2$Binned.Ratio))
          
          features[87] <- as.numeric(as.vector(splines2$Goodness.of.Fit))
          
          if (features[87] == Inf | is.nan(features[87])){
            
            features[87] <- 10
            
          }
          
          features[92] <- as.numeric(as.vector(splines2$Amplitude))
          
          splines2 <- splines2[7:106]
          
          splines2 <- as.numeric(as.vector(splines2))
          
          sdmag <- sd(splines2)
          
          features[88] <- sdmag
          
          n <- length(splines2)
          
          meanmag <- mean(splines2)
          
          features[89] <- skewness(splines2)
          
          features[90] <- ((n-1)/((n-2)*(n-3)))*((n+1)*(((sum((splines2 - meanmag)^4.0)/n)/sum(((splines2 - meanmag)^2.0)/n)^2.0)-3)+6)
          
          features[91] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((splines2 - meanmag) / sdmag) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
          
          features[93] <- (length(splines2[splines2 > (meanmag + sdmag)]) + length(splines2[splines2 < (meanmag - sdmag)])) / length(splines2)
          
          ics <- (1 / (n * sdmag)) * cumsum(splines2 - meanmag)
          
          features[94] <- max(ics) - min(ics)
          
          splines2 <- splines2[2:length(splines2)]
          
          spline2pca <- t(as.matrix(splines2)) %*% as.matrix(spca$rotation[,1:10])
          
          features[95:104] <- spline2pca
          
          intphase <- seq(-0.49, 0.49, by = 0.01) + 0.5
          
          locmax <- localMaxima(splines2)
          
          locmin <- localMinima(splines2)
          
          locmax <- locmax[which(!(locmax == 99))]
          
          locmin <- locmin[which(!(locmin == 99))]
          
          if (length(locmax) > 2){
            
            locmax <- locmax[order(splines2[locmax], decreasing = T)][1:2]
            
          }
          
          mulloc <- FALSE
          
          if (length(locmin) > 2){
            
            #locmin <- locmin[order(polyfit2[locmin], decreasing = F)][1:2]
            
            locmin <- locmin[which(!(locmin == 1 | locmin == 99))]
            
            if (length(locmin) > 2){
              
              locmin1 <- locmin[which(locmin > min(locmax) & locmin < max(locmax))]
              
              locmin2 <- locmin[which(!(locmin > min(locmax) & locmin < max(locmax)))]
              
              mulloc <- TRUE
              
              print(locmin2)
              
            }
            
          }
          
          splinemax <- cbind(intphase[locmax], splines2[locmax])
          
          splinemin <- cbind(intphase[locmin], splines2[locmin])
          
          if (mulloc == TRUE){
            
            locset1 <- splinemin[which(locmin %in% locmin1),]
            
            locset2 <- splinemin[which(locmin %in% locmin2),]
            
            if (!is.null(nrow(locset1))){
              
              locset1 <- colMeans(locset1, na.rm = T)
              
            }
            
            if (!is.null(nrow(locset2))){
              
              locset2 <- colMeans(locset2, na.rm = T)
              
            }
            
            splinemin <- rbind(locset1, locset2)
            
            print(splinemin)
            
          }
          
          splinepos <- rbind(splinemax, splinemin)
          
          splinepos <- splinepos[order((splinepos[,1] - splinepos[which.max(splinepos[,2]),1]) %% 1),]
          
          eclfeat1 <- try(abs(splinepos[1,2] - splinepos[3,2])*features[92]*2.0, TRUE)
          
          if (class(eclfeat1) == "try-error"){
            
            features[105] <- 0
            
          }
          else{
            
            features[105] <- eclfeat1
            
          }
          
          eclfeat2 <- try(abs(splinepos[2,2] - splinepos[4,2])*features[92]*2.0, TRUE)
          
          if (class(eclfeat2) == "try-error"){
            
            features[106] <- 0
            
          }
          else{
            
            features[106] <- eclfeat2
            
          }
          
          eclfeat3 <- try(3.0 * ((splinepos[1,1] - (splinepos[4,1] - 1)) / (splinepos[3,1] - (splinepos[4,1] - 1))), TRUE)
          
          if (class(eclfeat3) == "try-error"){
            
            features[107] <- 0
            
          }
          else{
            
            features[107] <- eclfeat3
            
          }
          
          eclfeat4 <- try(splinepos[4,1] - 1, TRUE)
          
          if (class(eclfeat4) == "try-error"){
            
            features[108] <- 0
            
          }
          else{
            
            features[108] <- eclfeat4
            
          }
          
        }
        else{
          
          features[86] <- as.numeric(as.vector(polyfit2$Act.Bins)) / as.numeric(as.vector(polyfit2$No.Bins))
          
          features[87] <- as.numeric(as.vector(polyfit2$Chi.Squared))
          
          features[92] <- as.numeric(as.vector(polyfit2$Int.Amplitude))
          
          normag <- polyfit2$Norm.Magnitude
          
          knotphase <- c(polyfit2$Knot.1, polyfit2$Knot.2, polyfit2$Knot.3, polyfit2$Knot.4)
          
          polyfit2 <- polyfit2[c(38:136)]
          
          polyfit2 <- as.numeric(as.vector(polyfit2))
          
          n <- length(polyfit2)
          
          meanmag <- mean(polyfit2)
          
          sdmag <- sd(polyfit2)
          
          features[88] <- sdmag
          
          features[89] <- skewness(polyfit2)
          
          features[90] <- ((n-1)/((n-2)*(n-3)))*((n+1)*(((sum((polyfit2 - meanmag)^4.0)/n)/sum(((polyfit2 - meanmag)^2.0)/n)^2.0)-3)+6)
          
          features[91] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((polyfit2 - meanmag) / sdmag) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
          
          features[93] <- (length(polyfit2[polyfit2 > (meanmag + sdmag)]) + length(polyfit2[polyfit2 < (meanmag - sdmag)])) / length(polyfit2)
          
          ics <- (1 / (n * sdmag)) * cumsum(polyfit2 - meanmag)
          
          features[94] <- max(ics) - min(ics)
          
          #polyfit2 <- polyfit2[c(1:49, 51:length(polyfit2))]
          
          polyfit2pca <- t(as.matrix(polyfit2)) %*% as.matrix(pca$rotation[,1:10])
          
          pcarecon2 <- as.vector(as.matrix(pca$rotation[,1:10]) %*% t(as.matrix(polyfit2pca)))
          
          lines(seq(-0.49, 0.49, by = 0.01), ((pcarecon2*features[92]*2.0) + normag), lwd = 2, col = "blue")
          
          features[95:104] <- polyfit2pca
          
          ## MAKE ECLIPSE TRACKING FEATURES USING THE MINIMUM (MAX) OF EACH OF THE PHASE REGIONS @ 2P.
          
          ## ALSO PERHAPS PREVENTING MINIMA IN ADJACENT BINS TO FIX THE KNOTS AT MINIMA PROBLEM.
          
          intphase <- seq(-0.49, 0.49, by = 0.01)
          
          locmax <- localMaxima(polyfit2)
          
          locmin <- localMinima(polyfit2)
          
          print(locmax)
          
          print(locmin)
          
          if (sum(locmax == 1 | locmax == 99) > 1){
            
            locmax <- locmax[which(!(locmax == 99))]
            
            locmin <- locmin[which(!(locmin == 99))]
            
          }
          
          if (length(locmax) > 2){
            
            locmax <- locmax[order(polyfit2[locmax], decreasing = T)][1:2]
            
          }
          
          mulloc <- FALSE
          
          if (length(locmin) > 2){
            
            #locmin <- locmin[order(polyfit2[locmin], decreasing = F)][1:2]
            
            locmin <- locmin[which(!(locmin == 1 | locmin == 99))]
            
            if (length(locmin) > 2){
              
              locmin1 <- locmin[which(locmin > min(locmax) & locmin < max(locmax))]
              
              locmin2 <- locmin[which(!(locmin > min(locmax) & locmin < max(locmax)))]
              
              mulloc <- TRUE
              
              print(locmin2)
              
            }
            
          }
          
          print(locmax)
          
          print(locmin)
          
          polymax <- cbind(intphase[locmax], polyfit2[locmax])
          
          polymin <- cbind(intphase[locmin], polyfit2[locmin])
          
          if (mulloc == TRUE){
            
            locset1 <- polymin[which(locmin %in% locmin1),]
            
            locset2 <- polymin[which(locmin %in% locmin2),]
            
            if (!is.null(nrow(locset1))){
              
              locset1 <- colMeans(locset1, na.rm = T)
              
            }
            
            if (!is.null(nrow(locset2))){
              
              locset2 <- colMeans(locset2, na.rm = T)
              
            }
            
            polymin <- rbind(locset1, locset2)
            
            print(polymin)
            
          }
          
          #polymod <- cbind(intphase, polyfit2)
          
          #polyknot1 <- polymod[which(polymod[,1] >= knotphase[1] & polymod[,1] < knotphase[2]),]
          
          #polyknot2 <- polymod[which(polymod[,1] >= knotphase[2] & polymod[,1] < knotphase[3]),]
          
          #polyknot3 <- polymod[which(polymod[,1] >= knotphase[3] & polymod[,1] < knotphase[4]),]
          
          #polyknot4a <- polymod[which(polymod[,1] < knotphase[1]),]
          
          #polyknot4b <- polymod[which(polymod[,1] >= knotphase[4]),]
          
          #polyknot4 <- rbind(polyknot4a, polyknot4b)
          
          #polymax <- rbind(polyknot1[which.max(polyknot1[,2]),], polyknot2[which.max(polyknot2[,2]),], polyknot3[which.max(polyknot3[,2]),], polyknot4[which.max(polyknot4[,2]),])
          
          #polymin <- rbind(polyknot1[which.min(polyknot1[,2]),], polyknot2[which.min(polyknot2[,2]),], polyknot3[which.min(polyknot3[,2]),], polyknot4[which.min(polyknot4[,2]),])
          
          #priloc <- which.max(polymax[,2])
          
          #polypos <- rbind(polymax[((priloc-1) %% 4)+1,], polymin[((priloc) %% 4)+1,], polymax[((priloc+1) %% 4)+1,], polymin[((priloc+2) %% 4)+1,])
          
          polypos <- rbind(polymax, polymin)
          
          polypos <- polypos[order((polypos[,1] - polypos[which.max(polypos[,2]),1]) %% 1),]
          
          print(polypos)
          
          eclfeat1 <- try(abs(polypos[1,2] - polypos[3,2])*features[92]*2.0, TRUE)
          
          if (class(eclfeat1) == "try-error"){
            
            features[105] <- 0
            
          }
          else{
            
            features[105] <- eclfeat1
            
          }
          
          eclfeat2 <- try(abs(polypos[2,2] - polypos[4,2])*features[92]*2.0, TRUE)
          
          if (class(eclfeat2) == "try-error"){
            
            features[106] <- 0
            
          }
          else{
            
            features[106] <- eclfeat2
            
          }
          
          if (polypos[3,1] >= 0){
            
            eclfeat3 <- try(3.0 * ((polypos[1,1] - polypos[4,1]) / (polypos[3,1] - polypos[4,1])), TRUE)
            
          }
          else{
            
            eclfeat3 <- try(3.0 * ((polypos[1,1] - polypos[4,1]) / (polypos[3,1]%%1 - polypos[4,1])), TRUE)
            
          }
          
          if (class(eclfeat3) == "try-error"){
            
            features[107] <- 0
            
          }
          else{
            
            features[107] <- eclfeat3
            
          }
          
          eclfeat4 <- try(polypos[4,1], TRUE)
          
          if (class(eclfeat4) == "try-error"){
            
            features[108] <- 0
            
          }
          else{
            
            features[108] <- eclfeat4
            
          }
          
        }
        
        features[109] <- resnorm
        
        madres <- median(abs(res - median(res)))
        
        features[110] <- madres / features[24]
        
        features[111] <- sum(diff(ts[,2])^2.0) / var(ts[,2])
        
        features[112] <- perratio
        
        featvec[i,] <- features
        
      }
      
    }
    
    if (nonper == FALSE){
      
      featvec <- as.data.frame(featvec)
      
      featvec <- cbind(objects$usnoref, objects$RA, objects$DEClin, featvec)
      
      colnames(featvec, do.NULL = FALSE)
      
      colnames(featvec) <- c("USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                             "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                             "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                             "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                             "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value", "Slope of Linear Trend",
                             "Period", "Double Period", "Period Confidence", "Daily Confidence", "Bestfit Stat", "Num Spurious", "Period Spurious Level", "Jitter Multiplier", "H1", "H2", "H3", "H4", "R21", "R31", "R41", "PH1", "PH2", "PH3", "PH4", "PH21", "PH31", "PH41", "Variance Ratio", "mp10", "mp90", "Psi-cs",
                             "Psi-Eta", "P Binned Ratio", "P Goodness of Fit", "P Int StD", "P Int Skewness", "P Int Kurtosis", "P Int Small Kurtosis", "P Int Amplitude",
                             "P Int Beyond 1 StD", "P Int cs", paste("P PCA", 1:10, sep=""), "P2 Binned Ratio", "P2 Goodness of Fit", "P2 Int StD", "P2 Int Skewness", "P2 Int Kurtosis", "P2 Int Small Kurtosis", "P2 Int Amplitude",
                             "P2 Int Beyond 1 StD", "P2 Int cs", paste("P2 PCA", 1:10, sep=""), "Eclipse Max Delta Mags", "Eclipse Min Delta Mags", "Eclipse Phase Ratio", "Reference Phase", "Residual Normality", "Residual Raw Scatter", "Squared Diffs over Variance", "Period Double Ratio")
      
      colnames(featvec) <- make.names(colnames(featvec), unique = T)
      
    }
    
    else{
      
      featvec <- as.data.frame(featvec)
      
      featvec <- featvec[,c(1:38)]
      
      featvec <- cbind(objects$usnoref, objects$RA, objects$DEClin, featvec)
      
      colnames(featvec, do.NULL = FALSE)
      
      colnames(featvec) <- c("USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                             "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                             "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                             "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                             "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value")
      
      colnames(featvec) <- make.names(colnames(featvec), unique = T)
      
    }
    
    close(channel)
    
    print(Sys.time() - start)
    
    gc()
    
    if (geo == TRUE){
      
      featvec <- featvec[1,]
      
    }
    
    featvec
    
  }
  
}