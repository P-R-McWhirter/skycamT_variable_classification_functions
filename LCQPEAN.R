LCQPEAN <- function(ra, declin, radii = 0.01, lim = 2, fit = "sin", redo = 5, seedno = c(20, 30, 40, 50, 60), pop = 20, ps = 20, nogen = 100, crov = 0.65, mut = 0.03, fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20, alldb = FALSE, geo = FALSE, quiet = FALSE, fulltest = TRUE) {
  
  closeAllConnections()
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  library(nonnest2)
  
  start <- Sys.time()
  
  if (quiet == FALSE){
    
    print("Executing Software...")
    
  }
  
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
    
    featvec <- matrix(0, length(objects$entries), 5)
    
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
      
      prelen <- length(ts[,2])
      
      ts <- unique(ts)
      
      ts <- ts[order(ts[,1]),]
      
      postlen <- length(ts[,2])
      
      difflen <- abs(prelen - postlen)
      
      if (quiet == FALSE){
        
        print(paste("There were ", difflen, " duplicate observations found within the data.", sep = ""))
        
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
      
      features <- rep(0, 5)
      
      features[1] <- objects$entries[i]
      
      meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
      
      amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.01))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.01)))) / 2)
      
      tsts <- ts
      
      per <- rep(0, redo)
      
      fstat <- rep(0, redo)
      
      models <- list()
      
      #print("Computing the Spurious Periods using a Lomb-Scargle Periodogram...")
      
      #sper <- spurperlsp(ts, from = freqmin, to = 4, ofac = 10, name = objects$usnoref[i], plot = !quiet)
      
      #spurs <- sper$scanned[order(sper$power, decreasing = TRUE)]
      
      #spurpwr <- sper$power[order(sper$power, decreasing = TRUE)]
      
      #spurs <- spurs[which(spurpwr > (max(spurpwr)*0.8))]
      
      #spurs <- spurs[1:30]
      
      #spurs <- c((1 / 29.53), spurs, 2*spurs)
      
      spurs <- c(1/0.9972, 1, 1/1.00274, 1/365.24)
      
      spurs <- c(spurs*2, spurs, spurs/2, 1/29.53)
      
      for (k in 1:redo){
        
        print(paste("Initiating Genetic Algorithm...       [", k, "/", redo, "]", sep = ""))
        
        ans <- LCGEN(ts, spurs = spurs, fit = fit, melc = meanmag, amplc = amplitude, from = freqmin, to = freqmax, seedno = seedno[k], pop = pop, pairups = ps, nogen = nogen, crossover = crov, mutation = mut, fdif = 0.6, dfrac = 0.7, clipcut = clipcut, bins = bins)
        
        per[k] <- as.numeric(ans$period)
        
        fstat[k] <- as.numeric(ans$fstat)
        
        models[[k]] <- ans$model
        
        loctest <- rep(0, 8)
        
        #for (z in 1:8){
        
        #loctest[z] <- vuongtest(ans$model, ans$othmodel[[z]])$LRTstat
        
        #}
        
        #if (min(loctest) < 0){
        
        #print(paste("Replacing suspected local minimum: ", per[k], ", with better candidate period: ", ans$othper[[which.min(loctest)]], ".", sep=""))
        
        #per[k] <- as.numeric(ans$othper[[which.min(loctest)]])
        
        #models[[k]] <- ans$othmodel[[which.min(loctest)]]
        
        #}
        
      }
      
      print(paste("Completed all ", redo, " runs...", sep = ""))
      
      print(paste("Chi-squared minimized Period found: ", as.character(per[which.min(fstat)]), ".", sep=""))
      
      LRT <- matrix(0, nrow = redo, ncol = redo)
      
      conf <- rep(0, redo)
      
      if (fulltest == TRUE){
        
        for (k in 1:redo){
          
          for (z in 1:redo){
            
            vuong <- vuongtest(models[[k]], models[[z]])
            
            LRT[k,z] <- vuong$LRTstat
            
            if (is.na(LRT[k,z])){
              
              LRT[k,z] <- 0
              
            }
            
          }
          
          conf[k] <- sum(LRT[k,])
          
        }
        
      }
      else{
        
        for (k in 1:redo){
          
          conf[k] <- vuongtest(models[[k]], ans$onemodel)$LRTstat
          
        }
        
      }
      
      vuongconf <- as.numeric(vuongtest(models[[which.max(conf)]], ans$onemodel)$LRTstat)
      
      print(paste("Vuong Closeness Test maximized Period found: ", as.character(per[which.max(conf)]), ".", sep=""))
      
      print(paste("Vuong Closeness Test LRT for this period against a constant model: ", as.character(vuongconf), sep=""))
      
      features[2] <- as.character(per[1])
      
      if (redo > 1){
        
        for (k in 2:redo){
          
          features[2] <- paste(features[2], as.character(per[k]), sep = " ")
          
        }
        
      }
      
      features[3] <- per[which.min(fstat)]
      
      features[4] <- per[which.max(conf)]
      
      features[5] <- vuongconf
      
      featvec[i,] <- features
      
    }
    
    featvec <- as.data.frame(featvec)
    
    featvec <- cbind(objects$usnoref, featvec)
    
    colnames(featvec, do.NULL = FALSE)
    
    colnames(featvec) <- c("USnoref", "#.Observations", "Output.Periods", "Genetic.Period", "Vuong.Period", "Vuong.Confidence")
    
    close(channel)
    
    print(Sys.time() - start)
    
    gc()
    
    if (geo == TRUE){
      
      featvec <- featvec[1,]
      
    }
    
    featvec
    
  }
  
}