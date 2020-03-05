LCPEAN <- function(ra, declin, radii = 0.01, lim = 2, fit = "lsp", jit = 0, noimod = 1, sigk = 2.5, seascoe, trendrem = TRUE, rangcut = rangcut, redo = 5, seedno = 20, pop = 20, ps = 20, nogen = 100, crov = 0.65, mut = 0.03, fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20, alldb = FALSE, geo = FALSE, quiet = FALSE, fulltest = TRUE, tenper = FALSE) {
  
  closeAllConnections()
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  library(nonnest2)
  library(DescTools)
  
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
    
    if (quiet == FALSE){
      
      if (alldb == TRUE){
        
        print(paste("Reading in all ", length(objects$usnoref), " object(s) with greater than ", lim, " observations.", sep=""))
        
      }
      
      else{
        
        print(paste("A total of ", length(objects$usnoref),
                    " object(s) identified within the search radius of ", radii, " degree(s) with greater than ", lim, " observations.", sep=""))
        
      }
      
    }
    
    featvec <- matrix(0, length(objects$entries), 9)
    
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
      
      # <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
      
      #ts <- unique(ts)
      
      #preav <- weighted.mean(ts[,2], (1 / ts[,3]))
      
      #ts[,2] <- ts[,2] - preav
      
      #tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
      
      #tsin <- rbind(tsin, ts)
      
      #tsin <- tsin[order(tsin[,1]),]
      
      #seascoe <- try(harm2fit(tsin, (1/365.24217), lambda = 1e-2), TRUE)
      
      #if (class(seascoe) == "try-error"){
      
      #seastrend <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/365.24217) + cos(2*pi*tsin[,1]/365.24217) +
      #                  sin(4*pi*tsin[,1]/365.24217) + cos(4*pi*tsin[,1]/365.24217), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      #seascoe <- seastrend$coefficients
      
      #seascoe[is.na(seascoe)] <- 0
      
      #}
      #else{
      
      #seastrend <- list()
      
      #seastrend$fitted.values <- seascoe[1] + seascoe[2]*tsin[,1] + seascoe[3]*sin(2*pi*tsin[,1]/365.24217) + seascoe[4]*cos(2*pi*tsin[,1]/365.24217) +
      #                            seascoe[5]*sin(4*pi*tsin[,1]/365.24217) + seascoe[6]*cos(4*pi*tsin[,1]/365.24217)
      
      #}
      
      #preav <- weighted.mean(ts[,2], (1 / ts[,3]))
      
      #tsin[,2] <- tsin[,2] - (seascoe[1] + seascoe[2]*tsin[,1] + seascoe[3]*sin(2*pi*tsin[,1]/365.24217) + seascoe[4]*cos(2*pi*tsin[,1]/365.24217) +
      #                      seascoe[5]*sin(4*pi*tsin[,1]/365.24217) + seascoe[6]*cos(4*pi*tsin[,1]/365.24217) - preav)
      
      #seastrend2 <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]) + cos(2*pi*tsin[,1]) +
      #                  sin(4*pi*tsin[,1]) + cos(4*pi*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      #seascoe2 <- seastrend2$coefficients
      
      #seascoe2[is.na(seascoe2)] <- 0
      
      #tsin[,2] <- tsin[,2] - (seascoe2[1] + seascoe2[2]*tsin[,1] + seascoe2[3]*sin(2*pi*tsin[,1]) + seascoe2[4]*cos(2*pi*tsin[,1]) +
      #                      seascoe2[5]*sin(4*pi*tsin[,1]) + seascoe2[6]*cos(4*pi*tsin[,1]) - preav)
      
      #maxspur <- max(ts[,1]) - min(ts[,1])
      
      #seastrend3 <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/maxspur) + cos(2*pi*tsin[,1]/maxspur) +
      #                  sin(4*pi*tsin[,1]/maxspur) + cos(4*pi*tsin[,1]/maxspur), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      #seascoe3 <- seastrend3$coefficients
      
      #seascoe3[is.na(seascoe3)] <- 0
      
      #seasamp <- sqrt(as.numeric(seascoe[3])^2.0 + as.numeric(seascoe[4])^2.0)
      
      #seasamp2 <- sqrt(as.numeric(seascoe[5])^2.0 + as.numeric(seascoe[6])^2.0)
      
      #set.seed(100)
      
      #info <- info[sample(1:nrow(info), ceiling(nrow(info)*0.88)),]
      
      #info <- info[order(info$Rcaterr),]
      
      #info <- info[1:ceiling(nrow(info)*0.9),]
      
      #ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
      
      #ts <- unique(ts)
      
      #ts[,2] <- ts[,2] - preav
      
      #maxtime <- min(c(500, max(ts[,1]) - min(ts[,1])))
      
      #if (maxtime > 80){
      
      #LPVsrh <- BGLS2(ts, plow = 80, phigh = maxtime, ofac = 5, dt = NULL, lent = NULL, spur = F, jit = jit*4, plot = F)
      
      #LPVguess <- 1 / LPVsrh$f[which.max(LPVsrh$p)]
      
      #print(LPVguess)
      
      #LPVfit <- try(harm2fit(tsin, (1/LPVguess), lambda = 1e-2), TRUE)
      
      #if (class(LPVfit) == "try-error"){
      
      #LPVfitted <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]/LPVguess) + cos(2*pi*tsin[,1]/LPVguess) +
      #                  sin(4*pi*tsin[,1]/LPVguess) + cos(4*pi*tsin[,1]/LPVguess), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
      #LPVfit <- LPVfitted$coefficients
      
      #LPVfit[is.na(LPVfit)] <- 0
      
      #}
      #else{
      
      #LPVfitted <- list()
      
      #LPVfitted$fitted.values <- LPVfit[1] + LPVfit[2]*tsin[,1] + LPVfit[3]*sin(2*pi*tsin[,1]/LPVguess) + LPVfit[4]*cos(2*pi*tsin[,1]/LPVguess) +
      #                            LPVfit[5]*sin(4*pi*tsin[,1]/LPVguess) + LPVfit[6]*cos(4*pi*tsin[,1]/LPVguess)
      
      #}
      
      #LPVamp <- sqrt(as.numeric(LPVfit[3])^2.0 + as.numeric(LPVfit[4])^2.0)
      
      #LPVamp2 <- sqrt(as.numeric(LPVfit[5])^2.0 + as.numeric(LPVfit[6])^2.0)
      
      #print(max(LPVfitted$fitted.values) - max(seastrend$fitted.values))
      
      #print(min(LPVfitted$fitted.values) - min(seastrend$fitted.values))
      
      #print(c(seasamp, seasamp2, LPVamp, LPVamp2))
      
      #ampsim <- as.vector(abs(outer(c(seasamp, seasamp2), c(LPVamp, LPVamp2), "-")))
      
      #print(ampsim)
      
      #ampratio <- seasamp2 / seasamp
      
      #ampratio2 <- LPVamp2 / LPVamp
      
      #print(c(ampratio, ampratio2))
      
      #ampcut <- ampratio / ampratio2
      
      #print(ampcut)
      
      #print(any(ampsim < 0.05))
      
      #print(c(seasamp, seasamp2))
      
      #trendcor <- cor(seastrend$fitted.values, LPVfitted$fitted.values)
      
      #print(trendcor)
      
      #}
      
      keep <- sigclip(info$Rcat, sig = sigk, tol = 0.000001)
      
      info <- info[keep,]
      
      ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
      
      #ts <- unique(ts)
      
      #plot(ts[,1:2], pch=19, ylim = c(max(ts[,2]), min(ts[,2])))
      
      #lines(tsin[,1], seastrend$fitted.values+preav, col = "red")
      
      #lines(tsin[,1], LPVfitted$fitted.values+preav, col = "blue")
      
      #plot(as.numeric(as.vector(info$MJD)), as.numeric(as.vector(info$Rcat)), pch = 19, ylim=c(max(as.numeric(as.vector(info$Rcat))), min(as.numeric(as.vector(info$Rcat)))))
      
      if (quiet == FALSE){
        
        print(paste("Object with reference id '", objects$usnoref[i],
                    "' successfully read in. Object has ", objects$entries[i], " observation(s).", sep=""))
        
      }
      
      #if (maxtime > 80){
      
      #if (any(ampsim < LPVcut[1]) & trendcor > LPVcut[2] & seasamp < LPVcut[3]){
      
      #print(paste("Amplitude Ratio between LPV and one year trend is : ", any(ampsim < 0.1), ". Prewhitening Light Curve with one year periodic model.", sep = ""))
      
      #ts[,2] <- ts[,2] - (seascoe[1] + seascoe[2]*ts[,1] + seascoe[3]*sin(2*pi*ts[,1]/365.24217) + seascoe[4]*cos(2*pi*ts[,1]/365.24217) +
      #                    seascoe[5]*sin(4*pi*ts[,1]/365.24217) + seascoe[6]*cos(4*pi*ts[,1]/365.24217) - preav)
      
      #ts[,2] <- ts[,2] - (seascoe2[1] + seascoe2[2]*ts[,1] + seascoe2[3]*sin(2*pi*ts[,1]) + seascoe2[4]*cos(2*pi*ts[,1]) +
      #                   seascoe2[5]*sin(4*pi*ts[,1]) + seascoe2[6]*cos(4*pi*ts[,1]) - preav)
      
      #ts[,2] <- ts[,2] - (seascoe3[1] + seascoe3[2]*ts[,1] + seascoe3[3]*sin(2*pi*ts[,1]/maxspur) + seascoe3[4]*cos(2*pi*ts[,1]/maxspur) +
      #                   seascoe3[5]*sin(4*pi*ts[,1]/maxspur) + seascoe3[6]*cos(4*pi*ts[,1]/maxspur) - preav)
      
      #}
      #else{
      
      #print(paste("Amplitude Ratio between LPV and one year trend is : ", any(ampsim < 0.05), ". Possible long period signal/alias, skipping yearly prewhitening.", sep = ""))
      
      #}
      
      #}
      #else{
      
      #print("Light Curve baseline too short for yearly trend removal.")
      
      #}
      
      #plot(ts[,1:2], pch = 19, ylim = c(max(ts[,2]), min(ts[,2])))
      
      prelen <- length(ts[,2])
      
      ts <- unique(ts)
      
      #ts[which(ts[,1] < 55200),2] <- ts[which(ts[,1] < 55200),2] - 0.1
      
      preav <- weighted.mean(ts[,2], (1 / ts[,3]))
      
      if (trendrem == TRUE){
        
        ts[,2] <- ts[,2] - (seascoe[1] + seascoe[2]*ts[,1] + seascoe[3]*sin(2*pi*ts[,1]/365.24217) + seascoe[4]*cos(2*pi*ts[,1]/365.24217) +
                              seascoe[5]*sin(4*pi*ts[,1]/365.24217) + seascoe[6]*cos(4*pi*ts[,1]/365.24217) - preav)
        
      }
      
      lcrpart <- try(lcrpartfit(ts, quiet = TRUE, cut = rangcut), TRUE)
      
      if (class(lcrpart) != "try-error"){
        
        ts <- lcrpart$ts
        
      }
      
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
      
      features <- rep(0, 9)
      
      features[1] <- length(ts[,1])
      
      meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
      
      #amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.4))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.4)))) / 2)
      
      amplitude <- abs(max(ts[,2]) - min(ts[,2]))/2.0
      
      tsts <- ts
      
      per <- rep(0, redo)
      
      fstat <- rep(0, redo)
      
      models <- list()
      
      print("Computing the Spurious Periods using a Lomb-Scargle Periodogram...")
      
      sper <- spurperlsp(ts, from = freqmin, to = 1, ofac = 5, name = objects$usnoref[i], plot = !quiet)
      
      spurs <- sper$scanned[order(sper$power, decreasing = TRUE)]
      
      spurpwr <- sper$power[order(sper$power, decreasing = TRUE)]
      
      spurs <- spurs[which(spurpwr > (max(spurpwr)*0.5))]
      
      numspurs <- length(spurs)
      
      #spurs <- c((1 / 29.53), spurs, 2*spurs)
      
      spurs <- c(1/0.9972, 1, 1/1.00274, spurs)
      
      spurs <- c(spurs*3, spurs*2, spurs, spurs/2, spurs/3, 1/29.53, abs((29.53^-1 + 1)^-1), abs((29.53^-1 - 1)^-1))
      
      spurs <- sort(spurs)
      
      print(paste("Number of detected Spurious Periods: ", numspurs, ".", sep=""))
      
      print("Initiating Period clustering Genetic Algorithm...")
      
      ansset <- LCMGEN(ts, fit = fit, melc = meanmag, amplc = amplitude, from = freqmin, to = freqmax, jit = jit, noimod = noimod, seedno = seedno, redo = redo, pop = pop, pairups = ps, nogen = nogen, crossover = crov, mutation = mut, fdif = 0.6, dfrac = 0.7, clipcut = clipcut, bins = bins)
      
      ansset2 <- LCMGEN(ts, fit = fit, melc = meanmag, amplc = amplitude, from = freqmin, to = freqmax, jit = jit, noimod = noimod, seedno = seedno+100, redo = redo, pop = pop, pairups = ps, nogen = nogen, crossover = crov, mutation = mut, fdif = 0.6, dfrac = 0.7, clipcut = clipcut, bins = bins)
      
      prlp <- rep(FALSE, redo)
      
      mnlp <- matrix(FALSE, nrow = redo, ncol = redo)
      
      wolp <- rep(0, redo)
      
      for (k in 1:redo){
        
        for (l in 1:redo){
          
          if (k > l){
            
            mnlp[k,l] <- c(ansset$permin[k], ansset$permax[k]) %overlaps% c(ansset$permin[l], ansset$permax[l])
            
          }
          
        }
        
        if (any(mnlp[k,] == TRUE)){
          
          wolp[k] <- which(mnlp[k,] == TRUE)[1]
          
        }
        else{
          
          wolp[k] <- NA
          
        }
        
      }
      
      for (k in 1:redo){
        
        ovlp <- rep(FALSE, redo)
        
        for (l in 1:redo){
          
          ovlp[l] <- c(ansset2$permin[k], ansset2$permax[k]) %overlaps% c(ansset$permin[l], ansset$permax[l])
          
        }
        
        prlp[k] <- any(ovlp)
        
      }
      
      nanum <- length(which(is.na(wolp)))
      
      if (any(prlp == FALSE) & nanum > 0){
        
        print("Replacing duplicate candidate period(s) with unknown candidate(s) from second optimisation...")
        
        newval <- which(prlp == FALSE)
        
        oldval <- order(wolp)
        
        oldval <- oldval[1:(redo-nanum)]
        
        swap <- min(c(length(oldval), length(newval)))
        
        for (k in 1:swap){
          
          ansset$periods[oldval[k]] <- ansset2$periods[newval[k]]
          
          ansset$permin[oldval[k]] <- ansset2$permin[newval[k]]
          
          ansset$permax[oldval[k]] <- ansset2$permax[newval[k]]
          
        }
        
      }
      
      print(paste("Printing cluster means for the ", redo, " Periods.", sep=""))
      
      print(ansset$periods)
      
      print("Using these cluster means for a fine resolution genetic optimisation...")
      
      for (k in 1:redo){
        
        print(paste("Initiating Genetic Algorithm...       [", k, "/", redo, "]", sep = ""))
        
        ans <- LCGEN(ts, fit = fit, melc = meanmag, amplc = amplitude, from = 1/ansset$permax[k], to = 1/ansset$permin[k], jit = jit, noimod = noimod, seedno = seedno, pop = pop, pairups = ps, nogen = 50, crossover = rep(0.65, 50), mutation = rep(0.03, 50), fdif = 0.6, dfrac = 0.7, clipcut = clipcut, bins = bins, tenper = tenper)
        
        per[k] <- as.numeric(ans$period)
        
        fstat[k] <- as.numeric(ans$fstat)
        
        models[[k]] <- ans$model
        
        loctest <- rep(0, length(ans$othmodel))
        
        for (z in 1:length(ans$othmodel)){
          
          loctest[z] <- try(vuongtest(ans$model, ans$othmodel[[z]])$LRTstat, TRUE)
          
          if (is.na(as.numeric(as.vector(loctest[z])))){
            
            loctest[z] <- 1e+12
            
          }
          
        }
        
        if (min(loctest) < 0){
          
          if (!any(abs(spurs - ans$othper[[which.min(loctest)]]) < 0.01 * spurs)){
            
            print(paste("Replacing suspected local minimum: ", per[k], ", with better candidate period: ", ans$othper[[which.min(loctest)]], ".", sep=""))
            
            per[k] <- as.numeric(ans$othper[[which.min(loctest)]])
            
            models[[k]] <- ans$othmodel[[which.min(loctest)]]
            
          }
          
        }
        
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
      
      vuongday <- as.numeric(vuongtest(models[[which.max(conf)]], ans$daymodel)$LRTstat)
      
      print(paste("Vuong Closeness Test maximized Period found: ", as.character(per[which.max(conf)]), ".", sep=""))
      
      print(paste("Vuong Closeness Test LRT for this period against a constant model: ", as.character(vuongconf), sep=""))
      
      print(paste("Vuong Closeness Test LRT for this period against a one-day periodic model: ", as.character(vuongday), sep=""))
      
      features[2] <- as.character(per[1])
      
      if (redo > 1){
        
        for (k in 2:redo){
          
          features[2] <- paste(features[2], as.character(per[k]), sep = " ")
          
        }
        
      }
      
      closespur <- spurs[which.min(abs(spurs - (1/per[which.max(conf)])))]
      
      features[3] <- -1.0*min(fstat)
      
      features[4] <- per[which.min(fstat)]
      
      features[5] <- per[which.max(conf)]
      
      features[6] <- vuongconf
      
      features[7] <- vuongday
      
      features[8] <- numspurs
      
      features[9] <- closespur
      
      featvec[i,] <- features
      
    }
    
    featvec <- as.data.frame(featvec)
    
    featvec <- cbind(objects$usnoref, featvec)
    
    colnames(featvec, do.NULL = FALSE)
    
    colnames(featvec) <- c("USnoref", "#.Observations", "Output.Periods", "Fit.Stat", "Genetic.Period", "Vuong.Period", "Vuong.Confidence", "Vuong.Day.Confidence", "Num.Spurious", "Nearest.Spurious")
    
    close(channel)
    
    print(Sys.time() - start)
    
    gc()
    
    if (geo == TRUE){
      
      featvec <- featvec[1,]
      
    }
    
    featvec
    
  }
  
}