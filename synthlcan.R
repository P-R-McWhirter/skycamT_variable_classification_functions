synthlcan_grape <- function(timeset, perset, snr = 2, type = "sin", fit = "bgls", redo = 5, lcseed = 10, seedno = 1, fulltest = TRUE, quiet = FALSE){
  
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
  
  amp <- 0.2*sqrt(2*snr)
  
  featvec <- matrix(0, nrow(timeset), 9)
  
  for (i in 1:nrow(timeset)){
    
    print(paste("Processing object ", i, "/", nrow(timeset), ".", sep=""))
    
    set.seed(lcseed)
    
    freqmin <- 1 / (max(timeset[i,], na.rm = T) - min(timeset[i,], na.rm = T))
    
    if (type == "non"){
      
      ts <- cbind(timeset[i,], 0.2*rnorm(length(timeset[i,]), 0, 1), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
    }
    else if (type == "sin"){
      
      ts <- cbind(timeset[i,], (amp*(sin(2*pi*timeset[i,]/perset[i]))+0.2*rnorm(length(timeset[i,]), 0, 1)), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
    }
    else if (type == "saw"){
      
      ts <- cbind(timeset[i,], (2*amp*((timeset[i,]/perset[i])-floor(timeset[i,]/perset[i])) - amp + 0.2*rnorm(length(timeset[i,]), 0, 1)), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
    }
    else if (type == "ecl"){
      
      prim <- 2*amp
      
      sec <- 0.5 * prim
      
      ts <- cbind(timeset[i,], 0.2*rnorm(length(timeset[i,]), 0, 1), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
      eclfold <- (timeset[i,]/perset[i]) %% 1
      
      eclstate <- eclfold >= 0.45 & eclfold <= 0.55
      
      eclstate2 <- eclfold <= 0.1 | eclfold >= 0.9
      
      ts[which(eclstate),2] <- ts[which(eclstate),2] - sec + abs(0.5 - eclfold[which(eclstate)])*(20*sec)
      
      ts[which(eclstate2 & eclfold <= 0.1),2] <- ts[which(eclstate2 & eclfold <= 0.1),2] - prim + abs(0.0 - eclfold[which(eclstate2 & eclfold <= 0.1)])*(10*prim)
      
      ts[which(eclstate2 & eclfold >= 0.9),2] <- ts[which(eclstate2 & eclfold >= 0.9),2] - prim + abs(1.0 - eclfold[which(eclstate2 & eclfold >= 0.9)])*(10*prim)
      
    }
    else if (type == "eeb"){
      
      prim <- 2*amp
      
      sec <- 0.5 * prim
      
      ts <- cbind(timeset[i,], 0.2*rnorm(length(timeset[i,]), 0, 1), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
      eclfold <- (timeset[i,]/perset[i]) %% 1
      
      eclstate <- eclfold >= 0.65 & eclfold <= 0.75
      
      eclstate2 <- eclfold <= 0.1 | eclfold >= 0.9
      
      ts[which(eclstate),2] <- ts[which(eclstate),2] - sec + abs(0.7 - eclfold[which(eclstate)])*(20*sec)
      
      ts[which(eclstate2 & eclfold <= 0.1),2] <- ts[which(eclstate2 & eclfold <= 0.1),2] - prim + abs(0.0 - eclfold[which(eclstate2 & eclfold <= 0.1)])*(10*prim)
      
      ts[which(eclstate2 & eclfold >= 0.9),2] <- ts[which(eclstate2 & eclfold >= 0.9),2] - prim + abs(1.0 - eclfold[which(eclstate2 & eclfold >= 0.9)])*(10*prim)
      
    }
    else{
      
      print(paste("'", type, "' is an invalid variability type."))
      
    }
    
    plot(ts[,1:2], pch=19, main = "Raw Light Curve", xlab = "Time (MJD)", ylab = "Magnitude")
    
    plot((ts[,1]/perset[i])%%1, ts[,2], pch=19, xlim = c(0, 1))
    
    features <- rep(0, 9)
    
    features[1] <- length(ts[,1])
    
    meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
    
    amplitude <- abs(max(ts[,2]) - min(ts[,2]))/2.0
    
    tsts <- ts
    
    per <- rep(0, redo)
    
    fstat <- rep(0, redo)
    
    models <- list()
    
    print("Computing the Spurious Periods using a Lomb-Scargle Periodogram...")
    
    sper <- spurperlsp(ts, from = freqmin, to = 1, ofac = 5, name = "", plot = F)
    
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
    
    ansset <- LCMGEN(ts, fit = fit, melc = meanmag, amplc = amplitude, from = freqmin, to = 20, jit = 0.2*2*amplitude, noimod = 1, seedno = seedno, redo = redo, pop = 200, pairups = 50, nogen = 100, crossover = rep(0.65, 100), mutation = seq(0.8, 0.0, length.out = 100), fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20)
    
    ansset2 <- LCMGEN(ts, fit = fit, melc = meanmag, amplc = amplitude, from = freqmin, to = 20, jit = 0.2*2*amplitude, noimod = 1, seedno = seedno+100, redo = redo, pop = 200, pairups = 50, nogen = 100, crossover = rep(0.65, 100), mutation = seq(0.8, 0.0, length.out = 100), fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20)
    
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
      
      ans <- LCGEN(ts, fit = fit, melc = meanmag, amplc = amplitude, from = 1/ansset$permax[k], to = 1/ansset$permin[k], jit = 0.2*2*amplitude, noimod = 1, seedno = seedno, pop = 200, pairups = 50, nogen = 50, crossover = rep(0.65, 50), mutation = rep(0.03, 50), fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20, tenper = F)
      
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
    
    nearspur <- spurs[which.min(abs(spurs - as.numeric(1/per[1])))]
    
    perspurlvl <- exp(-1.0 * ((as.numeric(1/per[1]) - nearspur)^2.0 / 2.0*(100^2.0)))
    
    perspurlvl <- -log10(perspurlvl)
    
    features[2] <- as.character(per[1])
    
    if (redo > 1){
      
      for (k in 2:redo){
        
        features[2] <- paste(features[2], as.character(per[k]), sep = " ")
        
      }
      
    }
    
    features[3] <- -1.0*min(fstat)
    
    features[4] <- per[which.min(fstat)]
    
    features[5] <- per[which.max(conf)]
    
    features[6] <- vuongconf
    
    features[7] <- vuongday
    
    features[8] <- numspurs
    
    features[9] <- perspurlvl
    
    featvec[i,] <- features
    
  }
  
  featvec <- as.data.frame(featvec)
  
  featvec <- cbind(perset, featvec)
  
  colnames(featvec, do.NULL = FALSE)
  
  colnames(featvec) <- c("Input.Period", "#.Observations", "Output.Periods", "Fit.Stat", "Genetic.Period", "Vuong.Period", "Vuong.Confidence", "Vuong.Day.Confidence", "Num.Spurious", "Period.Spurious.Level")
  
  runtime <- (Sys.time() - start)
  
  print(runtime)
  
  gc()
  
  out <- list(featvec = featvec, runtime = runtime)
  
  out
  
}






synthlcan_pergm <- function(timeset, perset, snr = 2, type = "sin", redo = 3, lcseed = 10, fulltest = TRUE, quiet = FALSE){
  
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
  
  amp <- 0.2*sqrt(2*snr)
  
  featvec <- matrix(0, nrow(timeset), 9)
  
  for (i in 1:nrow(timeset)){
    
    print(paste("Processing object ", i, "/", nrow(timeset), ".", sep=""))
    
    set.seed(lcseed)
    
    freqmin <- 1 / (max(timeset[i,], na.rm = T) - min(timeset[i,], na.rm = T))
    
    if (type == "non"){
      
      ts <- cbind(timeset[i,], 0.2*rnorm(length(timeset[i,]), 0, 1), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
    }
    else if (type == "sin"){
      
      ts <- cbind(timeset[i,], (amp*(sin(2*pi*timeset[i,]/perset[i]))+0.2*rnorm(length(timeset[i,]), 0, 1)), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
    }
    else if (type == "saw"){
      
      ts <- cbind(timeset[i,], (2*amp*((timeset[i,]/perset[i])-floor(timeset[i,]/perset[i])) - amp + 0.2*rnorm(length(timeset[i,]), 0, 1)), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
    }
    else if (type == "ecl"){
      
      prim <- 2*amp
      
      sec <- 0.5 * prim
      
      ts <- cbind(timeset[i,], 0.2*rnorm(length(timeset[i,]), 0, 1), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
      eclfold <- (timeset[i,]/perset[i]) %% 1
      
      eclstate <- eclfold >= 0.45 & eclfold <= 0.55
      
      eclstate2 <- eclfold <= 0.1 | eclfold >= 0.9
      
      ts[which(eclstate),2] <- ts[which(eclstate),2] - sec + abs(0.5 - eclfold[which(eclstate)])*(20*sec)
      
      ts[which(eclstate2 & eclfold <= 0.1),2] <- ts[which(eclstate2 & eclfold <= 0.1),2] - prim + abs(0.0 - eclfold[which(eclstate2 & eclfold <= 0.1)])*(10*prim)
      
      ts[which(eclstate2 & eclfold >= 0.9),2] <- ts[which(eclstate2 & eclfold >= 0.9),2] - prim + abs(1.0 - eclfold[which(eclstate2 & eclfold >= 0.9)])*(10*prim)
      
    }
    else if (type == "eeb"){
      
      prim <- 2*amp
      
      sec <- 0.5 * prim
      
      ts <- cbind(timeset[i,], 0.2*rnorm(length(timeset[i,]), 0, 1), rep(0.2, length(timeset[i,])))
      
      ts <- ts[complete.cases(ts),]
      
      eclfold <- (timeset[i,]/perset[i]) %% 1
      
      eclstate <- eclfold >= 0.65 & eclfold <= 0.75
      
      eclstate2 <- eclfold <= 0.1 | eclfold >= 0.9
      
      ts[which(eclstate),2] <- ts[which(eclstate),2] - sec + abs(0.7 - eclfold[which(eclstate)])*(20*sec)
      
      ts[which(eclstate2 & eclfold <= 0.1),2] <- ts[which(eclstate2 & eclfold <= 0.1),2] - prim + abs(0.0 - eclfold[which(eclstate2 & eclfold <= 0.1)])*(10*prim)
      
      ts[which(eclstate2 & eclfold >= 0.9),2] <- ts[which(eclstate2 & eclfold >= 0.9),2] - prim + abs(1.0 - eclfold[which(eclstate2 & eclfold >= 0.9)])*(10*prim)
      
    }
    else{
      
      print(paste("'", type, "' is an invalid variability type."))
      
    }
    
    plot(ts[,1:2], pch=19, main = "Light Curve", xlab = "Generation", ylab = "Magnitude")
    
    plot((ts[,1]/perset[i])%%1, ts[,2], pch=19, xlim = c(0, 1))
    
    features <- rep(0, 9)
    
    features[1] <- length(ts[,1])
    
    meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
    
    amplitude <- abs(max(ts[,2]) - min(ts[,2]))/2.0
    
    tsts <- ts
    
    per <- rep(0, redo)
    
    fstat <- rep(0, redo)
    
    models <- list()
    
    print("Computing the Spurious Periods using a Lomb-Scargle Periodogram...")
    
    sper <- spurperlsp(ts, from = freqmin, to = 1, ofac = 5, name = "", plot = F)
    
    spurs <- sper$scanned[order(sper$power, decreasing = TRUE)]
    
    spurpwr <- sper$power[order(sper$power, decreasing = TRUE)]
    
    spurs <- spurs[which(spurpwr > (max(spurpwr)*0.5))]
    
    numspurs <- length(spurs)
    
    #spurs <- c((1 / 29.53), spurs, 2*spurs)
    
    spurs <- c(1/0.9972, 1, 1/1.00274, spurs)
    
    spurs <- c(spurs*3, spurs*2, spurs, spurs/2, spurs/3, 1/29.53, abs((29.53^-1 + 1)^-1), abs((29.53^-1 - 1)^-1))
    
    spurs <- sort(spurs)
    
    print(paste("Number of detected Spurious Periods: ", numspurs, ".", sep=""))
    
    print("Initiating Bayesian Generalised Lomb-Scargle Periodogram...")
    
    ans <- BGLSgen(ts, plow = 0.05, phigh = 1/freqmin, ofac = 5, dt = NULL, lent = NULL, spur = F, jit = 0.2*2*amplitude, plot = F)
    
    per <- 1/ans$f[order(ans$logp, decreasing = T)][1:redo]
    
    fstat <- ans$logp[order(ans$logp, decreasing = T)][1:redo]
    
    print(per)
    
    print(fstat)
    
    for (z in 1:redo){
      
      ansref <- try(BGLSgen(ts, plow = 0.9*as.numeric(per[z]), phigh = 1.1*as.numeric(per[z]), ofac = 20, dt = NULL, lent = NULL, spur = F, jit = 0.2*2*amplitude, plot = F), TRUE)
      
      if (!(class(ansref) == "try-error")){
        
        per[z] <- 1/ansref$f[order(ansref$logp, decreasing = T)][1]
        
      }
      
    }
    
    print(per)
    
    print(fstat)
    
    vomds <- vuongmod(ts, fstat, per)
    
    for (k in 1:redo){
      
      models[[k]] <- vomds$model[[k]]
      
      loctest <- rep(0, length(vomds$othmodel[[k]]))
      
      for (z in 1:length(vomds$othmodel[[k]])){
        
        loctest[z] <- try(vuongtest(models[[k]], vomds$othmodel[[k]][[z]])$LRTstat, TRUE)
        
        if (is.na(as.numeric(as.vector(loctest[z])))){
          
          loctest[z] <- 1e+12
          
        }
        
      }
      
      if (min(loctest) < 0){
        
        if (!any(abs(spurs - vomds$othper[[k]][[which.min(loctest)]]) < 0.01 * spurs)){
          
          print(paste("Replacing suspected local minimum: ", per[k], ", with better candidate period: ", vomds$othper[[k]][[which.min(loctest)]], ".", sep=""))
          
          per[k] <- as.numeric(vomds$othper[[k]][[which.min(loctest)]])
          
          models[[k]] <- vomds$othmodel[[k]][[which.min(loctest)]]
          
        }
        
      }
      
    }
    
    print(paste("Completed all ", redo, " runs...", sep = ""))
    
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
        
        conf[k] <- vuongtest(models[[k]], vomds$onemodel)$LRTstat
        
      }
      
    }
    
    vuongconf <- as.numeric(vuongtest(models[[which.max(conf)]], vomds$onemodel)$LRTstat)
    
    vuongday <- as.numeric(vuongtest(models[[which.max(conf)]], vomds$daymodel)$LRTstat)
    
    print(paste("Vuong Closeness Test maximized Period found: ", as.character(per[which.max(conf)]), ".", sep=""))
    
    print(paste("Vuong Closeness Test LRT for this period against a constant model: ", as.character(vuongconf), sep=""))
    
    print(paste("Vuong Closeness Test LRT for this period against a one-day periodic model: ", as.character(vuongday), sep=""))
    
    nearspur <- spurs[which.min(abs(spurs - as.numeric(1/per[1])))]
    
    perspurlvl <- exp(-1.0 * ((as.numeric(1/per[1]) - nearspur)^2.0 / 2.0*(100^2.0)))
    
    perspurlvl <- -log10(perspurlvl)
    
    features[2] <- as.character(per[1])
    
    if (redo > 1){
      
      for (k in 2:redo){
        
        features[2] <- paste(features[2], as.character(per[k]), sep = " ")
        
      }
      
    }
    
    features[3] <- -1.0*min(fstat)
    
    features[4] <- per[1]
    
    features[5] <- per[which.max(conf)]
    
    features[6] <- vuongconf
    
    features[7] <- vuongday
    
    features[8] <- numspurs
    
    features[9] <- perspurlvl
    
    featvec[i,] <- features
    
  }
  
  featvec <- as.data.frame(featvec)
  
  featvec <- cbind(perset, featvec)
  
  colnames(featvec, do.NULL = FALSE)
  
  colnames(featvec) <- c("Input.Period", "#.Observations", "Output.Periods", "Fit.Stat", "Genetic.Period", "Vuong.Period", "Vuong.Confidence", "Vuong.Day.Confidence", "Num.Spurious", "Period.Spurious.Level")
  
  runtime <- (Sys.time() - start)
  
  print(runtime)
  
  gc()
  
  out <- list(featvec = featvec, runtime = runtime)
  
  out
  
}





vuongmod <- function(ts, logp, per){
  
  n <- length(logp)
  
  f <- 1/per
  
  SSTlm <- list()
  
  othmodel <- list()
  
  othper <- list()
  
  tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, ts)
  
  tsin <- tsin[order(tsin[,1]),]
  
  al <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3)
  
  mul <- 1:3
  
  othmodel <- list()
  
  othper <- list()
  
  for (k in 1:n){
    
    SSTlm[[k]] <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*f[k]*tsin[,1]) + cos(2*pi*f[k]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
    
    fs <- matrix(0, nrow = length(mul), ncol = length(al))
    
    for (i in 1:length(mul)){
      
      for (j in 1:length(al)){
        
        fs[i,j] <- abs(mul[i]*(0.99726957 / (0.99726957 / (f[k]^-1) + al[j])))^-1
        
      }
      
    }
    
    fs <- c(fs)
    
    othmodelk <- list()
    
    for (i in 1:length(fs)){
      
      othmodelk[[i]] <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*fs[i]*tsin[,1]) + cos(2*pi*fs[i]*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
      
    }
    
    othmodel[[k]] <- othmodelk
    
    othperk <- 1 / fs
    
    othper[[k]] <- othperk
    
  }
  
  onemodel <- lm(tsin[,2] ~ 1, data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  daymodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*(0.99726957^-1)*tsin[,1]) + cos(2*pi*(0.99726957^-1)*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  sp.out <- list(model = SSTlm, onemodel = onemodel, daymodel = daymodel, othmodel = othmodel, othper = othper)
  
}