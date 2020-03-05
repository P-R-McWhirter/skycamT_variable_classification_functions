lcg <- function(ra, declin, radii = 0.1, permin = 0.1, permax = 1000, lim = 2, ovsm = 2, rang = 2, rand = FALSE, ploterr = FALSE, alldb = FALSE, returnp = FALSE) {
  
  start <- Sys.time()
  
  library(RODBC)
  library(lomb)
  library(RobPer)
  
  set.seed(12)
  
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
    
    for (i in 1:length(objects$entries)){
      
      print(paste("Object with reference id '", objects$usnoref[i],
                  "' successfully read in. Object has ", objects$entries[i], " observation(s).", sep=""))
      
      info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                      objects$usnoref[i], "'", sep=""))
      
      ts <- matrix(c(info$MJD, (info$Rcat - mean(info$Rcat))), nrow = length(info$MJD))
      
      print("Removing Linear Trends...")
      
      fit <- lm(ts[,2] ~ ts[,1], data = as.data.frame(ts))
      
      coeff <- coefficients(fit)
      
      ts[,2] <- ts[,2] - (coeff[1] + coeff[2] * ts[,1])
      
      print("Fitting Lomb-Scargle Periodogram...")
      
      ans <- lsp(ts, from = permin, to = permax, type = "period", ofac = ovsm, alpha = 0.01, plot = FALSE)
      
      if (rand == TRUE){
        
        print("Using random dataset to remove sampling periodicities...")
        
        tsran <- matrix(c(info$MJD, rnorm(length(info$Rcat), mean = mean(info$Rcat), sd = sd(info$Rcat))), nrow = length(info$MJD))
        
        ansran <- lsp(tsran, from = permin, to = permax, type = "period", ofac = ovsm, plot = FALSE)
        
        pwr <- ans$power
        
        pwrran <- ansran$power[order(pwr, decreasing = T)]
        
        pwr <- pwr - pwrran
        
      }
      else
      {
        pwr <- ans$power
      }
      
      print("Displaying top 10 periods detected:")
      
      p <- ans$scanned[order(pwr, decreasing = T)]
      
      p <- p[findInterval(p, c(0.95,1.05)) != 1L & findInterval(p, c(0.45,0.55)) != 1L & findInterval(p, c(29.4,29.6)) != 1L]
      
      p <- p[1:10]
      
      effm <- 2 * length(ans$scanned) / ovsm
      
      level <- -log(1 - (1 - ans$alpha)^(1/effm))
      
      exPN <- rep.int(0, length(p))
      
      for (z in 1:length(p)){
        
        exPN[z] <- exp(-1.0 * pwr[ans$scanned == p[z]])
        
      }
      
      pv <- effm * exPN
      
      for (z in 1:length(pv)){
        
        if (pv[z] > 0.01){
          
          pv[z] <- 1 - (1 - exPN[z])^effm
          
        }
        
      }
      
      ts[,2] = info$Rcat
      
      print(p)
      
      print("Displaying corresponding p-values:")
      
      print(pv)
      
      print(paste("Standard Deviation of the Apparent Magnitude: ", sd(ts[,2]), sep=""))
      
      print("Plotting phase diagram of strongest period (min p-value)...")
      
      pos <- mean(info$Rcat) + rang * sd(info$Rcat)
      
      neg <- mean(info$Rcat) - rang * sd(info$Rcat)
      
      folded <- ((((info$MJD - info$MJD[which.max(info$MJD)]) / p[which.min(pv)]) + 1.25) %% 2) - 1
      
      plot(folded, info$Rcat, main = paste("", objects$usnoref[i], " Phase Plot", sep=""), xlab = paste("Phase at a period of ", p[which.min(pv)], " days", sep=""),
           ylab = "Apparent Magnitude", pch=19, ylim = c(pos, neg))
      
      if (ploterr == TRUE){
        
        magerr <- info$Rcaterr
        
        segments(folded, info$Rcat-magerr, folded, info$Rcat+magerr)
        
        epsilon = 0.01
        
        segments(folded-epsilon, info$Rcat-magerr, folded+epsilon, info$Rcat-magerr)
        
        segments(folded-epsilon, info$Rcat+magerr, folded+epsilon, info$Rcat+magerr)
        
      }
      
    }
    
  }
  
  close(channel)
  
  print(Sys.time() - start)
  
  if (returnp == TRUE){
    
    p
    
  }
  
}