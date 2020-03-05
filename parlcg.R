parlcg <- function(ra, declin, radii = 0.1, permin = 0.1, permax = 1000, lim = 2, ovsm = 1, rang = 2, rand = FALSE, ploterr = FALSE) {
  
  start <- Sys.time()
  
  library(RODBC)
  library(lomb)
  library(RobPer)
  library(foreach)
  library(doParallel)
  library(parallel)
  
  no_cores <- detectCores() - 2
  
  c1 <- makeCluster(no_cores)
  
  registerDoParallel(c1)
  
  set.seed(12)
  
  period <- seq(from = permin, to = permax, by = 1000)
  
  if (lim < 2){
    
    print("Limit too low, changing to '2' and continuing...")
    
    lim = 2
    
  }
  
  radra <- radii / abs(cos(declin))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries FROM objdat WHERE RA BETWEEN '", 
                                     ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                     declin+radii, "' and entries > '", lim, "'", sep=""))
  
  if (length(objects$entries) < 1){
    
    close(channel)
    stop(paste("No objects detected within search area of ", radii, " degree(s). Function will now exit.", sep=""))
    
  }
  
  else{
    
    print(paste("A total of ", length(objects$usnoref),
                " object(s) identified within the search radius of ", radii, " degree(s) with greater than ", lim, " observations.", sep=""))
    
    for (i in 1:length(objects$entries)){
      
      print(paste("Object with reference id '", objects$usnoref[i],
                  "' successfully read in. Object has ", objects$entries[i], " observation(s).", sep=""))
      
      print("Fitting Lomb-Scargle Periodogram...")
      
      info <- sqlQuery(channel, paste("SELECT * FROM obsdat WHERE usnoref = '", 
                                      objects$usnoref[i], "'", sep=""))
      
      ts <- matrix(c(info$MJD, (info$Rcat - mean(info$Rcat))), nrow = length(info$MJD))
      
      ans <- foreach(kn = 1:length(period)-1, .combine = c, .packages = "lomb", .multicombine = TRUE) %dopar% {
        
        lsp(ts, from = period[kn], to = period[kn+1], type = "period", ofac = ovsm, alpha = 0.01, plot = FALSE)
      
      }
      
      if (rand == TRUE){
        
        print("Using random dataset to remove sampling periodicities...")
        
        tsran <- matrix(c(info$MJD, rnorm(length(info$Rcat), mean = mean(info$Rcat), sd = sd(info$Rcat))), nrow = length(info$MJD))
        
        ansran <- foreach(kn = 1:length(period)-1, .combine = c, .packages = "lomb") %dopar% {
          
          lsp(tsran, from = period[kn], to = period[kn+1], type = "period", ofac = ovsm, plot = FALSE)
        
        }
        
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
  
  stopCluster(c1)
  
  print(Sys.time() - start)
  
  p
  
}