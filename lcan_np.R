lcan_np <- function(ra, declin, radii = 0.01, lim = 2) {
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  
  start <- Sys.time()
  
  set.seed(20)
  
  if (lim < 2){
    
    print("Limit too low, changing to '2' and continuing...")
    
    lim = 2
    
  }
  
  radra <- radii / abs(cos(declin))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries FROM objdat WHERE RA BETWEEN '", 
                                     ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                     declin+radii, "' AND entries > '", lim, "'", sep=""))
  
  if (length(objects$entries) < 1){
    
    close(channel)
    
    stop(paste("No objects detected within search area of ", radii, " degree(s) with greater than ", lim, " entries. Function will now exit.", sep=""))
    
  }
  
  else{
    
    featvec <- matrix(0, length(objects$entries), 26)
    
    for (i in 1:length(objects$entries)){
      
      info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                      objects$usnoref[i], "'", sep=""))
      
      ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
      
      ts <- ts[order(ts[,1]),]
      
      features <- rep.int(0, 26)
      
      features[1] <- objects$entries[i]
      
      stddev <- sd(ts[,2])
      
      features[8] <- stddev
      
      features[6] <- skewness(ts[,2])
      
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
      
      features[2] <- meanmag
      
      features[3] <- stddev / meanmag
      
      features[4] <- mean(w) * ((ts[(n-1),1] - ts[1,1]) ^ 2) * (wn / ((stddev ^ 2) * sum(w) * n^2))
      
      quan <- quantile(ts[,2])
      
      features[5] <- as.numeric(quan[4]) - as.numeric(quan[2])
      
      features[7] <- (((n*(n+1) / ((n-1)*(n-2)*(n-3))) * sum(((ts[,2] - meanmag) / stddev) ^ 4)) - ((3*(n-1)^2) / ((n-2)*(n-3))))
      
      cs <- (1 / (n * stddev)) * cumsum(ts[,2] - meanmag)
      
      features[9] <- max(cs) - min(cs)
      
      features[10] <- (length(ts[,2][ts[,2] > (meanmag + stddev)]) + length(ts[,2][ts[,2] < (meanmag - stddev)])) / length(ts[,2])
      
      sigmap <- (sqrt(n / (n - 1)) * (ts[,2] - meanmag) / (1 / ts[,3]))
      
      features[11] <- (1 / sqrt(n) * sum(abs(sigmap)) / sqrt(sum(sigmap ^ 2)))
      
      for (k in 2:n){
        
        slope <- (ts[k,2] - ts[(k-1),2]) / (ts[k,1] - ts[(k-1),1])
        
        if ((abs(slope) > abs(features[13])) & !is.na(slope)){
          
          features[13] <- slope
          
        }
        
      }
      
      amplitude <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
      
      features[14] <- amplitude
      
      features[15] <- median(abs(ts[,2] - median(ts[,2])))
      
      medianmag <- median(ts[,2])
      
      features[16] <- 1 - ((length(ts[,2][ts[,2] > (medianmag + (amplitude/10))]) + length(ts[,2][ts[,2] < (medianmag - (amplitude/10))])) / length(ts[,2]))
      
      if (n >= 30){
        
        datalast <- ts[,2][length(ts[,2]) - (29:0)]
        
        features[17] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / 30)
        
      }
      else
      {
        
        datalast <- ts[,2][length(ts[,2]) - ((n-1):0)]
        
        features[17] <- ((length(datalast[diff(datalast) > 0]) - length(datalast[diff(datalast) <= 0])) / n)
        
      }
      
      sortedmags <- sort(ts[,2])
      
      F_60_index <- ceiling(0.60 * n)
      
      F_40_index <- ceiling(0.40 * n)
      
      F_5_index <- ceiling(0.05 * n)
      
      F_95_index <- ceiling(0.95 * n)
      
      F_40_60 <- sortedmags[F_60_index] - sortedmags[F_40_index]
      
      F_5_95 <- sortedmags[F_95_index] - sortedmags[F_5_index]
      
      features[18] <- F_40_60 / F_5_95
      
      F_325_index <- ceiling(0.325 * n)
      
      F_675_index <- ceiling(0.675 * n)
      
      F_325_675 <- sortedmags[F_675_index] - sortedmags[F_325_index]
      
      features[19] <- F_325_675 / F_5_95
      
      F_25_index <- ceiling(0.25 * n)
      
      F_75_index <- ceiling(0.75 * n)
      
      F_25_75 <- sortedmags[F_75_index] - sortedmags[F_25_index]
      
      features[20] <- F_25_75 / F_5_95
      
      F_175_index <- ceiling(0.175 * n)
      
      F_825_index <- ceiling(0.825 * n)
      
      F_175_825 <- sortedmags[F_825_index] - sortedmags[F_175_index]
      
      features[21] <- F_175_825 / F_5_95
      
      F_10_index <- ceiling(0.10 * n)
      
      F_90_index <- ceiling(0.90 * n)
      
      F_10_90 <- sortedmags[F_90_index] - sortedmags[F_10_index]
      
      features[22] <- F_10_90 / F_5_95
      
      distmed <- abs(ts[,2] - medianmag)
      
      maxdist <- max(distmed)
      
      features[23] <- maxdist / medianmag
      
      features[24] <- F_5_95 / medianmag
      
      adnt <- ad.test(ts[,2])
      
      features[25] <- (1 / (1.0 + exp(-10.0 * (as.numeric(adnt$statistic) - 0.3))))
      
      lagval <- 100
      
      autocorlen <- 0
      
      while (autocorlen == 0){
        
        AC <- acf(ts[,2], lag.max = lagval, type = "correlation", plot = FALSE)
        
        racf <- as.numeric(AC$acf)
        
        autocorlen <- match(max(racf[racf < exp(-1)]), racf)
        
        features[26] <- autocorlen
        
        if (is.na(autocorlen)){
          
          autocorlen <- 0
          
        }
        
        lagval <- lagval + 100
        
      }
      
      nac <- length(racf)
      
      sigmapac <- sqrt(nac / (nac - 1)) * ((racf - mean(racf)) / sd(racf))
      
      features[12] <- (1 / sqrt(n) * sum(abs(sigmapac)) / sqrt(sum(sigmapac ^ 2)))
      
      featvec[i,] <- features
      
    }
    
  }
  
  featvec <- as.data.frame(featvec)
  
  featvec <- cbind(objects$usnoref, featvec)
  
  colnames(featvec, do.NULL = FALSE)
  
  colnames(featvec) <- c("USnoref","# Observations", "Mean Magnitude", "Mean Variance", "Variability Index",
                         "Q31", "Skewness", "Small Kurtosis", "StD", "Rcs", "Beyond 1 StD", "Stetson Kurtosis Measure",
                         "Stetson-K AC", "Maximum Slope", "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage",
                         "Pair Slope Trend", "Flux Percentile Ratio mid20", "Flux Percentile Ratio mid35",
                         "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                         "Percent Amplitude", "Percent Difference Flux Percentile", "Anderson-Darling Statistic", "Autocorrelation Length")
  
  close(channel)
  
  featvec
  
}