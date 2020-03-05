bandmethodmatch <- function(data, lim = 100, ep = 0.01, ofac = 10, nt = 150, b = 3, usenyq = FALSE) {
  
  start <- Sys.time()
  
  library(RODBC)
  
  eps <- ep
  
  data$Type <- factor(data$Type)
  
  data$Class <- factor(data$Class)
  
  par(mfrow=c(1, 1))
  
  print("Computing candidate periods for all objects using the band method")
  
  print("and cross-matching to their AAVSO periods...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  predfreq <- matrix(0, nrow = length(data$Name), ncol = nt*b)
  
  permatch <- rep(0, length(data$Name))
  
  for (ob in 1:length(data$Name)){
    
    print(paste("Processing object ", ob, " out of ", length(data$Name), ", named '", data$Name[ob], "' of class '", data$Type[ob], "'.", sep=""))
    
    radra <- radii / abs(cos(DEC[ob]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[ob]-radra, "' AND '", RA[ob]+radra, "' AND DEClin BETWEEN '", DEC[ob]-radii, "' AND '",
                                       DEC[ob]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[ob] - objects$DEClin)^2.0 + ((RA[ob] - objects$RA) / abs(cos(DEC[ob]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- matrix(ts, ncol = 3)
    
    ts <- unique(ts)
    
    ts <- ts[order(ts[,1]),]
    
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
    
    ovsm <- ofac
    
    sper <- spurper(ts, from = freqmin, to = freqmax, ofac = ovsm, name = "", plot = FALSE)
    
    spurs <- sper$scanned[order(sper$power, decreasing = TRUE)]
    
    spurs <- spurs[which(sper$power > (max(sper$power)/3))]
    
    filt <- c((1 / 29.53), spurs)
    
    freq <- bands(ts, filt, ofac = ofac, from = freqmin, to = freqmax, nt = nt, b = b)
    
    for (z in 1:length(freq)){
      
      predfreq[ob,z] <- freq[z]
      
    }
    
    for (z in (length(freq)+1):(nt*b)){
      
      predfreq[ob,z] <- freq[1]
      
    }
    
    if (length(which(any((abs(data$Period[ob] - (1 / predfreq[ob,]))) < ep * data$Period[ob]))) > 0){
      
      permatch[ob] <- 1
      
      print(paste("The AAVSO Period of ", data$Period[ob], " days IS present to within +/-", ep*100, "% tolerance.", sep = ""))
      
    }
    else{
      
      permatch[ob] <- 0
      
      print(paste("The AAVSO Period of ", data$Period[ob], " days IS NOT present to within +/-", ep*100, "% tolerance.", sep = ""))
      
    }
    
  }
  
  print("--------------------------------------------")
  
  print("Matching periods from AAVSO to those from bands-method candidates...")
  
  print(paste("A period will be accepted if within +/-", ep*100, "% of the AAVSO value.", sep=""))
  
  print("--------------------------------------------")
  
  for (i in 1:length(levels(data$Class))){
    
    print(paste("Computing matches for the ", levels(data$Class)[i], " class with AAVSO period data.", sep=""))
    
    totnum <- length(data$Class[which(data$Class == levels(data$Class)[i])])
    
    print(paste("There is a total number of ", totnum, " objects."))
    
    matchnum <- 0
    
    multnum <- 0
    
    for (j in 1:length(data$Name)){
      
      curdata <- data[j,]
      
      if (curdata$Class == levels(data$Class)[i]){
        
        matchnum <- matchnum + length(which(any((abs(curdata$Period - (1 / predfreq[j,]))) < ep * curdata$Period)))
        
        multnum <- multnum + length(which(((1 / predfreq[j,]) > curdata$Period) &
                                            (abs(((1 / predfreq[j,]) / curdata$Period) - floor((1 / predfreq[j,]) / curdata$Period)) < ep))) + 
          length(which(((1 / predfreq[j,]) < curdata$Period) &
                         (abs((curdata$Period / (1 / predfreq[j,])) - floor(curdata$Period / (1 / predfreq[j,]))) < ep)))
        
      }
      
    }
    
    multnum <- multnum - matchnum
    
    #matchnum <- length(data$Period[which(((abs(data$Period - (1 / predfreq[i,]))) < ep * data$Period) & (data$Class == levels(data$Class)[i]))])
    
    #multnum <- length(data$Period[which(((1 / predfreq[i,]) > data$Period) &
    #                                      (abs(((1 / predfreq[i,]) / data$Period) - floor((1 / predfreq[i,]) / data$Period)) < ep) & (data$Type == levels(data$Type)[i]))]) + 
    #  length(data$Period[which(((1 / data$`#1 Freq`) < data$Period) &
    #                             (abs((data$Period / (1 / data$`#1 Freq`)) - floor(data$Period / (1 / data$`#1 Freq`))) < ep) & (data$Type == levels(data$Type)[i]))]) - matchnum
    
    print(paste("The number of objects fulfilling the +/-", ep*100, "% requirement is ", matchnum, ".", sep=""))
    
    #print(paste("The number of multiples fulfilling the +/-", ep*100, "% requirement is ", multnum, ".", sep=""))
    
    permat <- (matchnum / totnum) * 100
    
    permult <- (multnum / totnum) * 100
    
    perall <- ((matchnum + multnum) / totnum) * 100
    
    print(paste("The matches are ", permat, "% of the total number.", sep=""))
    
    #(paste("The multiples are ", permult, "% of the total number.", sep=""))
    
    #print(paste("Combined this is ", perall, "% of the total number.", sep=""))
    
    #print(paste("The miss rate is ", 100 - perall, "% of the total number.", sep=""))
    
    print(paste("The miss rate is ", 100 - permat, "% of the total number.", sep=""))
    
    print("--------------------------------------------")
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  permatch
  
}