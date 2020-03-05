fsffold <- function(data, lim = 100) {
  
  start <- Sys.time()
  
  par(mfrow=c(1, 1))
  
  print("Folding all objects with the AAVSO period...")
  
  print("Objects without an AAVSO period have been removed.")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  radii <- 0.1
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  per <- as.numeric(as.vector(data$Period))
  
  ansmat <- matrix(0, nrow = length(data$Name), ncol = 26)
  
  #ansmat <- as.data.frame(ansmat)
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Class[i], "' and subtype '", data$Type[i], "'.", sep=""))
    
    radra <- radii / abs(cos(DEC[i]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       RA[i]-radra, "' AND '", RA[i]+radra, "' AND DEClin BETWEEN '", DEC[i]-radii, "' AND '",
                                       DEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    geodist <- sqrt((DEC[i] - objects$DEClin)^2.0 + ((RA[i] - objects$RA) / abs(cos(DEC[i]*pi/180)))^2.0)
    
    mindist <- which.min(geodist)
    
    objects <- objects[mindist,]
    
    info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                    objects$usnoref[1], "'", sep=""))
    
    ts <- cbind(info$MJD, info$Rcat, info$Rcaterr)
    
    ts <- unique(ts)
    
    melc <- weighted.mean(ts[,2], (1/ts[,3]))
    
    amplc <- abs((median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05))) - median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))) / 2)
    
    ts[,2] <- ts[,2] - melc
    
    ts[,2] <- ts[,2] / amplc
    
    ts[,3] <- ts[,3] / amplc 
    
    ts <- ts[which(ts[,2] < 1),]
    
    ts <- ts[which(ts[,2] > -1),]
    
    ts[,1] <- (((ts[,1] / (per[i])) + 0.25) %% 1)
    
    ts <- ts[order(ts[,1]),]
    
    fold <- ts[,1] - 1
    
    folded <- c(fold, ts[,1])
    
    foldts <- c(ts[,2], ts[,2])
    
    magerr <- c(ts[,3], ts[,3])
    
    fts <- cbind(folded, foldts, magerr)
    
    fts <- fts[order(fts[,1]),]
    
    ftsin <- matrix(c(seq(from = -1, to = 1, length.out = 15768), rep(0, times = 15768), rep(1/0, times = 15768)), nrow = 15768)
    
    ftsin <- rbind(ftsin, fts)
    
    ftsin <- ftsin[order(ftsin[,1]),]
    
    #model <- lm(ftsin[,2] ~ 0 + sin(2*pi*ftsin[,1]) + cos(2*pi*ftsin[,1]) +
    #              sin(4*pi*ftsin[,1]) + cos(4*pi*ftsin[,1]) + sin(6*pi*ftsin[,1]) +
    #              cos(6*pi*ftsin[,1]) + sin(8*pi*ftsin[,1]) + cos(8*pi*ftsin[,1]) +
    #              sin(10*pi*ftsin[,1]) + cos(10*pi*ftsin[,1]) + sin(12*pi*ftsin[,1]) + cos(12*pi*ftsin[,1]) +
    #              sin(14*pi*ftsin[,1]) + cos(14*pi*ftsin[,1]) + sin(16*pi*ftsin[,1]) + cos(16*pi*ftsin[,1]), data = as.data.frame(ftsin[,1:2]), weights = 1/ftsin[,3])
    
    #coeff <- model$coefficients
    
    #coeff[is.na(coeff)] <- 0
    
    Z <- cbind(sin(2*pi*ftsin[,1]), cos(2*pi*ftsin[,1]), sin(4*pi*ftsin[,1]), cos(4*pi*ftsin[,1]), sin(6*pi*ftsin[,1]), cos(6*pi*ftsin[,1]), sin(8*pi*ftsin[,1]), cos(8*pi*ftsin[,1]),
               sin(10*pi*ftsin[,1]), cos(10*pi*ftsin[,1]), sin(12*pi*ftsin[,1]), cos(12*pi*ftsin[,1]), sin(14*pi*ftsin[,1]), cos(14*pi*ftsin[,1]), sin(16*pi*ftsin[,1]), cos(16*pi*ftsin[,1]))
    
    Z <- as.matrix(Z)
    
    lambda <- 10^(-3)
    
    M <- t(Z)%*%Z + diag(rep(16,1))*lambda
    
    coeff <- solve(M)%*%t(Z)%*%(ftsin[,2])
    
    #if (abs((1 / freq[k]) - 1.00) < 0.01){
    
    #plot(fts[,1], fts[,2], main = paste("", 1/freq[k], " Folded Light Curve spread from -1.0 to 1.0", sep=""), xlab = paste("Phase at ", (1 / freq[k]), " days", sep=""),
    #   ylab = "Scaled Magnitude", pch=19, ylim = c(1, -1))
    
    #lines(sort(ftsin[,1]), model$fitted.values, col = 2)
    
    #(paste("Period : ", 1/freq[k], " days", sep = ""))
    
    #print(coeff)
    
    #}
    
    A1 <- (coeff[1]^2.0 + coeff[2]^2.0)^0.5
    
    A2 <- (coeff[3]^2.0 + coeff[4]^2.0)^0.5
    
    A3 <- (coeff[5]^2.0 + coeff[6]^2.0)^0.5
    
    A4 <- (coeff[7]^2.0 + coeff[8]^2.0)^0.5
    
    A5 <- (coeff[9]^2.0 + coeff[10]^2.0)^0.5
    
    A6 <- (coeff[11]^2.0 + coeff[12]^2.0)^0.5
    
    A7 <- (coeff[13]^2.0 + coeff[14]^2.0)^0.5
    
    A8 <- (coeff[15]^2.0 + coeff[16]^2.0)^0.5
    
    NPH1 <- atan2(coeff[2], coeff[1])
    
    NPH2 <- atan2(coeff[4], coeff[3]) - NPH1
    
    NPH3 <- atan2(coeff[6], coeff[5]) - NPH1
    
    NPH4 <- atan2(coeff[8], coeff[7]) - NPH1
    
    NPH5 <- atan2(coeff[10], coeff[9]) - NPH1
    
    NPH6 <- atan2(coeff[12], coeff[11]) - NPH1
    
    NPH7 <- atan2(coeff[14], coeff[13]) - NPH1
    
    NPH8 <- atan2(coeff[16], coeff[15]) - NPH1
    
    PH1 <- 0
    
    PH2 <- atan2(sin(NPH2), cos(NPH2))
    
    PH3 <- atan2(sin(NPH3), cos(NPH3))
    
    PH4 <- atan2(sin(NPH4), cos(NPH4))
    
    PH5 <- atan2(sin(NPH5), cos(NPH5))
    
    PH6 <- atan2(sin(NPH6), cos(NPH6))
    
    PH7 <- atan2(sin(NPH7), cos(NPH7))
    
    PH8 <- atan2(sin(NPH8), cos(NPH8))
    
    ansmat[i,] <- c(data$Class[i], as.character(data$Type[i]), per[i], as.numeric(A1), as.numeric(A2), as.numeric(A3), as.numeric(A4),
                    as.numeric(A5), as.numeric(A6), as.numeric(A7), as.numeric(A8),
                    as.numeric(NPH1), as.numeric(NPH2), as.numeric(NPH3), as.numeric(NPH4),
                    as.numeric(NPH5), as.numeric(NPH6), as.numeric(NPH7), as.numeric(NPH8),
                    as.numeric(PH2), as.numeric(PH3), as.numeric(PH4), as.numeric(PH5),
                    as.numeric(PH6), as.numeric(PH7), as.numeric(PH8))
    
  }
  
  colnames(ansmat) <- c("Class", "Type", "Period", "A1", "A2", "A3", "A4", "A5", "A6", "A7",
                        "A8", "NPH1", "NPH2", "NPH3", "NPH4", "NPH5", "NPH6",
                        "NPH7", "NPH8", "PH2", "PH3", "PH4",
                        "PH5", "PH6", "PH7", "PH8")
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  ansmat
  
}