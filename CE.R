CE <- function(x, ofac = 1, from = 0.001, to = 20, bins = 20, plot = TRUE) {
  
  set.seed(10)
  type = "frequency"
  alpha = 0.01
  level = 0
  
  eps = 1e-3
  
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  yb <- y
  
  n <- length(y)
  tspan <- t[n] - t[1]
  
  df <- 1 / (tspan * ofac)
  
  freq <- seq(from, to, by = df)
  
  m = length(freq)
  
  print(paste("Number of trial periods extracted: ", m, ".", sep=""))
  
  pwr <- rep(0, m)
  
  x[,2] <- x[,2] - mean(x[,2])
  
  me <- mean(y)
  
  tot <- length(x[,1])
  
  for (k in 1:m) {
    
    p <- matrix(0, nrow = bins, ncol = bins)
    
    xts <- x
    
    xts[,1] <- (xts[,1] * freq[k]) %% 1
    
    xts[,2] <- (xts[,2] - min(xts[,2]))
    
    xts[,2] <- xts[,2] / (max(xts[,2]) - min(xts[,2]))
    
    #memg <- as.numeric(tapply(xts[,2], cut(xts[,1], bins), mean))
    
    #meph <- seq(0, 1, length.out = bins+1)
    
    #meph <- meph[1:bins] + diff(meph)/2
    
    #xts <- cbind(meph, memg)
    
    #xts <- xts[complete.cases((xts)),]
    
    ab = matrix(c(0, 0, 1, 1), 2, 2)
    
    datbin <- bin2(cbind(xts[,1], xts[,2]), ab = ab, nbin = c(bins, bins))
    
    p <- datbin$nc / tot
    
    pj <- colSums(p)
    
    po <- log(pj / p)
    
    po[which(po == Inf | po == -Inf | is.nan(po))] <- 0
    
    pwr[k] <- sum(colSums(p * po))
    
  }
  
  #pwr <- pwr * ((n - 1) / (2 * n))
  
  pwr <- max(pwr) - pwr
  
  # Fit a exponential function to cancel out the red noise
  
  #rfit <- 10^seq(log10(from), log10(to), length.out = 51)
  
  #redf <- rep(0, 50)
  
  #for (j in 1:(length(rfit)-1)){
    
    #while (redf[j] == 0 | any(freq == redf[j])){
      
      #redf[j] <- runif(1, rfit[j], rfit[j+1])
      
    #}
    
  #}
  
  #redn <- rep(0, length(redf))
  
  #for (k in 1:length(redf)) {
    
    #x <- xb
    
    #x1 <- xb
    
    #p <- (t / (1 / redf[k])) %% 1
    
    #x[, 1] <- p
    
    #x <- x[order(p),]
    
    #x1 <- x
    
    #for (i in 2:length(x[,1])) {
      
      #x1[i,] <- x[i-1,]
      
    #}
    
    #x1[1,] <- x[length(x[,1]),]
    
    #redn[k] <- sum((x1[,2] - x[,2])^2) / sum((x[,2] - me)^2)
    
  #}
  
  #redn <- redn * ((n - 1) / (2 * n))
  
  #redn <- 1 - redn
  
  #arang <- seq(from = 0, to = 10, by = 0.01)
  
  #mseold <- 10000
  
  #for (k in 1:length(arang)){
    
    #rprof <- lm(redn ~ exp(-1.0*arang[k]*redf), data = as.data.frame(cbind(redf, redn)))
    
    #mse <- sum((rprof$residuals)^2)
    
    #if (mse < mseold){
      
      #mseold <- mse
      
      #besta <- arang[k]
      
    
  #}
  
  #print("Red Noise fitted using exponential profil
      #rnco <- rprof$coefficientse, k*exp(-At) + C.")
      
      #rnco[is.na(rnco)] <- 0
    #}
  
  #print(paste("k = ", rnco[1], ", A = ", besta, ", C 
  #    = ", rnco[2], ".", sep=""))
  
  #pwr <- pwr - (rnco[1] + rnco[2] * exp(-1.0*besta*freq))
  
  peak.at = c(freq[which.max(pwr)], 1/freq[which.max(pwr)])
  
  p = 0
  
  if (plot == TRUE){
    
    plot(freq, pwr, pch=19, type="n", main="Conditional Entropy", xlab = "Frequency (Cycles/day)", ylab = "Conditional Entropy", ylim = c(min(pwr), max(pwr)))
    
    lines(freq, pwr, type="l")
    
  }
  
  sp.out <- list(scanned = freq, power = pwr, data = datanames, n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = max(pwr), peak.at = peak.at, p.value = p)
}