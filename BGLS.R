BGLS <- function(data, radii = 0.001, lim = 100, plow = 0.5, phigh = 100, ofac = 1, spur = F, jit = 0){
  
  starttime <- Sys.time()
  
  ra <- as.numeric(as.vector(data$RA))
  
  declin <- as.numeric(as.vector(data$DEC))
  
  radra <- radii / abs(cos(declin*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                     ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                     declin+radii, "' AND entries > '", lim, "'", sep=""))
  
  geodist <- sqrt((declin - objects$DEClin)^2.0 + (ra - objects$RA)^2.0)
  
  mindist <- which.min(geodist)
  
  objects <- objects[mindist,]
  
  info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                  objects$usnoref[1], "'", sep=""))
  
  info <- unique(info)
  
  #info <- info[which(info[,1] <= min(info[,1])+500),]
  
  #info <- info[sample(1:nrow(info), nrow(info)),]
  
  #info <- info[1:5000,]
  
  #info <- info[order(info[,1]),]
  
  t <- info$MJD
  
  y <- info$Rcat
  
  err <- info$Rcaterr
  
  n_steps <- as.integer(ofac*length(t)*(1/plow - 1/phigh))
  f <- seq(from = 1/phigh, to = 1/plow, length.out = n_steps)
  f <- sort(f)
  
  y <- y - mean(y)
  
  nor <- max(abs(y))
  
  y <- y / nor
  
  err <- err / nor
  
  err2 <- err*err + jit*jit
  
  if (spur == T){
    
    y <- rep(0, length(t))
    
  }
  
  w <- 1/err2
  W <- sum(w)
  
  bigY <- sum(w*y)
  
  logp <- rep(0, length(f))
  p <- rep(0, length(f))
  constants <- rep(0, length(f))
  exponents <- rep(0, length(f))
  
  for (i in 1:length(f)){
    
    wi <- 2 * pi * f[i]
    
    theta <- 0.5 * atan2(sum(w*sin(2*wi*t)), sum(w*cos(2*wi*t)))
    
    st <- (wi * t) - theta
    
    cosx <- cos(st)
    sinx <- sin(st)
    wcosx <- w*cosx
    wsinx <- w*sinx
    
    C <- sum(wcosx)
    S <- sum(wsinx)
    
    YCh <- sum(y*wcosx)
    YSh <- sum(y*wsinx)
    CCh <- sum(wcosx*cosx)
    SSh <- sum(wsinx*sinx)
    
    if (CCh != 0 & SSh != 0){
      
      K <- (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
      
      L <- (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
      
      M <- (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
      
      constants[i] <- 1/(sqrt(CCh*SSh*abs(K)))
      
    }
    else if (CCh == 0){
      
      K <- (S*S - W*SSh)/(2*SSh)
      
      L <- (bigY*SSh - S*YSh)/(SSh)
      
      M <- (YSh*YSh)/(2*SSh)
      
      constants[i] <- 1/(sqrt(SSh*abs(K)))
      
    }
    else if (SSh == 0){
      
      K <- (C*C - W*CCh)/(2*CCh)
      
      L <- (bigY*CCh - C*YCh)/(CCh)
      
      M <- (YCh*YCh)/(2*CCh)
      
      constants[i] <- 1/(sqrt(CCh*abs(K)))
      
    }
    
    if (K > 0){
      
      print("K is positive, This should not happen.")
      
    }
    
    exponents[i] <- (M - ((L*L)/(4*K)))
    
    logp[i] <- log10(constants[i]) + (exponents[i] * log10(exp(1)))
    
    p[i] <- 10^(as.numeric(logp[i]))
    
    if (p[i] < 0.00001){
      
      p[i] <- 0
      
    }
    
  }
  
  p <- p / max(p)
  
  plot(1/f, p, type = "l")
  
  print(1 / f[order(p, decreasing = T)][1:50])
  
  print(p[order(p, decreasing = T)][1:50])
  
  Sys.time() - starttime
  
}




BGLS2 <- function(ts, plow = 0.5, phigh = 100, ofac = 1, dt = NULL, lent = NULL, spur = F, jit = 0, plot = F){
  
  starttime <- Sys.time()
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  if (is.null(dt)){
    
    dt <- max(t) - min(t)
    
  }
  
  if (is.null(lent)){
    
    lent <- length(t)
    
  }
  
  n_steps <- as.integer(ofac*lent*(1/plow - 1/phigh))
  f <- seq(from = 1/phigh, to = 1/plow, length.out = n_steps)
  f <- sort(f)
  
  y <- y - mean(y)
  
  nor <- max(abs(y))
  
  y <- y / nor
  
  err <- err / nor
  
  err2 <- err*err + jit*jit
  
  if (spur == T){
    
    y <- rep(0, length(t))
    
  }
  
  w <- 1/err2
  W <- sum(w)
  
  bigY <- sum(w*y)
  
  logp <- rep(0, length(f))
  p <- rep(0, length(f))
  constants <- rep(0, length(f))
  exponents <- rep(0, length(f))
  
  for (i in 1:length(f)){
    
    wi <- 2 * pi * f[i]
    
    theta <- 0.5 * atan2(sum(w*sin(2*wi*t)), sum(w*cos(2*wi*t)))
    
    st <- (wi * t) - theta
    
    cosx <- cos(st)
    sinx <- sin(st)
    wcosx <- w*cosx
    wsinx <- w*sinx
    
    C <- sum(wcosx)
    S <- sum(wsinx)
    
    YCh <- sum(y*wcosx)
    YSh <- sum(y*wsinx)
    CCh <- sum(wcosx*cosx)
    SSh <- sum(wsinx*sinx)
    
    if (CCh != 0 & SSh != 0){
      
      K <- (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
      
      L <- (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
      
      M <- (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
      
      constants[i] <- 1/(sqrt(CCh*SSh*abs(K)))
      
    }
    else if (CCh == 0){
      
      K <- (S*S - W*SSh)/(2*SSh)
      
      L <- (bigY*SSh - S*YSh)/(SSh)
      
      M <- (YSh*YSh)/(2*SSh)
      
      constants[i] <- 1/(sqrt(SSh*abs(K)))
      
    }
    else if (SSh == 0){
      
      K <- (C*C - W*CCh)/(2*CCh)
      
      L <- (bigY*CCh - C*YCh)/(CCh)
      
      M <- (YCh*YCh)/(2*CCh)
      
      constants[i] <- 1/(sqrt(CCh*abs(K)))
      
    }
    
    if (K > 0){
      
      print("K is positive, This should not happen.")
      
    }
    
    exponents[i] <- (M - ((L*L)/(4*K)))
    
    logp[i] <- log10(constants[i]) + (exponents[i] * log10(exp(1)))
    
    p[i] <- 10^(as.numeric(logp[i]))
    
    if (p[i] < 0.00001){
      
      p[i] <- 0
      
    }
    
  }
  
  norm <- max(p)
  
  p <- p / max(p)
  
  if (plot == T){
    
    plot((1/f), p, type = "l")
    
  }
  
  #print(logp[which.max(p)])
  
  Sys.time() - starttime
  
  sp.out <- list(norm = norm, f = f, p = p, logp = log10(p))
  
}


BGLSgen <- function(ts, plow = 0.5, phigh = 100, ofac = 1, dt = NULL, lent = NULL, spur = F, jit = 0, plot = F){
  
  starttime <- Sys.time()
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  if (is.null(dt)){
    
    dt <- max(t) - min(t)
    
  }
  
  if (is.null(lent)){
    
    lent <- length(t)
    
  }
  
  n_steps <- as.integer(ofac*lent*(1/plow - 1/phigh))
  f <- seq(from = 1/phigh, to = 1/plow, length.out = n_steps)
  f <- sort(f)
  
  y <- y - mean(y)
  
  nor <- max(abs(y))
  
  y <- y / nor
  
  err <- err / nor
  
  err2 <- err*err + jit*jit
  
  if (spur == T){
    
    y <- rep(0, length(t))
    
  }
  
  w <- 1/err2
  W <- sum(w)
  
  bigY <- sum(w*y)
  
  logp <- rep(0, length(f))
  p <- rep(0, length(f))
  constants <- rep(0, length(f))
  exponents <- rep(0, length(f))
  
  for (i in 1:length(f)){
    
    wi <- 2 * pi * f[i]
    
    theta <- 0.5 * atan2(sum(w*sin(2*wi*t)), sum(w*cos(2*wi*t)))
    
    st <- (wi * t) - theta
    
    cosx <- cos(st)
    sinx <- sin(st)
    wcosx <- w*cosx
    wsinx <- w*sinx
    
    C <- sum(wcosx)
    S <- sum(wsinx)
    
    YCh <- sum(y*wcosx)
    YSh <- sum(y*wsinx)
    CCh <- sum(wcosx*cosx)
    SSh <- sum(wsinx*sinx)
    
    if (CCh != 0 & SSh != 0){
      
      K <- (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
      
      L <- (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
      
      M <- (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
      
      constants[i] <- 1/(sqrt(CCh*SSh*abs(K)))
      
    }
    else if (CCh == 0){
      
      K <- (S*S - W*SSh)/(2*SSh)
      
      L <- (bigY*SSh - S*YSh)/(SSh)
      
      M <- (YSh*YSh)/(2*SSh)
      
      constants[i] <- 1/(sqrt(SSh*abs(K)))
      
    }
    else if (SSh == 0){
      
      K <- (C*C - W*CCh)/(2*CCh)
      
      L <- (bigY*CCh - C*YCh)/(CCh)
      
      M <- (YCh*YCh)/(2*CCh)
      
      constants[i] <- 1/(sqrt(CCh*abs(K)))
      
    }
    
    if (K > 0){
      
      print("K is positive, This should not happen.")
      
    }
    
    exponents[i] <- (M - ((L*L)/(4*K)))
    
    logp[i] <- log10(constants[i]) + (exponents[i] * log10(exp(1)))
    
    p[i] <- 10^(as.numeric(logp[i]))
    
    if (p[i] < 0.00001){
      
      p[i] <- 0
      
    }
    
  }
  
  norm <- max(p)
  
  p <- p / max(p)
  
  if (plot == T){
    
    plot((1/f), p, type = "l")
    
  }
  
  finfreq <- f[which.max(p)]
  
  #print(logp[which.max(p)])
  
  tsin <- matrix(c(seq(from = min(ts[,1]), to = max(ts[,1]), length.out = 15768), rep(weighted.mean(ts[,2], (1 / ts[,3])), times = 15768), rep(1/0, times = 15768)), nrow = 15768)
  
  tsin <- rbind(tsin, ts)
  
  tsin <- tsin[order(tsin[,1]),]
  
  SSTlm <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*finfreq*tsin[,1]) + cos(2*pi*finfreq*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  onemodel <- lm(tsin[,2] ~ 1, data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  daymodel <- lm(tsin[,2] ~ tsin[,1] + sin(2*pi*tsin[,1]) + cos(2*pi*tsin[,1]), data = as.data.frame(tsin[,1:2]), weights = 1/tsin[,3])
  
  Sys.time() - starttime
  
  sp.out <- list(norm = norm, f = f, p = p, logp = log10(p), model = SSTlm, onemodel = onemodel, daymodel = daymodel)
  
}


SBGLS <- function(data, radii = 0.001, lim = 100, plow = 0.5, phigh = 100, ofac = 1, jit = 0, obsstart = 50, obsstep = 1){
  
  starttime <- Sys.time()
  
  ra <- as.numeric(as.vector(data$RA))
  
  declin <- as.numeric(as.vector(data$DEC))
  
  radra <- radii / abs(cos(declin*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                     ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                     declin+radii, "' AND entries > '", lim, "'", sep=""))
  
  geodist <- sqrt((declin - objects$DEClin)^2.0 + (ra - objects$RA)^2.0)
  
  mindist <- which.min(geodist)
  
  objects <- objects[mindist,]
  
  info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                  objects$usnoref[1], "'", sep=""))
  
  info <- unique(info)
  
  t <- info$MJD
  
  y <- info$Rcat
  
  err <- info$Rcaterr
  
  n <- length(t)
  
  timespan <- max(t) - min(t)
  
  ts <- cbind(t, y, err)
  
  n_steps <- as.integer(ofac*n*(1/plow - 1/phigh))
  f <- seq(from = 1/phigh, to = 1/plow, length.out = n_steps)
  f <- sort(f)
  
  num <- seq(obsstart, n, by = obsstep)
  
  rownum <- length(num)
  
  logp <- matrix(0, nrow = rownum, ncol = n_steps)
  
  k <- 1
  
  for (i in num){
    
    ans <- BGLS2(ts[1:i,], plow = plow, phigh = phigh, ofac = ofac, dt = timespan, lent = n, spur = F, jit = jit)
    
    logp[k,] <- ans$p
    
    k <- k + 1
    
  }
  
  image(x = sort(1/f), y = num, t(logp)[nrow(t(logp)):1,], axes = T, col = grey(seq(0, 1, length = 256)))
  
  Sys.time() - starttime
  
}



SBGLS2 <- function(ts, plow = 0.5, phigh = 100, ofac = 1, jit = 0, obsstart = 50, obsstep = 1){
  
  starttime <- Sys.time()
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  n <- length(t)
  
  timespan <- max(t) - min(t)
  
  ts <- cbind(t, y, err)
  
  n_steps <- as.integer(ofac*n*(1/plow - 1/phigh))
  f <- seq(from = 1/phigh, to = 1/plow, length.out = n_steps)
  f <- sort(f)
  
  num <- seq(obsstart, n, by = obsstep)
  
  rownum <- length(num)
  
  logp <- matrix(0, nrow = rownum, ncol = n_steps)
  
  k <- 1
  
  for (i in num){
    
    if (i %% 10 == 0){print(i)}
    
    ans <- BGLS2(ts[1:i,], plow = plow, phigh = phigh, ofac = ofac, dt = timespan, lent = n, spur = F, jit = jit)
    
    logp[k,] <- ans$p
    
    k <- k + 1
    
  }
  
  image(x = sort(1/f), y = num, t(logp)[nrow(t(logp)):1,], axes = T, col = grey(seq(0, 1, length = 256)))
  
  Sys.time() - starttime
  
}