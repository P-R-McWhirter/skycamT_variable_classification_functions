fLSPw <- function(x, y, f, p0 = 0.05, iofac = 4) {
  
  library(pracma)
  
  library(stats)
  
  type = "frequency"
  alpha = 0.01
  level = 0
  
  nfreq <- 64
  
  macc <- 4
  
  f <- sort(f)
  
  fmax <- f[length(f)]
  
  fmin <- f[1]
  
  n <- length(y)
  
  ave <- mean(y)
  
  vr <- var(y)
  
  xmin <- min(x)
  
  xmax <- max(x)
  
  xdif <- xmax - xmin
  
  ofac <- max(iofac, (1/(xdif*fmin)))
  
  df <- 1 / (xdif*ofac)
  
  fc <- n / (2*xdif)
  
  hifac <- f/fc
  
  nout <- rep.int(0, length(hifac))
  
  for (k in 1:length(hifac)){
    
    if ((0.5*ofac*hifac[k]*n) < 0){
      
      nout[k] <- ceiling(0.5*ofac*hifac[k]*n)
      
    }
    else
    {
      
      nout[k] <- floor(0.5*ofac*hifac[k]*n)
      
    }
    
  }
  
  noutmax <- nout[length(nout)]
  
  nfreqt <- 2*noutmax*macc
  
  if (nfreq < nfreqt){
    
    nfreq <- 2^(nextpow2(nfreqt))
    
  }
  
  ndim <- nfreq
  
  wk1 <- rep.int(0, ndim)
  
  wk2 <- wk1
  
  fac <- ndim / (xdif * ofac)
  
  fndim <- ndim
  
  ck <- 1 + rem((x-xmin)*fac, fndim)
  
  ckk <- 1 + rem(2*(ck-1), fndim)
  
  for (j in 1:n){
    
    print(j)
    
    wk1 <- spread(wk1, (y[j] - ave), ndim, ck[j], macc)
    
    wk2 <- spread(wk2, 1, ndim, ckk[j], macc)
    
  }
  
  tmp1 <- fft(wk1[1:nfreq])
  
  tmp2 <- fft(wk2[1:nfreq])
  
  #wk1 <- fliplr(tmp1[(length(tmp1)/2+2):length(tmp1)])
  
  wk1 <- tmp1[length(tmp1):(length(tmp1)/2+2)]
  
  wk1 <- c(Re(wk1), Im(wk1))
  
  #wk2 <- fliplr(tmp2[(length(tmp2)/2+2):length(tmp2)])
  
  wk2 <- tmp2[length(tmp2):(length(tmp2)/2+2)]
  
  wk2 <- c(Re(wk2), Im(wk2))
  
  k <- seq(from = 1, to = (2*noutmax-1), by = 2)
  
  kp1 <- k + 1
  
  hypo <- sqrt(wk2[k]^2 + wk2[kp1]^2)
  
  hc2wt <- 0.5 * wk2[k] / hypo
  
  hs2wt <- 0.5 * wk2[kp1] / hypo
  
  cwt <- sqrt(0.5 + hc2wt)
  
  swt <- abs(sqrt(0.5 - hc2wt)) * sign(hs2wt)
  
  den <- 0.5 * n + hc2wt * wk2[k] + hs2wt * wk2[kp1]
  
  cterm <- (cwt * wk1[k] + swt * wk1[kp1])^2 / den
  
  sterm <- (cwt * wk1[kp1] - swt * wk1[k])^2 / (n - den)
  
  P <- (cterm + sterm) / (2 * vr)
  
  F <- (1:noutmax) * df
  
  fr <- list()
  
  Pr <- list()
  
  A <- rep.int(0, (length(f)-1))
  
  z0 <- A
  
  A0 <- A
  
  for (n in 1:(length(f)-1)){
    
    fr[[n]] <- F[(nout[n]+1):nout[n+1]]
    
    Pr[[n]] <- P[(nout[n]+1):nout[n+1]]
    
    A[n] <- df * sum(Pr[[n]])
    
    z0[n] <- log(length(fr[[n]])) - log(p0)
    
    A0[n] <- qgamma((1-p0), shape = length(fr[[n]]), scale = df)
    
  }
  
  #peak.at = c(1/fr[[which.max(Pr)]], fr[[which.max(Pr)]])
  
  p = 0
  
  sp.out <- list(scanned = fr, power = Pr, data = "empty", n = length(y), type = type, ofac = ofac, n.out = n, alpha = alpha, sig.level = level, 
                 peak = 1, peak.at = 1, p.value = p)
  
}