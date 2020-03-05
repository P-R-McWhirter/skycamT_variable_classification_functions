lsp_amp <- function (x, times = NULL, from = NULL, to = NULL, type = c("frequency", 
                                                            "period"), ofac = 1, alpha = 0.01, plot = TRUE, ...) 
{
  type <- match.arg(type)
  if (ofac != floor(ofac)) {
    ofac <- floor(ofac)
    warning("ofac coerced to integer")
  }
  if (ofac < 1) {
    ofac <- 1
    warning("ofac must be integer >=1. Set to 1")
  }
  if (!is.null(times)) {
    if (!is.vector(times)) 
      stop("no multivariate methods available")
    if (length(x) != length(times)) 
      stop("Length of data and times vector must be equal")
    names <- c(deparse(substitute(times)), deparse(substitute(x)))
  }
  if (is.null(times) && is.null(ncol(x))) {
    names <- c("Time", deparse(substitute(x)))
    times <- 1:length(x)
  }
  if (is.matrix(x) || is.data.frame(x)) {
    if (ncol(x) > 2) 
      stop("no multivariate methods available")
    if (ncol(x) == 2) {
      names <- colnames(x)
      times <- x[, 1]
      x <- x[, 2]
    }
  }
  times <- times[!is.na(x)]
  x <- x[!is.na(x)]
  nobs <- length(x)
  if (nobs < 2) 
    stop("time series must have at least two observations")
  times <- as.numeric(times)
  start <- min(times)
  end <- max(times)
  av.int <- mean(diff(times))
  o <- order(times)
  times <- times[o]
  x <- x[o]
  y <- cbind(times, x)
  colnames(y) <- names
  datanames <- colnames(y)
  t <- y[, 1]
  y <- y[, 2]
  n <- length(y)
  tspan <- t[n] - t[1]
  fr.d <- 1/tspan
  step <- 1/(tspan * ofac)
  if (type == "period") {
    hold <- from
    from <- to
    to <- hold
    if (!is.null(from)) 
      from <- 1/from
    if (!is.null(to)) 
      to <- 1/to
  }
  if (is.null(to)) {
    f.max <- floor(0.5 * n * ofac) * step
  }
  else {
    f.max <- to
  }
  freq <- seq(fr.d, f.max, by = step)
  if (!is.null(from)) 
    freq <- freq[freq >= from]
  n.out <- length(freq)
  if (n.out == 0) 
    stop("erroneous frequency range specified ")
  x <- t * 2 * pi
  y <- y - mean(y)
  norm <- 1/(2 * var(y))
  w <- 2 * pi * freq
  PN <- rep(0, n.out)
  for (i in 1:n.out) {
    wi <- w[i]
    tau <- 0.5 * atan2(sum(sin(wi * t)), sum(cos(wi * t)))/wi
    arg <- wi * (t - tau)
    cs <- cos(arg)
    sn <- sin(arg)
    A <- (sum(y * cs))^2
    B <- sum(cs * cs)
    C <- (sum(y * sn))^2
    D <- sum(sn * sn)
    PN[i] <- A/B + C/D
  }
  PN <- norm * PN
  PN.max <- max(PN)
  peak.freq <- freq[PN == PN.max]
  if (type == "period") 
    peak.at <- c(1/peak.freq, peak.freq)
  else peak.at <- c(peak.freq, 1/peak.freq)
  scanned <- if (type == "frequency") 
    freq
  else 1/freq
  if (type == "period") {
    scanned <- scanned[n.out:1]
    PN <- PN[n.out:1]
  }
  effm <- 2 * n.out/ofac
  level <- -log(1 - (1 - alpha)^(1/effm))
  exPN <- exp(-PN.max)
  p <- effm * exPN
  if (p > 0.01) 
    p <- 1 - (1 - exPN)^effm
  PN <- sqrt(PN)
  level <- sqrt(level)
  PN.max = sqrt(PN.max)
  #scanned <- log10(scanned)
  sp.out <- list(scanned = scanned, power = PN, data = datanames, 
                 n = n, type = type, ofac = ofac, n.out = n.out, alpha = alpha, 
                 sig.level = level, peak = PN.max, peak.at = peak.at, 
                 p.value = p)
  class(sp.out) <- "lsp"
  if (plot) {
    plot(sp.out, ylab = "sqrt(normalised power)", ...)
    return(invisible(sp.out))
  }
  else return(sp.out)
}



spurper_amp <- function(x, from = 0, to = 20, ofac = 1, name = "", plot = TRUE) {
  
  from <- 0
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  
  fmin <- from
  fmax <- to
  
  n <- length(t)
  
  ttot <- max(t) - min(t)
  
  df <- 1 / (ttot * ofac)
  
  f <- seq(fmin, fmax, by = df)
  
  m <- length(f)
  
  p <- rep(0, m)
  
  print("Calculating Spurious Frequencies...")
  
  for (j in 1:m){
    
    p[j] <- abs(sum(exp(1i * 2.0 * pi * f[j] * t)))^2 / n
    
  }
  
  p <- sqrt(p)
  
  peak.at = c(f[which.max(p)], 1/f[which.max(p)])
  
  if (plot == TRUE){
    
    plot(c(rev(-f), f), c(rev(p), p), pch=19, type="n", main=paste("Window Function for ", name, "", sep=""), xlab = "Frequency (Cycles/days)", ylab = "Spectral Power", ylim = c(min(p), max(p)))
    
    lines(c(rev(-f), f), c(rev(p), p), type="l")
    
  }
  
  sp.out <- list(scanned = f, power = p, data = datanames, n = n, ofac = ofac, peak = max(p), peak.at = peak.at)
  
}





spurper_amp_offset <- function(x, to = 20, offset = 0, xlim = NULL, ofac = 1, name = "", plot = TRUE) {
  
  from <- 0
  xb <- x
  datanames <- colnames(x)
  t <- x[, 1]
  tb <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  
  fmin <- from
  fmax <- to
  
  if (is.null(xlim)){
    
    xlim = c(offset-to, to+offset)
    
  }
  
  n <- length(t)
  
  ttot <- max(t) - min(t)
  
  df <- 1 / (ttot * ofac)
  
  f <- seq(fmin, fmax, by = df)
  
  m <- length(f)
  
  p <- rep(0, m)
  
  print("Calculating Spurious Frequencies...")
  
  for (j in 1:m){
    
    p[j] <- abs(sum(exp(1i * 2.0 * pi * f[j] * t)))^2 / n
    
  }
  
  p <- sqrt(p)
  
  peak.at = c(f[which.max(p)], 1/f[which.max(p)])
  
  if (plot == TRUE){
    
    plot(c(rev(-f+offset), f+offset), c(rev(p), p), pch=19, type="n", main=paste("Window Function for ", name, "", sep=""), xlab = "Frequency (Cycles/days)", ylab = "Spectral Power", xlim = xlim, ylim = c(min(p), max(p)))
    
    lines(c(rev(-f+offset), f+offset), c(rev(p), p), type="l")
    
  }
  
  sp.out <- list(scanned = f, power = p, data = datanames, n = n, ofac = ofac, peak = max(p), peak.at = peak.at)
  
}