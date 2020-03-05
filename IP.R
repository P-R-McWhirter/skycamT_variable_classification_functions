InPot <- function(x) {
  
  datanames <- colnames(x)
  t <- x[, 1]
  y <- x[, 2]
  err <- x[, 3]
  
  n <- length(y)
  
  sigy <- median(err)
  
  IP <- 0
  
  dif <- outer(y, y, '-')
  
  Gx <- (1 / (sqrt(2 * pi) * sigy)) * exp(-0.5 * dif^2 / sigy^2)
  
  IP <- mean(mean(Gx))
  
  IP
  
}