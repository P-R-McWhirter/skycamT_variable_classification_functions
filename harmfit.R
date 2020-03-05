harm8fit <- function(ts, f, lambda = 1e-2){
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  n <- length(t)
  
  Z <- as.matrix(cbind(1, t, sin(2*pi*f*t), cos(2*pi*f*t), sin(4*pi*f*t), cos(4*pi*f*t), sin(6*pi*f*t), cos(6*pi*f*t), sin(8*pi*f*t), cos(8*pi*f*t), sin(10*pi*f*t), cos(10*pi*f*t), sin(12*pi*f*t), cos(12*pi*f*t), sin(14*pi*f*t), cos(14*pi*f*t), sin(16*pi*f*t), cos(16*pi*f*t)))
  
  W <- as.matrix(diag(err^-1))
  
  M <- t(Z)%*%W%*%Z + diag(c(1, 1, 1, 1, 16, 16, 81, 81, 256, 256, 625, 625, 1296, 1296, 2401, 2401, 4096, 4096)*n)*lambda
  
  coeff <- solve(M)%*%t(Z)%*%W%*%y
  
  coeff[is.na(coeff)] <- 0
  
  coeff
  
}


harm4fit <- function(ts, f, lambda = 1e-2){
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  n <- length(t)
  
  Z <- as.matrix(cbind(1, t, sin(2*pi*f*t), cos(2*pi*f*t), sin(4*pi*f*t), cos(4*pi*f*t), sin(6*pi*f*t), cos(6*pi*f*t), sin(8*pi*f*t), cos(8*pi*f*t)))
  
  W <- as.matrix(diag(err^-1))
  
  M <- t(Z)%*%W%*%Z + diag(c(1, 1, 1, 1, 16, 16, 81, 81, 256, 256)*n)*lambda
  
  coeff <- solve(M)%*%t(Z)%*%W%*%y
  
  coeff[is.na(coeff)] <- 0
  
  coeff
  
}


harm2fit <- function(ts, f, lambda = 1e-2){
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  n <- length(t)
  
  Z <- as.matrix(cbind(1, t, sin(2*pi*f*t), cos(2*pi*f*t), sin(4*pi*f*t), cos(4*pi*f*t)))
  
  W <- as.matrix(diag(err^-1))
  
  M <- t(Z)%*%W%*%Z + diag(c(1, 1, 1, 1, 1, 1)*n)*lambda
  
  coeff <- solve(M)%*%t(Z)%*%W%*%y
  
  coeff[is.na(coeff)] <- 0
  
  coeff
  
}

harm8fit_lin <- function(ts, f, lambda = 1e-2){
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  n <- length(t)
  
  Z <- as.matrix(cbind(1, sin(2*pi*f*t), cos(2*pi*f*t), sin(4*pi*f*t), cos(4*pi*f*t), sin(6*pi*f*t), cos(6*pi*f*t), sin(8*pi*f*t), cos(8*pi*f*t), sin(10*pi*f*t), cos(10*pi*f*t), sin(12*pi*f*t), cos(12*pi*f*t), sin(14*pi*f*t), cos(14*pi*f*t), sin(16*pi*f*t), cos(16*pi*f*t)))
  
  W <- as.matrix(diag(err^-1))
  
  M <- t(Z)%*%W%*%Z + diag(c(1, 1, 1, 16, 16, 81, 81, 256, 256, 625, 625, 1296, 1296, 2401, 2401, 4096, 4096)*n)*lambda
  
  coeff <- solve(M)%*%t(Z)%*%W%*%y
  
  coeff[is.na(coeff)] <- 0
  
  coeff
  
}



harm16fit_lin <- function(ts, f, lambda = 1e-2){
  
  t <- ts[,1]
  
  y <- ts[,2]
  
  err <- ts[,3]
  
  n <- length(t)
  
  Z <- as.matrix(cbind(1, sin(2*pi*f*t), cos(2*pi*f*t), sin(4*pi*f*t), cos(4*pi*f*t), sin(6*pi*f*t), cos(6*pi*f*t), sin(8*pi*f*t), cos(8*pi*f*t), sin(10*pi*f*t), cos(10*pi*f*t), sin(12*pi*f*t), cos(12*pi*f*t), sin(14*pi*f*t), cos(14*pi*f*t), sin(16*pi*f*t), cos(16*pi*f*t), sin(18*pi*f*t), cos(18*pi*f*t), sin(20*pi*f*t), cos(20*pi*f*t), sin(22*pi*f*t), cos(22*pi*f*t), sin(24*pi*f*t), cos(24*pi*f*t), sin(26*pi*f*t), cos(26*pi*f*t), sin(28*pi*f*t), cos(28*pi*f*t), sin(30*pi*f*t), cos(30*pi*f*t), sin(32*pi*f*t), cos(32*pi*f*t)))
  
  W <- as.matrix(diag(err^-1))
  
  M <- t(Z)%*%W%*%Z + diag(c(1, 1, 1, 16, 16, 81, 81, 256, 256, 625, 625, 1296, 1296, 2401, 2401, 4096, 4096, 6561, 6561, 10000, 10000, 14641, 14641, 20736, 20736, 28561, 28561, 38416, 38416, 50625, 50625, 65536, 65536)*n)*lambda
  
  coeff <- solve(M)%*%t(Z)%*%W%*%y
  
  coeff[is.na(coeff)] <- 0
  
  coeff
  
}