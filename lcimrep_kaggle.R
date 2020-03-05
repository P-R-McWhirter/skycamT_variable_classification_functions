lcimrep_kaggle <- function(ts, id, band, sigk = 4, errs = FALSE, record = FALSE) {
    
    keep <- sigclip(ts[,2], sig = sigk, tol = 0.000001)
    
    ts <- ts[keep,]
    
    ts <- unique(ts)
    
    ts <- ts[order(ts[,1]),]
    
    t <- ts[,1]
    
    timdif <- max(t) - min(t)
    
    magerr <- ts[,3]
    
    dm <- matrix(0, nrow = length(ts[,1]), ncol = length(ts[,1]))
    
    dt <- matrix(0, nrow = length(ts[,1]), ncol = length(ts[,1]))
    
    for (j in 1:length(ts[,1])){
      
      for (k in 1:length(ts[,1])){
        
        if (j != k & j < k){
          
          dm[j,k] <- ts[k,2] - ts[j,2]
          
          dt[j,k] <- ts[k,1] - ts[j,1]
          
        }
        
      }
      
    }
  
  mbin <- c(-5,-3,-2.5,-2,-1.5,-1,-0.5,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,0.5,1,1.5,2,2.5,3,5)
  
  tbin <- log10(c(0, (1/145),(2/145),(3/145),(4/145),(1/25),(2/25),(3/25),1.5,2.5,3.5,4.5,5.5,7,10,20,30,60,90,120,240,600,960,2000))
  
  dm <- c(dm)
  
  dt <- c(dt)
  
  imag <- matrix(0, nrow = length(mbin)-1, ncol = length(tbin)-1)
  
  imag <- matrix(0, nrow = 64, ncol = 64)
  
  del <- cbind(dt, dm)
  
  del <- del[which(del[,1] > 0),]
  
  del[,1] <- log10(del[,1])
  
  del[,2] <- log10(abs(del[,2]))
  
  ab = matrix(c(0, -5.3, 3, 6.4), 2, 2)
  
  datbin <- bin2(del, ab = ab, nbin = c(64, 64))
  
  imag <- datbin$nc / (nrow(del)*(nrow(del)-1))
  
  #for (j in 1:length(tbin)-1){
  
  #  for (k in 1:length(mbin)-1){
  
  #    imag[j,k] <- nrow(del[which(del[,1] > tbin[j] & del[,1] <= tbin[j+1] & del[,2] > mbin[k] & del[,2] <= mbin[k+1]),])
  
  #    print(nrow(del[which(del[,1] > tbin[j] & del[,1] <= tbin[j+1] & del[,2] > mbin[k] & del[,2] <= mbin[k+1]),]))
  
  #  }
  
  #}
  
  dm <- dm[dm != 0]
  
  dt <- dt[dt != 0]
  
  #hist(dm, breaks = 50)
  
  #hist(dt, breaks = 50)
  
  #hist(dt[dt < 365], breaks = 50)
  
  #hist(dt[dt < 10], breaks = 50)
  
  #hist(dt[dt < 1], breaks = 50)
  
  #hist(dt[dt < 5/24/60], breaks = 50)
  
  out <- list(dm = dm, dt = dt, imag = imag)
  
  if (record == T){
  
    imag <- imag * 255
    
    imag <- imag[, rev(seq_len(ncol(imag)))]
  
    try(png(paste("F:/Documents/Kaggle/Abstract/", id, "-", band, ".png", sep="")), TRUE)
  
    par(mar = rep(0, 4))
  
    image(imag/max(imag), axes = FALSE, col = my.colors(256))
    
    dev.off()
    
    image(imag, axes = FALSE, col = my.colors(256))
    
  }
  
  out
  
}




lcimrepcom_kaggle <- function(ts, id, sigk = 4, record = F){
  
  out <- matrix(0, nrow = 6, ncol = 4096)
  
  bands <- c("u", "g", "r", "i", "z", "y")
  
  for (i in 1:6){
    
    ans <- try(c(lcimrep_kaggle(ts[which(ts[,4] == i-1),1:3], id = id, band = bands[i], sigk = sigk, record = F)$imag), TRUE)
    
    if (class(ans) == "try-error"){
      
      out[i,] <- rep(NA, 4096)
      
    }
    else{
      
      out[i,] <- ans
      
    }
    
    gc()
    
  }
  
  out <- as.data.frame(out)
  
  out
  
}






lcimrepall_kaggle <- function(ts, sigk = 4, record = F){
  
  out <- list()
  
  ts <- ts[,c(1,2,4,5,3)]
  
  id <- unique(ts[,1])
  
  for (i in 1:length(id)){
    
    print(paste("Processing Light Curve ", i, " out of ", length(id), ".", sep = ""))
    
    out[[i]] <- lcimrepcom_kaggle(ts[which(ts[,1] == id[i]),-1], id = id[i], record = F)
    
  }
  
  names(out) <- as.character(as.vector(id))
  
  out
  
}