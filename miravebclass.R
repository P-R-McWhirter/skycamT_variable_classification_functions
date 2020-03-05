miravebclass <- function(x) {
  
  y <- x[,1]
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  p <- matrix(0, length(y), 10)
  
  for (i in 1:length(y)){
    
    print(paste("Object number ", i, " commencing period feature vector generation...", sep=""))
    
    p[i,] <- lcg(x[i,2], x[i,3], 0.01, 0.01, 1000, 100, 2, rang = 2, rand = FALSE, ploterr = TRUE, aalldb = FALSE, returnp = TRUE)
    
  }
  
  p
  
}