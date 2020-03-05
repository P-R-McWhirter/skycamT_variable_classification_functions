loadtimepos <- function(data){
  
  closeAllConnections()
  
  library(RODBC)
  library(lomb)
  library(nortest)
  library(e1071)
  
  start <- Sys.time()
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  result <- matrix(NA, nrow = nrow(data), ncol = max(data$entries))
  
  for (i in 1:nrow(data)){
    
    print(paste("Computing object number ", i, ": ", data$usnoref[i], ". Number of entries: ", data$entries[i], ".", sep = ""))
  
    timepos <- sqlQuery(channel, paste("SELECT MJD FROM obsdat100split WHERE usnoref = '", 
                                    data$usnoref[i], "'", sep=""))
    
    timepos <- as.vector(as.matrix(timepos))
    
    result[i,1:length(timepos)] <- timepos
  
  }
  
  close(channel)
  
  print(Sys.time() - start)
  
  gc()
  
  result
  
}