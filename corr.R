corr <- function(directory, threshold = 0){
  
  files_full <- list.files(directory, full.names=TRUE)
  
  data <- complete("specdata")
  
  gdid <- c()
  
  for (i in 1:length(data[,1])){
    
    if (data[i,2] >= threshold){
      
      gdid <- c(gdid, data[i,1])
      
    }
    
  }
  
  res <- c()
  
  for (i in gdid) {
    
    dat <- read.csv(files_full[i])
    
    nit_set <- dat[,"nitrate"]
    
    sul_set <- dat[,"sulfate"]
    
    good <- complete.cases(dat)
    
    nit_test <- nit_set[good]
    
    sul_test <- sul_set[good]
    
    res <- c(res, cor(nit_test, sul_test))
    
  }
  
  res
  
}