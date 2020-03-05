complete <- function(directory, id = 1:332) {
  
  files_full <- list.files(directory, full.names=TRUE)
  
  dat <- data.frame()
  
  nobs <- c()
  
  for (i in id) {
    
    hold <- read.csv(files_full[i])
    
    good <- complete.cases(hold)
    
    hogd <- hold[good,1]
    
    hosi <- length(hogd)
    
    nobs <- c(nobs, hosi)
    
  }
  
  res <- data.frame(id, nobs)
  
  res
}