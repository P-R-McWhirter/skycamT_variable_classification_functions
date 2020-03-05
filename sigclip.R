sigclip <- function(y, sig = 3, tol = 0.01){
  
  done <- FALSE
  
  n <- length(y)
  
  medy <- median(y)
  
  sdy <- sd(y)
  
  i <- 1
  
  while (done == FALSE){
    
    oldsdy <- sdy
    
    oldmedy <- medy
    
    testy <- y[which(y <= (medy + (sig * sdy)) & y >= (medy - (sig * sdy)))]
    
    nn <- length(testy)
    
    result <- which(y <= (medy + (sig * sdy)) & y >= (medy - (sig * sdy)), arr.ind = TRUE)
    
    medy <- median(testy)
    
    sdy <- sd(testy)
    
    if (((oldsdy - sdy) / sdy) <= tol){
      
      done <- TRUE
      
    }
    else{
      
      i <- i + 1
      
    }
    
  }
  
  result
  
}