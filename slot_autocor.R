slot_autocor <- function(data, time, t = 4, k = 1, second_round = FALSE, k1 = 100) {
  
  slots <- matrix(0, nrow = k, ncol = 1)
  
  i <- as.integer(1)
  
  # make time start from 0
  time <- time - min(time)
  
  # subtract mean from mag values
  data <- data - mean(data)
  
  prod <- matrix(0, nrow = k, ncol = 1)
  
  pairs <- outer(time, time, FUN = "-")
  
  pairs[lower.tri(pairs)] <- 10000000
  
  ks <- (floor(abs(pairs) / t + 0.5))
  
  #We calculate the slotted autocorrelation for k=0 seperately
  idx <- which(ks == 0, arr.ind = TRUE)
  
  prod[1] <- ((sum(data ^ 2) + sum(data[idx[,1]] * data[idx[,2]])) / (length(idx[,1]) + length(data)))
  
  slots[1] <- 0
  
  # We calculate it dor the rest of the ks
  if (second_round == FALSE){
    
    for (z in seq.int(from = 1, to = k, by = 1)){
      
      idx <- which(ks == z, arr.ind = TRUE)
      
      if (length(idx[,1]) != 0){
        
        prod[z+1] <- sum(data[idx[,1]] * data[idx[,2]]) / (length(idx[,1]))
        
        slots[i+1] <- z
        
        i <- i + 1
        
      }
      else{
        
        prod[z+1] <- Inf
        
      }
      
    }
    
  }
  else{
    
    for (z in seq.int(from = k1, to = k, by = 1)){
      
      idx <- which(ks == z, arr.ind = TRUE)
      
      if (length(idx[,1]) != 0){
        
        prod[z+1] <- sum(data[idx[,1]] * data[idx[,2]]) / (length(idx[,1]))
        
        slots[i] <- z
        
        i <- i + 1
        
      }
      else{
        
        prod[z+1] <- Inf
        
      }
      
    }
    
    if (length(prod[which(prod != 0)]) > 0){
      
      prod <- prod[1 : max(which(prod != 0))]
      
    }
    
  }
  
  if (length(slots[which(slots != 0)]) > 0){
    
    slots <- slots[1 : max(which(slots != 0))]
    
  }
  
  sp.out <- list(SAC = prod / prod[1], slots = as.integer(as.vector(slots)))
  
}