sacf <- function(data, T = -99, plot = FALSE) {
  
  magnitude <- data[,2]
  
  time <- data[,1]
  
  n = length(time)
  
  if (T == -99){
    
    deltaT <- time[2:n] - time[1:(n-1)]
    
    sorted_deltaT <- sort(deltaT)
    
    T <- sorted_deltaT[as.integer(n * 0.05) + 1]
    
  }
  
  K = 100
  
  sacf_res <- slot_autocor(magnitude, time, T, K)
  
  SAC <- sacf_res$SAC
  
  SACT <- SAC
  
  slots <- sacf_res$slots
  
  SAC2 <- SAC[slots]
  
  for (i in 1:length(SACT)){
    
    index <- i
    
    value <- SACT[i]
    
    if (length(SACT) > 0){
      
      if (value < exp(-1.0) & value != 0){
        
        autocorlen <- i
        
        break
        
      }
      else{
        
        autocorlen <- NA
        
      }
      
    }
    else{
      
      autocorlen <- NA
      
    }
    
  }
  
  while (is.na(autocorlen)){
    
    K <- K + K
    
    if (K > ((max(time) - min(time)) / T)){
      
      break
      
    }
    else{
      
      sacf_res <- slot_autocor(magnitude, time, T, K, second_round = TRUE, k1 = K/2)
      
      SAC <- sacf_res$SAC
      
      SACT <- c(SACT, SAC)
      
      slots <- sacf_res$slots
      
      SAC2 <- SAC[slots]
      
      for (i in 1:length(SACT)){
        
        index <- i
        
        value <- SACT[i]
        
        if (length(SACT) > 0){
          
          if (value < exp(-1.0) & value != 0){
            
            autocorlen <- i
            
            break
            
          }
          else{
            
            autocorlen <- NA
            
          }
          
        }
        else{
          
          autocorlen <- NA
          
        }
        
      }
      
    }
    
  }
  
  SACT <- SACT[which(SACT != 0 & SACT != 1)]
  
  SACT <- c(1, SACT)
  
  autocorlen <- min(which(SACT < exp(-1.0), arr.ind = TRUE))
  
  ans <- (autocorlen) * T
  
  SACT <- SACT[1:as.integer(autocorlen + 1)]
  
  SACT[which(SACT == Inf)] <- max(SACT[which(SACT != Inf)])
  
  if (plot == TRUE){
    
    plot((0:(length(SACT)-1)), SACT, pch=19, type="n", main="Slotted Autocorrelation Function", xlab = "Slots", ylab = "Slotted Autocorrelation")
    
    lines((0:(length(SACT)-1)), SACT, type="l")
    
  }
  
  sp.out <- list(Length = ans, ACL = autocorlen, SAC = SACT, slots = slots, T = T)
  
}