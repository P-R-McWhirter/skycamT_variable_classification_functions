newclass <- function(data, lim = 100, type = "LSP", usenyq = FALSE, trueper = FALSE, quiet = TRUE) {
  
  start <- Sys.time()
  
  if (type == "LSP"){
    
    typename <- "Lomb-Scargle Periodogram"
    
  }
  else if (type == "GLS"){
    
    typename <- "Generalised Lomb-Scargle Periodogram"
    
  }
  else if (type == "VRP"){
    
    typename <- "Variance Ratio Periodogram"
    
  }
  else if (type == "SLLK"){
    
    typename <- "String Length Lafler-Kinman Periodogram"
    
  }
  else if (type == "CKP"){
    
    typename <- "Correntropy Kernelised Periodogram"
    
  }
  else
  {
    
    stop(paste("'", type, "' is an invalid Periodogram. Program will now exit.", sep=""))
    
  }
  
  print(paste("Periodogram type '", type, "' selected.", sep=""))
  
  print(paste("Periodic Feature Extraction will use the ", typename, ".", sep=""))
  
  result <- data.frame()
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  Period <- as.numeric(as.vector(data$Period))
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    if (trueper == FALSE){
    
      object <- lcxpl(RA[i], DEC[i], 0.001, lim, 4, run = 50000, type = type, rang = 2, ploterr = TRUE, alldb = FALSE, usenyq = usenyq, geo = TRUE, nonper = FALSE, trueper = FALSE, actper = 1.0, quiet = quiet)
    
    }
    else{
      
      object <- lcxpl(RA[i], DEC[i], 0.001, lim, 4, run = 50000, type = type, rang = 2, ploterr = TRUE, alldb = FALSE, usenyq = usenyq, geo = TRUE, nonper = FALSE, trueper = TRUE, actper = Period[i], quiet = quiet)
      
    }
    object <- as.data.frame(object)
    
    object <- cbind(data$Name[i], data$Type[i], data$Period[i], object[1,])
    
    result <- rbind(result, object[1,])
    
  }
  
  colnames(result) <- c("Name", "Type", "AAVSO Period", "USnoref", "LS.Period", "Psi-eta", "Psi-cs",
                        "R21", "R31", "PH21", "PH31", "Skewness", "Kurtosis", "Stetson K", "Q31", "A", "H1",
                        "W", "mp10", "mp90")
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}