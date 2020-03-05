firstclass <- function(data, lim = 100) {
  
  start <- Sys.time()
  
  result <- data.frame()
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    object <- lcan(RA[i], DEC[i], 0.001, lim, ovsm = 2, run = 50000, type = "LSP", rang = 2, ploterr = FALSE, alldb = FALSE, usenyq = TRUE, geo = TRUE, nonper = TRUE, quiet = TRUE)
    
    object <- as.data.frame(object)
    
    object <- cbind(data$Name[i], data$Type[i], data$Period[i], object[1,])
    
    result <- rbind(result, object[1,])
    
  }
  
  colnames(result) <- c("Name", "Type", "Period", "USnoref", "# Observations", "Mean Magnitude", "Mean Variance", "Variability Index",
                        "Q31", "Skewness", "Small Kurtosis", "StD", "Renyi's Quadratic Entropy", "Rcs", "Beyond 1 StD",
                        "Stetson Kurtosis Measure", "Stetson-K AC", "Maximum Slope", "Amplitude", "Median Absolute Deviation",
                        "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20", "Flux Percentile Ratio mid35",
                        "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                        "Percent Amplitude", "Percent Difference Flux Percentile", "Anderson-Darling Statistic",
                        "Autocorrelation Length", "Slotted Autocorrelation Length")
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}