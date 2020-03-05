secondclass <- function(data, radii = 0.001, lim = 100, type = "LSP", usenyq = FALSE, quiet = TRUE) {
  
  start <- Sys.time()
  
  if (type == "LSP"){
    
    typename <- "Lomb-Scargle Periodogram"
    
  }
  else if (type == "BGLS"){
    
    typename <- "Bayesian Generalised Lomb-Scargle Periodogram"
    
  }
  else if (type == "LSPL"){
    
    typename <- "Hybrid-Log Lomb-Scargle Periodogram"
    
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
  else if (type == "BKR"){
    
    typename <- "Blum-Kiefer-Rosenblatt Periodogram"
    
  }
  else if (type == "CKP"){
    
    typename <- "Correntropy Kernelised Periodogram"
    
  }
  else if (type == "IFA"){
    
    typename <- "Image Feature Analysis Periodogram"
    
  }
  else if (type == "CNN"){
    
    typename <- "Convolutional Neural Network"
    
  }
  else if (type == "FSF"){
    
    typename <- "Folded Light Curve Shape-Based Features"
    
  }
  else if (type == "GEN"){
    
    typename <- "Genetic Algorithm"
    
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
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    object <- lcan(RA[i], DEC[i], radii = radii, lim, 5, run = 50000, type = type, rang = 2, ploterr = TRUE, alldb = FALSE, usenyq = usenyq, geo = TRUE, nonper = FALSE, quiet = quiet)
    
    object <- as.data.frame(object)
    
    object <- cbind(data$Name[i], data$Type[i], data$Period[i], object[1,])
    
    result <- rbind(result, object[1,])
    
  }
  
  colnames(result) <- c("Name", "Type", "Period", "USnoref", "#.Observations", "Regression.Flag", "Slope.of.Linear.Trend", "#1.Freq", "#1.p-value", "#2.Freq",
                        "#2.p-value", "#3.Freq", "#3.p-value", "A11", "A12", "A13", "A14", "A21", "A22", "A23",
                        "A24", "A31", "A32", "A33", "A34", "PH12", "PH13", "PH14", "PH21", "PH22", "PH23", "PH24",
                        "PH31", "PH32", "PH33", "PH34", "Variance.Ratio", "Mean.Magnitude", "Mean.Variance", "Variability.Index",
                        "Folded.Variability.Index", "Q31", "Skewness", "Small.Kurtosis", "StD", "Renyi's.Quadratic.Entropy", "Rcs", "Psi-cs", "Beyond.1.StD",
                        "Stetson.Kurtosis.Measure", "Stetson-K.AC", "Maximum.Slope", "Amplitude", "Median.Absolute.Deviation",
                        "Median.Buffer.Range.Percentage", "Pair.Slope.Trend", "Flux.Percentile.Ratio.mid20", "Flux.Percentile.Ratio .mid35",
                        "Flux.Percentile.Ratio.mid50", "Flux.Percentile.Ratio.mid65", "Flux.Percentile.Ratio.mid80",
                        "Percent.Amplitude", "Percent.Difference.Flux.Percentile", "Anderson-Darling.Statistic",
                        "Autocorrelation.Length", "Slotted.Autocorrelation.Length")
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}