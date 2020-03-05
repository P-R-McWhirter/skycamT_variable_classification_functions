createTrainset <- function(varcatalog){
  
  data <- data.frame()
  
  for (i in 1:nrow(varcatalog)){
    
    print(paste("Computing object ", i, " of ", nrow(varcatalog), ".", sep = ""))
    
    object <- try(sclc(as.numeric(as.vector(varcatalog$RA))[i], as.numeric(as.vector(varcatalog$DEC))[i], radii = 0.1, lim = 100, geo = FALSE, pca = polyfit2.pca, spca = spline.pca, sigk = 4, lctrenddata = lctrenddata, seedno = 20, redo = 5, quiet = F, nonper = F), TRUE)
    
    if (class(object) != "try-error"){
      
      object <- cbind(rep(varcatalog$Name[i], nrow(object)), object)
      
      data <- rbind(data, object)
      
    }
    
  }
  
  colnames(data) <- c("Name", "USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                        "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                        "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                        "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                        "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value", "Slope of Linear Trend",
                        "Period", "Double Period", "Period Confidence", "Daily Confidence", "Bestfit Stat", "Num Spurious", "Period Spurious Level", "Jitter Multiplier", "H1", "H2", "H3", "H4", "R21", "R31", "R41", "PH1", "PH2", "PH3", "PH4", "PH21", "PH31", "PH41", "Variance Ratio", "mp10", "mp90", "Psi-cs",
                        "Psi-Eta", "P Binned Ratio", "P Goodness of Fit", "P Int StD", "P Int Skewness", "P Int Kurtosis", "P Int Small Kurtosis", "P Int Amplitude",
                        "P Int Beyond 1 StD", "P Int cs", paste("P PCA", 1:10, sep=""), "P2 Binned Ratio", "P2 Goodness of Fit", "P2 Int StD", "P2 Int Skewness", "P2 Int Kurtosis", "P2 Int Small Kurtosis", "P2 Int Amplitude",
                        "P2 Int Beyond 1 StD", "P2 Int cs", paste("P2 PCA", 1:10, sep=""), "Eclipse Max Delta Mags", "Eclipse Min Delta Mags", "Eclipse Phase Ratio", "Reference Phase", "Residual Normality", "Residual Raw Scatter", "Squared Diffs over Variance", "Period Double Ratio")
  
  colnames(data) <- make.names(colnames(data), unique = T)
  
  data
  
}




createTrainset_nonper <- function(varcatalog){
  
  data <- data.frame()
  
  for (i in 1:nrow(varcatalog)){
    
    print(paste("Computing object ", i, " of ", nrow(varcatalog), ".", sep = ""))
    
    object <- try(sclc(as.numeric(as.vector(varcatalog$RA))[i], as.numeric(as.vector(varcatalog$DEC))[i], radii = 0.1, lim = 100, geo = FALSE, pca = polyfit2.pca, spca = spline.pca, sigk = 4, lctrenddata = lctrenddata, seedno = 20, redo = 5, quiet = F, nonper = T), TRUE)
    
    if (class(object) != "try-error"){
      
      object <- cbind(rep(varcatalog$Name[i], nrow(object)), object)
      
      data <- rbind(data, object)
      
    }
    
  }
  
  colnames(data) <- c("USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                      "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                      "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                      "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                      "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value")
  
  colnames(data) <- make.names(colnames(data), unique = T)
  
  data
  
}