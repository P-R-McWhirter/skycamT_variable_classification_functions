featgen <- function(data, radii = 0.001, lim = 100, sigk = 2.5, lctrenddata, seed = 20, redo = 5, nonper = FALSE, quiet = TRUE) {
  
  library(randomForest)
  
  start <- Sys.time()
  
  result <- data.frame()
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  Period <- as.numeric(as.vector(data$Period))
  
  print(paste("Using Seed of ", seed, ". Varying this can effect the fine tuning", sep = ""))
  
  print("of the period estimation and the Polyfit algorithm.")
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    object <- try(sclc(RA[i], DEC[i], radii = radii, lim, geo = TRUE, pca = polyfit2.pca, spca = spline.pca, sigk = sigk, lctrenddata = lctrenddata, seedno = seed, redo = redo, quiet = quiet, nonper = nonper), TRUE)
    
    if (class(object) == "try-error"){
      
      if (nonper == FALSE){
        
        object <- as.data.frame(matrix(NA, nrow = 1, ncol = 115))
        
        colnames(object, do.NULL = FALSE)
        
        colnames(object) <- c("USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                              "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                              "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                              "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                              "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value", "Slope of Linear Trend",
                              "Period", "Double Period", "Period Confidence", "Daily Confidence", "Bestfit Stat", "Num Spurious", "Period Spurious Level", "Jitter Multiplier", "H1", "H2", "H3", "H4", "R21", "R31", "R41", "PH1", "PH2", "PH3", "PH4", "PH21", "PH31", "PH41", "Variance Ratio", "mp10", "mp90", "Psi-cs",
                              "Psi-Eta", "P Binned Ratio", "P Goodness of Fit", "P Int StD", "P Int Skewness", "P Int Kurtosis", "P Int Small Kurtosis", "P Int Amplitude",
                              "P Int Beyond 1 StD", "P Int cs", paste("P PCA", 1:10, sep=""), "P2 Binned Ratio", "P2 Goodness of Fit", "P2 Int StD", "P2 Int Skewness", "P2 Int Kurtosis", "P2 Int Small Kurtosis", "P2 Int Amplitude",
                              "P2 Int Beyond 1 StD", "P2 Int cs", paste("P2 PCA", 1:10, sep=""), "Eclipse Max Delta Mags", "Eclipse Min Delta Mags", "Eclipse Phase Ratio", "Reference Phase", "Residual Normality", "Residual Raw Scatter", "Squared Diffs over Variance", "Period Double Ratio")
        
        colnames(object) <- make.names(colnames(object), unique = T)
        
      }
      else{
        
        object <- as.data.frame(matrix(NA, nrow = 1, ncol = 41))
        
        colnames(object, do.NULL = FALSE)
        
        colnames(object) <- c("USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                              "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                              "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                              "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                              "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value")
        
        colnames(object) <- make.names(colnames(object), unique = T)
        
      }
      
    }
    
    object <- as.data.frame(object)
    
    object <- cbind(data$Name[i], data$Type[i], data$Period[i], object[1,])
    
    result <- rbind(result, object[1,])
    
  }
  
  if (nonper == FALSE){
    
    colnames(result) <- c("Name", "Type", "AAVSO Period", "USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                          "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                          "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                          "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                          "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value", "Slope of Linear Trend",
                          "Period", "Double Period", "Period Confidence", "Daily Confidence", "Bestfit Stat", "Num Spurious", "Period.Spurious.Level", "Jitter Multiplier", "H1", "H2", "H3", "H4", "R21", "R31", "R41", "PH1", "PH2", "PH3", "PH4", "PH21", "PH31", "PH41", "Variance Ratio", "mp10", "mp90", "Psi-cs",
                          "Psi-Eta", "P Binned Ratio", "P Goodness of Fit", "P Int StD", "P Int Skewness", "P Int Kurtosis", "P Int Small Kurtosis", "P Int Amplitude",
                          "P Int Beyond 1 StD", "P Int cs", paste("P PCA", 1:10, sep=""), "P2 Binned Ratio", "P2 Goodness of Fit", "P2 Int StD", "P2 Int Skewness", "P2 Int Kurtosis", "P2 Int Small Kurtosis", "P2 Int Amplitude",
                          "P2 Int Beyond 1 StD", "P2 Int cs", paste("P2 PCA", 1:10, sep=""), "Eclipse Max Delta Mags", "Eclipse Min Delta Mags", "Eclipse Phase Ratio", "Reference Phase", "Residual Normality", "Residual Raw Scatter", "Squared Diffs over Variance", "Period Double Ratio")
    
  }
  else{
    
    colnames(result) <- c("Name", "Type", "AAVSO Period", "USnoref", "Right Ascension", "Declination", "# Observations", "Observation Window", "Match Distance", "Mean Magnitude", "Mean Variance", "Median Magnitude", "Robust Median Statistic", "Eta", "Inverse Von Neumann Ratio", "Modified Eta",
                          "Q31", "StD", "Skewness", "Kurtosis", "Small Kurtosis", "Rcs", "Con", "Beyond 1 StD", "Stetson I", "Stetson J", "Stetson K", "Stetson L", "Maximum Slope",
                          "Amplitude", "Median Absolute Deviation", "Median Buffer Range Percentage", "Pair Slope Trend", "Flux Percentile Ratio mid20",
                          "Flux Percentile Ratio mid35", "Flux Percentile Ratio mid50", "Flux Percentile Ratio mid65", "Flux Percentile Ratio mid80",
                          "Percent Amplitude", "Percent Difference Flux Percentile", "Autocorrelation Length", "B-R Colour", "Quick LSP Period", "Quick LSP P Value")
    
  }
  
  colnames(result) <- make.names(colnames(result), unique = T)
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}






classgen <- function(data, model, traindata, showprobs = FALSE){
  
  library(randomForest)
  
  start <- Sys.time()
  
  clsprb <- predict(model, data, type = "prob")
  
  clsnm <- colnames(clsprb)
  
  precls <- rep("", nrow(clsprb))
  
  preval <- rep(0, nrow(clsprb))
  
  tracls <- rep("", nrow(clsprb))
  
  disc <- rep(0, nrow(clsprb))
  
  for (i in 1:nrow(clsprb)){
    
    precls[i] <- clsnm[which.max(clsprb[i,])]
    
    preval[i] <- max(clsprb[i,])
    
  }
  
  prox <- predict(model, rbind(data, traindata), proximity = TRUE)$proximity
  
  for (i in 1:nrow(clsprb)){
    
    pro <- max(prox[i, (nrow(prox)-nrow(traindata)+1):nrow(prox)])
    
    disc[i] <- (1 - pro)/pro
    
  }
  
  traloc <- which(data$Name %in% traindata$Name)
  
  tracls[traloc] <- as.character(as.vector(data$Type[traloc]))
  
  result <- cbind(data[,c(1, 4, 9, 5:6)], precls, preval, disc, data[,2], tracls, data[,c(3, 46, 48, 7, 10, 79)], clsprb)
  
  colnames(result) <- c("Name", "USnoref", "Match Distance", "Right Ascension", "Declination", "Predicted Class", "Predicted Probability", "Anomaly Score", "Class", "Train Class", "AAVSO Period", "Period", "Period Confidence", "No of Epochs", "Mean Magnitude", "Amplitude", paste(clsnm, " prob", sep=""))
  
  colnames(result) <- make.names(colnames(result), unique = T)
  
  if (showprobs == FALSE){
    
    result <- result[,1:16]
    
  }
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}







vargen <- function(data, model, traindata, cutoffs = NULL, showprobs = FALSE){
  
  library(randomForest)
  
  start <- Sys.time()
  
  clsprb <- predict(model, data, type = "prob")
  
  clsnm <- colnames(clsprb)
  
  precls <- rep("", nrow(clsprb))
  
  preval <- rep(0, nrow(clsprb))
  
  tracls <- rep("", nrow(clsprb))
  
  if (is.null(cutoffs)){
  
    for (i in 1:nrow(clsprb)){
    
      precls[i] <- clsnm[which.max(clsprb[i,])]
    
      preval[i] <- max(clsprb[i,])
    
    }
    
  }
  else{
    
    cutprb <- matrix(0, nrow = nrow(clsprb), ncol = ncol(clsprb))
    
    for (i in 1:nrow(clsprb)){
      
      cutprb[i,] <- clsprb[i,] - cutoffs
      
      precls[i] <- clsnm[which.max(cutprb[i,])]
      
      preval[i] <- clsprb[i,which.max(cutprb[i,])]
      
    }
    
  }
  
  traloc <- which(data$Name %in% traindata$Name)
  
  tracls[traloc] <- as.character(as.vector(data$Type[traloc]))
  
  result <- cbind(data[,c(1, 4, 9, 5:6)], precls, preval, data[,2], tracls, data[,c(7, 10)], clsprb)
  
  colnames(result) <- c("Name", "USnoref", "Match Distance", "Right Ascension", "Declination", "Predicted Class", "Predicted Probability", "Class", "Train Class", "No of Epochs", "Mean Magnitude", paste(clsnm, " prob", sep=""))
  
  colnames(result) <- make.names(colnames(result), unique = T)
  
  if (showprobs == FALSE){
    
    result <- result[,1:11]
    
  }
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}








selcutset <- function(data, cutoff, clsnm, doprobs = FALSE){
  
  probs <- data[,17:ncol(data)]
  
  cutprobs <- probs
  
  precls <- as.character(as.vector(data$Predicted.Class))
  
  preval <- data$Predicted.Probability
  
  for (i in 1:nrow(data)){
    
    cutprobs[i,] <- probs[i,] - cutoff[i,]
    
    precls[i] <- clsnm[which.max(cutprobs[i,])]
    
    if (max(cutprobs[i,]) < 0){
      
      precls[i] <- "----"
      
    }
    
    if (doprobs == FALSE){
    
      preval[i] <- probs[i,which.max(cutprobs[i,])]
    
    }
    else{
      
      preval[i] <- max(cutprobs[i,])
      
      probs[i,] <- cutprobs[i,]
      
    }
    
  }
  
  data[,17:ncol(data)] <- probs
  
  data$Predicted.Class <- precls
  
  data$Predicted.Probability <- preval
  
  data
  
}