statplots <- function(featvec) {
  
  print("Plotting Box Plots for design matrix...")
  
  featvec$Type <- factor(featvec$Type)
  
  boxplot(`Mean.Magnitude`~Type, data=featvec, main="Box Plot of Mean Magnitude", xlab = "Class", ylab="Mean Magnitude")
  
  boxplot(`Mean.Variance`~Type, data=featvec, main="Box Plot of Mean Variance", xlab = "Class", ylab="Mean Variance")
  
  boxplot(`Variability.Index`~Type, data=featvec, main="Box Plot of Variability Index", xlab = "Class", ylab="Variability Index")
  
  boxplot(`Q31`~Type, data=featvec, main="Box Plot of Q31", xlab = "Class", ylab="Q31")
  
  boxplot(`Skewness`~Type, data=featvec, main="Box Plot of Skewness", xlab = "Class", ylab="Skewness")
  
  boxplot(`Small.Kurtosis`~Type, data=featvec, main="Box Plot of Small Kurtosis", xlab = "Class", ylab="Small Kurtosis")
  
  boxplot(`StD`~Type, data=featvec, main="Box Plot of StD", xlab = "Class", ylab="StD")
  
  boxplot(`Rcs`~Type, data=featvec, main="Box Plot of Rcs", xlab = "Class", ylab="Rcs")
  
  boxplot(`Beyond.1.StD`~Type, data=featvec, main="Box Plot of Beyond 1 StD", xlab = "Class", ylab="Beyond 1 StD")
  
  boxplot(`Stetson.Kurtosis.Measure`~Type, data=featvec, main="Box Plot of Stetson Kurtosis Measure", xlab = "Class", ylab="Stetson Kurtosis Measure")
  
  boxplot(`Stetson.K.AC`~Type, data=featvec, main="Box Plot of Stetson-K AC", xlab = "Class", ylab="Stetson-K AC")
  
  boxplot(`Maximum.Slope`~Type, data=featvec, main="Box Plot of Maximum Slope", xlab = "Class", ylab="Maximum Slope")
  
  boxplot(`Amplitude`~Type, data=featvec, main="Box Plot of Amplitude", xlab = "Class", ylab="Amplitude")
  
  boxplot(`Median.Absolute.Deviation`~Type, data=featvec, main="Box Plot of Median Absolute Deviation", xlab = "Class", ylab="Median Absolute Deviation")
  
  boxplot(`Median.Buffer.Range.Percentage`~Type, data=featvec, main="Box Plot of Median Buffer Range Percentage", xlab = "Class", ylab="Median Buffer Range Percentage")
  
  boxplot(`Pair.Slope.Trend`~Type, data=featvec, main="Box Plot of Pair Slope Trend", xlab = "Class", ylab="Pair Slope Trend")
  
  boxplot(`Flux.Percentile.Ratio.mid20`~Type, data=featvec, main="Box Plot of Flux Percentile Ratio mid20", xlab = "Class", ylab="Flux Percentile Ratio mid20")
  
  boxplot(`Flux Percentile Ratio mid35`~Type, data=featvec, main="Box Plot of Flux Percentile Ratio mid35", xlab = "Class", ylab="Flux Percentile Ratio mid35")
  
  boxplot(`Flux Percentile Ratio mid50`~Type, data=featvec, main="Box Plot of Flux Percentile Ratio mid50", xlab = "Class", ylab="Flux Percentile Ratio mid50")
  
  boxplot(`Flux Percentile Ratio mid65`~Type, data=featvec, main="Box Plot of Flux Percentile Ratio mid65", xlab = "Class", ylab="Flux Percentile Ratio mid65")
  
  boxplot(`Flux Percentile Ratio mid80`~Type, data=featvec, main="Box Plot of Flux Percentile Ratio mid80", xlab = "Class", ylab="Flux Percentile Ratio mid80")
  
  boxplot(`Percent Amplitude`~Type, data=featvec, main="Box Plot of Percent Amplitude", xlab = "Class", ylab="Percent Amplitude")
  
  boxplot(`Percent Difference Flux Percentile`~Type, data=featvec, main="Box Plot of Percent Difference Flux Percentile", xlab = "Class", ylab="Percent Difference Flux Percentile")
  
  boxplot(`Anderson-Darling Statistic`~Type, data=featvec, main="Box Plot of Anderson-Darling Statistic", xlab = "Class", ylab="Anderson-Darling Statistic")
  
  boxplot(`Autocorrelation Length`~Type, data=featvec, main="Box Plot of Autocorrelation Length", xlab = "Class", ylab="Autocorrelation Length")
  
  boxplot(`Slotted Autocorrelation Length`~Type, data=featvec, main="Box Plot of Slotted Autocorrelation Length", xlab = "Class", ylab="Slotted Autocorrelation Length")
  
  if (length(featvec) > 32){
    
    boxplot(`Slope of Linear Trend`~Type, data=featvec, main="Box Plot of Slope of Linear Trend", xlab = "Class", ylab="Slope of Linear Trend")
    
    boxplot(`#1 Freq`~Type, data=featvec, main="Box Plot of #1 Freq", xlab = "Class", ylab="#1 Freq")
    
    boxplot(`#1 p-value`~Type, data=featvec, main="Box Plot of #1 p-value", xlab = "Class", ylab="#1 p-value")
    
    boxplot(`#2 Freq`~Type, data=featvec, main="Box Plot of #2 Freq", xlab = "Class", ylab="#2 Freq")
    
    boxplot(`#2 p-value`~Type, data=featvec, main="Box Plot of #2 p-value", xlab = "Class", ylab="#2 p-value")
    
    boxplot(`#3 Freq`~Type, data=featvec, main="Box Plot of #3 Freq", xlab = "Class", ylab="#3 Freq")
    
    boxplot(`#3 p-value`~Type, data=featvec, main="Box Plot of #3 p-value", xlab = "Class", ylab="#3 p-value")
    
    boxplot(`A11`~Type, data=featvec, main="Box Plot of A11", xlab = "Class", ylab="A11")
    
    boxplot(`A12`~Type, data=featvec, main="Box Plot of A12", xlab = "Class", ylab="A12")
    
    boxplot(`A13`~Type, data=featvec, main="Box Plot of A13", xlab = "Class", ylab="A13")
    
    boxplot(`A14`~Type, data=featvec, main="Box Plot of A14", xlab = "Class", ylab="A14")
    
    boxplot(`A21`~Type, data=featvec, main="Box Plot of A21", xlab = "Class", ylab="A21")
    
    boxplot(`A22`~Type, data=featvec, main="Box Plot of A22", xlab = "Class", ylab="A22")
    
    boxplot(`A23`~Type, data=featvec, main="Box Plot of A23", xlab = "Class", ylab="A23")
    
    boxplot(`A24`~Type, data=featvec, main="Box Plot of A24", xlab = "Class", ylab="A24")
    
    boxplot(`A31`~Type, data=featvec, main="Box Plot of A31", xlab = "Class", ylab="A31")
    
    boxplot(`A32`~Type, data=featvec, main="Box Plot of A32", xlab = "Class", ylab="A32")
    
    boxplot(`A33`~Type, data=featvec, main="Box Plot of A33", xlab = "Class", ylab="A33")
    
    boxplot(`A34`~Type, data=featvec, main="Box Plot of A34", xlab = "Class", ylab="A34")
    
    boxplot(`PH12`~Type, data=featvec, main="Box Plot of PH12", xlab = "Class", ylab="PH12")
    
    boxplot(`PH13`~Type, data=featvec, main="Box Plot of PH13", xlab = "Class", ylab="PH13")
    
    boxplot(`PH14`~Type, data=featvec, main="Box Plot of PH14", xlab = "Class", ylab="PH14")
    
    boxplot(`PH21`~Type, data=featvec, main="Box Plot of PH21", xlab = "Class", ylab="PH21")
    
    boxplot(`PH22`~Type, data=featvec, main="Box Plot of PH22", xlab = "Class", ylab="PH22")
    
    boxplot(`PH23`~Type, data=featvec, main="Box Plot of PH23", xlab = "Class", ylab="PH23")
    
    boxplot(`PH24`~Type, data=featvec, main="Box Plot of PH24", xlab = "Class", ylab="PH24")
    
    boxplot(`PH31`~Type, data=featvec, main="Box Plot of PH31", xlab = "Class", ylab="PH31")
    
    boxplot(`PH32`~Type, data=featvec, main="Box Plot of PH32", xlab = "Class", ylab="PH32")
    
    boxplot(`PH33`~Type, data=featvec, main="Box Plot of PH33", xlab = "Class", ylab="PH33")
    
    boxplot(`PH34`~Type, data=featvec, main="Box Plot of PH34", xlab = "Class", ylab="PH34")
    
    boxplot(`Variance Ratio`~Type, data=featvec, main="Box Plot of Variance Ratio", xlab = "Class", ylab="Variance Ratio")
    
    boxplot(`Folded Variability Index`~Type, data=featvec, main="Box Plot of Folded Variability Index", xlab = "Class", ylab="Folded Variability Index")
    
    boxplot(`Psi-cs`~Type, data=featvec, main="Box Plot of Psi-cs", xlab = "Class", ylab="Psi-cs")
    
  }
  
  print("Plotting Complete.")
  
}




statplots_neuro <- function(featvec) {
  
  print("Plotting Box Plots for design matrix...")
  
  featvec$Type <- factor(featvec$Type)
  
  boxplot(`Mean.Magnitude`~Type, data=featvec, main="Box Plot of Mean Magnitude", xlab = "Class", ylab="Mean Magnitude")
  
  boxplot(`Q31`~Type, data=featvec, main="Box Plot of Q31", xlab = "Class", ylab="Q31")
  
  boxplot(`Rcs`~Type, data=featvec, main="Box Plot of Rcs", xlab = "Class", ylab="Rcs")
  
  boxplot(`Stetson.Kurtosis.Measure`~Type, data=featvec, main="Box Plot of Stetson Kurtosis Measure", xlab = "Class", ylab="Stetson Kurtosis Measure")
  
  boxplot(`Autocorrelation.Length`~Type, data=featvec, main="Box Plot of Autocorrelation Length", xlab = "Class", ylab="Autocorrelation Length")
  
  boxplot(`X.1.Freq`~Type, data=featvec, main="Box Plot of #1 Freq", xlab = "Class", ylab="#1 Freq")
  
  boxplot(`A11`~Type, data=featvec, main="Box Plot of A11", xlab = "Class", ylab="A11")
    
  boxplot(`Psi.cs`~Type, data=featvec, main="Box Plot of Psi-cs", xlab = "Class", ylab="Psi-cs")
  
  boxplot(`Renyi.s.Quadratic.Entropy`~Type, data=featvec, main="Box Plot of the Renyi Quadratic Entropy", xlab = "Class", ylab="Renyi Quadratic Entropy")
  
  print("Plotting Complete.")
  
}