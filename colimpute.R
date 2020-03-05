colimpute <- function(data, maxiter = 10, ntree = 100){
  
  library(missForest)
  
  library(foreach)
  
  library(doParallel)
  
  data$B.R.Colour[which(data$B.R.Colour == 1000)] <- NA
  
  n <- nrow(data)
  
  sortedRA <- sort(data$Right.Ascension)
  
  RA1 <- sortedRA[ceiling(0.2 * n)]
  
  RA2 <- sortedRA[ceiling(0.4 * n)]
  
  RA3 <- sortedRA[ceiling(0.6 * n)]
  
  RA4 <- sortedRA[ceiling(0.8 * n)]
  
  RA <- c(min(data$Right.Ascension, na.rm = T)-1, RA1, RA2, RA3, RA4, max(data$Right.Ascension, na.rm = T)+1)
  
  sortedDEC <- sort(data$Declination)
  
  DEC1 <- sortedDEC[ceiling(0.2 * n)]
  
  DEC2 <- sortedDEC[ceiling(0.4 * n)]
  
  DEC3 <- sortedDEC[ceiling(0.6 * n)]
  
  DEC4 <- sortedDEC[ceiling(0.8 * n)]
  
  DEC <- c(min(data$Declination, na.rm = T)-1, DEC1, DEC2, DEC3, DEC4, max(data$Declination, na.rm = T)+1)
  
  print(RA)
  
  print(diff(RA))
  
  print(DEC)
  
  print(diff(DEC))
  
  registerDoParallel(cores=3)
  
  impdata <- as.data.frame(matrix(0, nrow = 0, ncol = ncol(data)))
  
  colnames(impdata) <- colnames(data)
  
  for (i in 1:5){
    
    for (j in 1:5){
    
      print(paste("Imputing Percentile Area [", i, "/5,", j, "/5].", sep = ""))
      
      objnum <- nrow(data[which(data$Right.Ascension > RA[i] & data$Right.Ascension <= RA[i+1] & data$Declination > DEC[j] & data$Declination <= DEC[j+1]),])
      
      print(paste("Number of Objects in this Area: ", objnum, ".", sep = ""))
      
      objnumnacol <- length(which(is.na(data$B.R.Colour[which(data$Right.Ascension > RA[i] & data$Right.Ascension <= RA[i+1] & data$Declination > DEC[j] & data$Declination <= DEC[j+1])])))
      
      print(paste("Number of Objects with unknown colour in this Area: ", objnumnacol, ".", sep = ""))
  
      print(paste("This is ", 100*(objnumnacol/objnum), "% of the Objects.", sep = ""))
      
      data.mis <- missForest(data[which(data$Right.Ascension > RA[i] & data$Right.Ascension <= RA[i+1] & data$Declination > DEC[j] & data$Declination <= DEC[j+1]),10:ncol(data)], maxiter = maxiter, ntree = ntree, verbose = F, parallelize = 'forests')$ximp
    
      data.mis <- cbind(data[which(data$Right.Ascension > RA[i] & data$Right.Ascension <= RA[i+1] & data$Declination > DEC[j] & data$Declination <= DEC[j+1]),1:9], data.mis)
    
      impdata <- rbind(impdata, data.mis)
    
    }
    
  }
  
  print("Applying Imputed Colours into original dataset...")
  
  impdata <- impdata[order(impdata$Name),]
  
  print("Operation Complete.")
  
  registerDoSEQ()
  
  impdata
  
}