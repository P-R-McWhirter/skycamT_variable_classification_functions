periodmatch <- function(data, ep = 0.01) {
  
  datalen <- length(data[,1])
  
  print(paste("There are ", datalen, " samples.", sep=""))
  
  data$Type <- factor(data$Type)
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  data <- data[which(!is.na(data$Period)),]
  
  postrem <- length(data[,1])
  
  removed <- datalen - postrem
  
  print("Filtering out objects lacking AAVSO periods.")
  
  print(paste("This removes ", removed, " objects leaving ", postrem, " for future operations.", sep=""))
  
  print("Matching periods from AAVSO to those from local Periodogram...")
  
  print(paste("A period will be accepted if within +/-", ep*100, "% of the AAVSO value.", sep=""))
  
  print("Objects without an AAVSO period will be ignored.")
  
  print("--------------------------------------------")
  
  for (i in 1:length(levels(data$Type))){
    
    print(paste("Computing matches for the ", levels(data$Type)[i], " class with AAVSO period data.", sep=""))
    
    totnum <- length(data$Type[which(data$Type == levels(data$Type)[i])])
    
    print(paste("There is a total number of ", totnum, " objects."))
    
    matchnum <- length(data$Period[which((abs(data$Period - (1 / data$`X.1.Freq`)) < ep * data$Period) & (data$Type == levels(data$Type)[i]))])
    
    multnum <- length(data$Period[which(((1 / data$`X.1.Freq`) > data$Period) &
                                          (abs(((1 / data$`X.1.Freq`) / data$Period) - floor((1 / data$`X.1.Freq`) / data$Period)) < ep) & (data$Type == levels(data$Type)[i]))]) + 
                                          length(data$Period[which(((1 / data$`X.1.Freq`) < data$Period) &
                                          (abs((data$Period / (1 / data$`X.1.Freq`)) - floor(data$Period / (1 / data$`X.1.Freq`))) < ep) & (data$Type == levels(data$Type)[i]))]) - matchnum
    
    print(paste("The number of objects fulfilling the +/-", ep*100, "% requirement is ", matchnum, ".", sep=""))
    
    print(paste("The number of multiples fulfilling the +/-", ep*100, "% requirement is ", multnum, ".", sep=""))
    
    permat <- (matchnum / totnum) * 100
    
    permult <- (multnum / totnum) * 100
    
    perall <- ((matchnum + multnum) / totnum) * 100
    
    print(paste("The matches are ", permat, "% of the total number.", sep=""))
    
    print(paste("The multiples are ", permult, "% of the total number.", sep=""))
    
    print(paste("Combined this is ", perall, "% of the total number.", sep=""))
    
    print(paste("The miss rate is ", 100 - perall, "% of the total number.", sep=""))
    
    print("Next, we grant a match if the correct period was in the top 3 scanned periods.")
    
    matchnum2 <- length(data$Period[which((abs(data$Period - (1 / data$`X.2.Freq`)) < ep * data$Period) & (data$Type == levels(data$Type)[i]))])
    
    multnum2 <- length(data$Period[which(((1 / data$`X.2.Freq`) > data$Period) &
                                            (abs(((1 / data$`X.2.Freq`) / data$Period) - floor((1 / data$`X.2.Freq`) / data$Period)) < ep) & (data$Type == levels(data$Type)[i]))]) + 
                                            length(data$Period[which(((1 / data$`X.2.Freq`) < data$Period) &
                                            (abs((data$Period / (1 / data$`X.2.Freq`)) - floor(data$Period / (1 / data$`X.2.Freq`))) < ep) & (data$Type == levels(data$Type)[i]))]) - matchnum2
    
    matchnum3 <- length(data$Period[which((abs(data$Period - (1 / data$`X.3.Freq`)) < ep * data$Period) & (data$Type == levels(data$Type)[i]))])
    
    multnum3 <- length(data$Period[which(((1 / data$`X.3.Freq`) > data$Period) &
                                            (abs(((1 / data$`X.3.Freq`) / data$Period) - floor((1 / data$`X.3.Freq`) / data$Period)) < ep) & (data$Type == levels(data$Type)[i]))]) + 
                                            length(data$Period[which(((1 / data$`X.3.Freq`) < data$Period) &
                                            (abs((data$Period / (1 / data$`X.3.Freq`)) - floor(data$Period / (1 / data$`X.3.Freq`))) < ep) & (data$Type == levels(data$Type)[i]))]) - matchnum3
    
    matchtot <- min(c((matchnum + matchnum2 + matchnum3), totnum))
    
    multtot <- min(c((multnum + multnum2 + multnum3), (totnum - matchtot)))
    
    print(paste("The number of matches for #1: ", matchnum, ", for period #2: ", matchnum2, " and for period #3: ", matchnum3, ".", sep=""))
    
    print(paste("This provides an overall match number of ", matchtot, ".", sep=""))
    
    print(paste("The number of multiples for #1: ", multnum, ", for period #2: ", multnum2, " and for period #3: ", multnum3, ".", sep=""))
    
    print(paste("This provides an overall multiple number of ", multtot, ".", sep=""))
    
    permat2 <- (matchtot / totnum) * 100
    
    permult2 <- (multtot / totnum) * 100
    
    perall2 <- ((matchtot + multtot) / totnum) * 100
    
    print(paste("The matches are ", permat2, "% of the total number.", sep=""))
    
    print(paste("The multiples are ", permult2, "% of the total number.", sep=""))
    
    print(paste("Combined this is ", perall2, "% of the total number.", sep=""))
    
    print(paste("The miss rate is ", 100 - perall2, "% of the total number.", sep=""))
    
    print("--------------------------------------------")
    
  }
  
}





periodmodes <- function(data, ep = 0.01, per = F) {
  
  datalen <- length(data[,1])
  
  print(paste("There are ", datalen, " samples.", sep=""))
  
  data$Type <- factor(data$Type)
  
  name <- as.character(as.vector(data$Name))
  
  type <- as.character(as.vector(data$Type))
  
  if (per == F){
  
    data$Period <- as.numeric(as.vector(data$Period))
  
    data <- data[which(!is.na(data$Period)),]
  
  }
  else{
    
    data <- as.data.frame(cbind(as.numeric(as.vector(data$AAVSO.Period)), 1/as.numeric(as.vector(data$Period))))
    
    colnames(data) <- c("Period", "X.1.Freq")
    
    data <- data[which(!is.na(data$Period)),]
    
  }
  
  postrem <- length(data[,1])
  
  removed <- datalen - postrem
  
  print("Filtering out objects lacking AAVSO periods.")
  
  print(paste("This removes ", removed, " objects leaving ", postrem, " for future operations.", sep=""))
  
  print("Matching periods from AAVSO to those from local Periodogram...")
  
  print(paste("A period will be accepted if within +/-", ep*100, "% of the AAVSO value.", sep=""))
  
  print("Objects without an AAVSO period will be ignored.")
  
  print("--------------------------------------------")
  
  result <- matrix(0, nrow = postrem, ncol = 6)
  
  for (i in 1:postrem){
    
    result[i,6] <- 1
    
    if (abs(data$Period[i] - (1 / data$`X.1.Freq`[i])) < ep * data$Period[i]){
      
      result[i,1] <- 1
      
      result[i,6] <- 0
      
    }
    
    if ((((1 / data$`X.1.Freq`[i]) > data$Period[i]) & (floor((1 / data$`X.1.Freq`[i]) / data$Period[i]) <= 3) & 
         (abs(((1 / data$`X.1.Freq`[i]) / data$Period[i]) - floor((1 / data$`X.1.Freq`[i]) / data$Period[i])) < ep)) |
        (((1 / data$`X.1.Freq`[i]) > data$Period[i]) & (ceiling((1 / data$`X.1.Freq`[i]) / data$Period[i]) <= 4) & 
         (abs(((1 / data$`X.1.Freq`[i]) / data$Period[i]) - ceiling((1 / data$`X.1.Freq`[i]) / data$Period[i])) < ep))){
      
      if (result[i,1] == 0){
        
        result[i,2] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
    if ((((1 / data$`X.1.Freq`[i]) < data$Period[i]) & (floor(data$Period[i] / (1 / data$`X.1.Freq`[i])) <= 3) & 
         (abs((data$Period[i] / (1 / data$`X.1.Freq`[i])) - floor(data$Period[i] / (1 / data$`X.1.Freq`[i]))) < ep)) | 
        (((1 / data$`X.1.Freq`[i]) < data$Period[i]) & (ceiling(data$Period[i] / (1 / data$`X.1.Freq`[i])) <= 4) & 
         (abs((data$Period[i] / (1 / data$`X.1.Freq`[i])) - ceiling(data$Period[i] / (1 / data$`X.1.Freq`[i]))) < ep))){
      
      if (result[i,1] == 0){
        
        result[i,3] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
    if ((abs(abs(data$Period[i]/(1 - data$Period[i])) - (1 / data$`X.1.Freq`[i])) < ep) |
        (abs(abs(data$Period[i]/(1 + data$Period[i])) - (1 / data$`X.1.Freq`[i])) < ep)){
      
      if (result[i,1] == 0){
        
        result[i,4] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
    if ((abs(abs(data$Period[i]/(1 - 2*data$Period[i])) - (1 / data$`X.1.Freq`[i])) < ep) |
        (abs(abs(data$Period[i]/(1 + 2*data$Period[i])) - (1 / data$`X.1.Freq`[i])) < ep)){
      
      if (result[i,1] == 0){
        
        result[i,5] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
  }
  
  print("Operation Complete... Printing out final statistics...")
  
  print(paste("Percentage of objects with a hit match: ", sum(result[,1])*100/postrem, "%", sep=""))
  
  print(paste("Percentage of objects with a multiple match: ", sum(result[,2])*100/postrem, "%", sep=""))
  
  print(paste("Percentage of objects with a submultiple match: ", sum(result[,3])*100/postrem, "%", sep=""))
  
  print(paste("Percentage of objects with a one-day alias match: ", sum(result[,4])*100/postrem, "%", sep=""))
  
  print(paste("Percentage of objects with a half-day alias match: ", sum(result[,5])*100/postrem, "%", sep=""))
  
  print(paste("Percentage of objects with an unknown failure: ", sum(result[,6])*100/postrem, "%", sep=""))
  
  result <- as.data.frame(result)
  
  if (per == F){
  
    result <- cbind(data$Name, data$Type, data$Period, 1/data$`X.1.Freq`, result)
  
  }
  else{
    
    result <- cbind(name, type, data$Period, 1/data$`X.1.Freq`, result)
    
  }
  
  colnames(result) <- c("Name", "Type", "AAVSO.Period", "Per.Period", "Hit", "Multiple", "Sub.Multiple", "Alias.1", "Alias.0.5", "Unknown")
  
  result
  
}




synthpermodes <- function(data, ep = 0.01, grape = T, quiet = F) {
  
  datalen <- length(data[,1])
  
  if (quiet == F){
    
    print(paste("There are ", datalen, " samples.", sep=""))
    
  }
  
  input <- data$Input.Period
  
  if (grape == T){
    
    data$Vuong.Period <- as.numeric(as.vector(data$Vuong.Period))
    
    data <- data[which(!is.na(data$Vuong.Period)),]
    
    pred <- data$Vuong.Period
    
  }
  else{
    
    data$BGLS.Period <- as.numeric(as.vector(data$BGLS.Period))
    
    data <- data[which(!is.na(data$BGLS.Period)),]
    
    pred <- data$BGLS.Period
    
  }
  
  if (quiet == F){
    
    print("Matching periods from Reference to those from local Periodogram...")
    
    print(paste("A period will be accepted if within +/-", ep*100, "% of the AAVSO value.", sep=""))
    
    print("--------------------------------------------")
    
  }
  
  result <- matrix(0, nrow = datalen, ncol = 6)
  
  for (i in 1:datalen){
    
    result[i,6] <- 1
    
    if (abs(input[i] - (pred[i])) < ep * input[i]){
      
      result[i,1] <- 1
      
      result[i,6] <- 0
      
    }
    
    if ((((pred[i]) > input[i]) & (floor((pred[i]) / input[i]) <= 3) & 
         (abs(((pred[i]) / input[i]) - floor((pred[i]) / input[i])) < ep)) |
        (((pred[i]) > input[i]) & (ceiling((pred[i]) / input[i]) <= 4) & 
         (abs(((pred[i]) / input[i]) - ceiling((pred[i]) / input[i])) < ep))){
      
      if (result[i,1] == 0){
        
        result[i,2] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
    if ((((pred[i]) < input[i]) & (floor((input[i]) / pred[i]) <= 3) & 
         (abs((input[i] / (pred[i])) - floor(input[i] / (pred[i]))) < ep)) | 
        (((pred[i]) < input[i]) & (ceiling((input[i]) / pred[i]) <= 4) & 
         (abs((input[i] / (pred[i])) - ceiling(input[i] / (pred[i]))) < ep))){
      
      if (result[i,1] == 0){
        
        result[i,3] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
    if ((abs(abs(input[i]/(1 - input[i])) - (pred[i])) < ep) |
        (abs(abs(input[i]/(1 + input[i])) - (pred[i])) < ep)){
      
      if (result[i,1] == 0){
        
        result[i,4] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
    if ((abs(abs(input[i]/(1 - 2*input[i])) - (pred[i])) < ep) |
        (abs(abs(input[i]/(1 + 2*input[i])) - (pred[i])) < ep)){
      
      if (result[i,1] == 0){
        
        result[i,5] <- 1
        
        result[i,6] <- 0
        
      }
      
    }
    
  }
  
  if (quiet == F){
    
    print("Operation Complete... Printing out final statistics...")
    
    print(paste("Percentage of objects with a hit match: ", sum(result[,1])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with a multiple match: ", sum(result[,2])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with a submultiple match: ", sum(result[,3])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with a one-day alias match: ", sum(result[,4])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with a half-day alias match: ", sum(result[,5])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with an unknown failure: ", sum(result[,6])*100/datalen, "%", sep=""))
    
  }
  
  result <- as.data.frame(result)
  
  result <- cbind(input, pred, result)
  
  colnames(result) <- c("Input.Period", "Per.Period", "Hit", "Multiple", "Sub.Multiple", "Alias.1", "Alias.0.5", "Unknown")
  
  result
  
}







bandspermodes <- function(data, ep = 0.01, quiet = F) {
  
  datalen <- length(data[,1])
  
  if (quiet == F){
    
    print(paste("There are ", datalen, " samples.", sep=""))
    
  }
  
  input <- data$Input.Period
    
  data$Vuong.Period <- as.numeric(as.vector(data$Vuong.Period))
    
  data <- data[which(!is.na(data$Vuong.Period)),]
    
  pred <- data$Vuong.Period
  
  pred2 <- data$Vuong.Period2
  
  if (quiet == F){
    
    print("Matching periods from Reference to those from local Periodogram...")
    
    print(paste("A period will be accepted if within +/-", ep*100, "% of the AAVSO value.", sep=""))
    
    print("--------------------------------------------")
    
  }
  
  result <- matrix(0, nrow = datalen, ncol = 3)
  
  for (i in 1:datalen){
    
    result[i,3] <- 1
    
    if (abs(input[i] - (pred[i])) < ep * input[i]){
      
      result[i,1] <- 1
      
      result[i,2] <- 0
      
      result[i,3] <- 0
      
    }
    
    if (abs(input[i] - (pred2[i])) < ep * input[i]){
      
      result[i,1] <- 0
      
      result[i,2] <- 1
      
      result[i,3] <- 0
      
    }
    
  }
  
  if (quiet == F){
    
    print("Operation Complete... Printing out final statistics...")
    
    print(paste("Percentage of objects with a hit match: ", sum(result[,1])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with a half-period match: ", sum(result[,2])*100/datalen, "%", sep=""))
    
    print(paste("Percentage of objects with an unknown failure: ", sum(result[,3])*100/datalen, "%", sep=""))
    
  }
  
  result <- as.data.frame(result)
  
  result <- cbind(input, pred, result)
  
  colnames(result) <- c("Input.Period", "Per.Period", "Hit", "Half.Hit", "Unknown")
  
  result
  
}






runtime_calc <- function(data) {
  
  runtime <- c(mean(data$Runtime[which(data$Observations >= 100 & data$Observations < 200)]), mean(data$Runtime[which(data$Observations >= 200 & data$Observations < 500)]), mean(data$Runtime[which(data$Observations >= 500 & data$Observations < 1000)]), mean(data$Runtime[which(data$Observations >= 1000 & data$Observations < 2000)]), mean(data$Runtime[which(data$Observations >= 2000 & data$Observations < 3000)]), mean(data$Runtime[which(data$Observations >= 3000 & data$Observations < 4000)]), mean(data$Runtime[which(data$Observations >= 4000 & data$Observations < 5000)]), mean(data$Runtime[which(data$Observations >= 5000 & data$Observations < 8000)]), mean(data$Runtime[which(data$Observations >= 8000 & data$Observations < 11000)]), mean(data$Runtime[which(data$Observations >= 11000 & data$Observations < 15000)]))
  
  runtime
  
}