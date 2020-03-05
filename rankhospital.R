rankhospital <- function(state, outcome, num) {
  
  data <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
  
  valid <- c("heart attack", "heart failure", "pneumonia")
  
  if(!(state %in% data$State && outcome %in% valid)) {
    
    if (state %in% data$State){
      err = "invalid outcome"
    }
    else{
      err = "invalid state"
    }
    
    stop(err, call. = TRUE, domain = NULL)
    geterrmessage()
  }
  
  data <- data[data$State == state,]
  
  if (outcome == "heart attack"){
    pos = 11
  }
  else if (outcome == "heart failure"){
    pos = 17
  }
  else {
    pos = 23
  }
  
  data[, pos] <- as.numeric(data[, pos])
  
  bad <- is.na(data[,pos])
  
  data <- data[!bad,]
  
  range(data[, pos], na.rm = TRUE)
  
  if (num == "worst"){
    num = length(data[, 2])
  }
  else if (num == "best"){
    num = 1
  }
  
  if (num > length(data[, 2])){
    return(NA)
  }
  
  data <- data[order(data[,pos], data[,2]),]
  
  bes <- which.min(data[, pos])
  
  data[num, 2]
  
}