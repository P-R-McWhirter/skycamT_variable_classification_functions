best <- function(state, outcome) {
  
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
  
  range(data[, pos], na.rm = TRUE)
  
  bes <- which.min(data[, pos])
  
  data[bes, 2]
}