sellimobj <- function(lim){
  
  library(RODBC)
  
  print(paste("Extracting all Skycam-T objects with at least ", lim, " entries...", sep=""))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  objects <- sqlQuery(channel, paste("Select * from objdat where entries >= '", lim, "'", sep = ""))
  
  n <- nrow(objects)
  
  objdatset <- as.data.frame(cbind(as.character(as.vector(objects$usnoref)), rep("--", n), rep("NEW", n), as.numeric(as.vector(objects$RA)), as.numeric(as.vector(objects$DEClin)), rep("", n), rep("--", n)))
  
  colnames(objdatset) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude")
  
  objdatset <- unique(objdatset)
  
  n <- nrow(objdatset)
  
  objdatset <- objdatset[order(objdatset$Name),]
  
  rownames(objdatset) <- 1:n
  
  print(paste("Operation complete. ", n, " Unique objects selected.", sep=""))
  
  objdatset
  
}