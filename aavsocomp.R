aavsocomp <- function(radii = 0.01, lim = 100) {
  
  library(RODBC)
  library(splitstackshape)
  
  print("Reading in AAVSO dataset for comparison to STILT data...")
  
  AAVSO_vars <- read.csv("F:/Documents/PhD/AAVSO/var.csv")
  
  print("Operation successful...")
  
  print("Converting RA and DEC into degrees...")
  
  AAVSO_vars <- cSplit(AAVSO_vars, "Coords", " ")
  
  dRA <- abs(as.numeric(AAVSO_vars$Coords_1)) + (as.numeric(AAVSO_vars$Coords_2) / 60) + (as.numeric(AAVSO_vars$Coords_3) / 3600)
  
  dDEC <- sign(AAVSO_vars$Coords_4) * (abs(as.numeric(AAVSO_vars$Coords_4)) + (as.numeric(AAVSO_vars$Coords_5) / 60) + (as.numeric(AAVSO_vars$Coords_6) / 3600))
  
  dRA <- (dRA * 360) / 24
  
  AAVSO_vars <- cbind(as.character(AAVSO_vars$Name), as.character(AAVSO_vars$AUID), as.character(AAVSO_vars$Type), as.numeric(dRA), as.numeric(dDEC), as.character(AAVSO_vars$Period), as.character(AAVSO_vars$Mag))
  
  AAVSO_vars <- as.data.frame(AAVSO_vars)
  
  colnames(AAVSO_vars) <- c("Name","AUID", "Type", "RA", "DEC", "Period", "Magnitude")
  
  print("Conversion Complete...")
  
  print("Comparing coordinates of AAVSO catalogue to STILT database...")
  
  print(paste("Selecting only variable stars with more than ", lim, " observations...", sep = ""))
  
  radra <- radii / abs(cos(dDEC*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  present <- rep.int(0, length(AAVSO_vars$Name))
  
  for (i in 1:length(AAVSO_vars$Name)){
    
    if ((i %% 1000) == 0){
      
      print(i)
      
    }
    
    objcount <- sqlQuery(channel, paste("SELECT DISTINCT count(*) FROM objdat WHERE RA BETWEEN '", 
                                        dRA[i]-radra[i], "' AND '", dRA[i]+radra[i], "' AND DEClin BETWEEN '", dDEC[i]-radii, "' AND '",
                                        dDEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    objcount <- as.numeric(objcount)
    
    if (objcount > 0){
      
      present[i] <- 1
      
    }
    
  }
  
  AAVSO_vars[present == 1,]
  
}





rcdscomp <- function(data, radii = 0.01, lim = 100) {
  
  library(RODBC)
  library(splitstackshape)
  
  print("Reading in Richards dataset for comparison to STILT data...")
  
  print("Operation successful...")
  
  print("Comparing coordinates of Richards catalogue to STILT database...")
  
  print(paste("Selecting only variable stars with more than ", lim, " observations...", sep = ""))
  
  dRA <- data$RA
  
  dDEC <- data$DEC
  
  radra <- radii / abs(cos(dDEC*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  present <- rep.int(0, length(dRA))
  
  for (i in 1:length(dRA)){
    
    if ((i %% 1000) == 0){
      
      print(i)
      
    }
    
    objcount <- sqlQuery(channel, paste("SELECT DISTINCT count(*) FROM objdat WHERE RA BETWEEN '", 
                                        dRA[i]-radra[i], "' AND '", dRA[i]+radra[i], "' AND DEClin BETWEEN '", dDEC[i]-radii, "' AND '",
                                        dDEC[i]+radii, "' AND entries > '", lim, "'", sep=""))
    
    objcount <- as.numeric(objcount)
    
    if (objcount > 0){
      
      present[i] <- 1
      
    }
    
  }
  
  data[present == 1,]
  
}