varcomp <- function(x, lim = 100, samp = 10000) {
  
  library(RODBC)
  library(splitstackshape)
  
  set.seed(20)
  
  radii <- 0.01
  
  x$RA <- as.numeric(as.vector(paste(x$RA)))
  
  x$DEC <- as.numeric(as.vector(paste(x$DEC)))
  
  print("Using AAVSO matched dataset for the generation of a non-variable star dataset...")
  
  print(paste("Selecting all objects from Skycam database with greater than ", lim, " observations..."))
  
  radra <- radii / abs(cos(x$DEC*pi/180))
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  nonvars <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, RA, DEClin FROM objdat WHERE entries >= '", lim, "'", sep=""))
  
  totnum <- length(nonvars[,1])
  
  print("Operation Complete.")
  
  print(paste("A total of ", totnum, " objects identified.", sep = ""))
  
  nosamp <- FALSE
  
  if (samp > totnum){
    
    print("Less objects selected than the requested sample size...")
    
    print("Program will not subset the matching results.")
    
    samp <- totnum
    
    nosamp <- TRUE
    
  }
  
  nonvars$RA <- as.numeric(as.vector(paste(nonvars$RA)))
  
  nonvars$DEClin <- as.numeric(as.vector(paste(nonvars$DEClin)))
  
  if (nosamp == FALSE){
  
    print(paste("Sampling ", samp, " of these Skycam objects", sep = ""))
  
    samplen <- length(nonvars$RA)
  
    sampind <- sample(1:samplen, samp)
  
    nonvars <- nonvars[sampind,]
  
  }
  
  print("Removing all known variable objects from the observation-limited Skycam subset...")
  
  for (i in 1:length(x$RA)){
    
    if (i %% 1000 == 0){
      
      print(paste("Matching object number #", i, ".", sep = ""))
      
    }
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                        x$RA[i]-radra[i], "' AND '", x$RA[i]+radra[i], "' AND DEClin BETWEEN '", x$DEC[i]-radii, "' AND '",
                                        x$DEC[i]+radii, "' AND entries >= '", lim, "'", sep=""))
    
    if (length(objects[,1]) > 0){
    
      objects$RA <- as.numeric(as.vector(paste(objects$RA)))
    
      objects$DEC <- as.numeric(as.vector(paste(objects$DEC)))
    
      geodist <- sqrt((x$DEC[i] - objects$DEClin)^2.0 + ((x$RA[i] - objects$RA) / abs(cos(x$DEC[i]*pi/180)))^2.0)
    
      mindist <- which.min(geodist)
    
      objects <- objects[mindist,]
    
      nonvars <- nonvars[which(as.character(nonvars$usnoref) != as.character(objects$usnoref[1])),]
    
    }
  
  }
  
  print("Operation Complete.")
  
  newnum <- length(nonvars[,1])
  
  print(paste("A total of ", samp-newnum, " variable objects have been removed from this dataset.", sep = ""))
  
  print(paste("This leaves ", newnum, " assumed non-variable object(s).", sep = ""))
  
  result <- as.data.frame(cbind(as.character(nonvars$usnoref), rep("", newnum), rep("NV", newnum), as.numeric(nonvars$RA), as.numeric(nonvars$DEClin), runif(newnum, 0.05, 50), rep("", newnum), rep("NV", newnum)))
  
  colnames(result) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude", "Class")
  
  result
  
}