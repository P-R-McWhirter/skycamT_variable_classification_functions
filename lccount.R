lccount <- function(data, lim = 100, radii = 0.1) {
  
  start <- Sys.time()
  
  library(RODBC)
  
  print("Analysing all objects...")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  aper <- as.numeric(as.vector(data$Period))
  
  for (k in 1:length(data$Name)){
    
    print(paste("Processing object ", k, " out of ", length(data$Name), ", named '", data$Name[k], "' of class '", data$Type[k], "'.", sep=""))
    
    ra <- RA[k]
    
    declin <- DEC[k]
    
    radra <- radii / abs(cos(DEC[k]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                       declin+radii, "' AND entries > '", lim, "'", sep=""))
    
    if (length(objects$entries) < 1){
      
      close(channel)
      
      stop(paste("No objects detected within search area of ", radii, " degree(s) with greater than ", lim, " entries. Function will now exit.", sep=""))
      
    }
    
    else{
      
      minmag <- as.numeric(sqlQuery(channel, "SELECT max(r2mag) FROM catdat"))
      
      maxmag <- as.numeric(sqlQuery(channel, "SELECT min(r2mag) FROM catdat"))
      
      brcol <- rep(0, length(objects$usnoref))
      
      rmag <- rep(maxmag, length(objects$usnoref))
      
      namecol <- rep("#00FF00", length(objects$usnoref))
      
      for (i in 1:length(objects$usnoref)){
        
        colandmag <- sqlQuery(channel, paste("SELECT BRcolour, r2mag FROM catdat WHERE usnoref = '", objects$usnoref[i], "'", sep=""))
        
        brcol[i] <- colandmag$BRcolour[1]
        
        rmag[i] <- colandmag$r2mag[1]
        
        if (is.na(rmag[i])){
          
          rmag[i] <- 12
          
        }
        
        if (is.na(brcol[i])){
          
          brcol[i] <- 1000
          
        }
        
        if (brcol[i] <= -0.45){
          
          namecol[i] <- "#0000FF"
          
        }
        else if (brcol[i] > -0.45 & brcol[i] <= -0.2){
          
          namecol[i] <- "#00FFFF"
          
        }
        else if (brcol[i] > -0.2 & brcol[i] <= 0.2){
          
          namecol[i] <- "#FFFFFF"
          
        }
        else if (brcol[i] > 0.2 & brcol[i] <= 0.64){
          
          namecol[i] <- "#FFFF80"
          
        }
        else if (brcol[i] > 0.64 & brcol[i] <= 1.06){
          
          namecol[i] <- "#FFFF00"
          
        }
        else if (brcol[i] > 1.06 & brcol[i] <= 1.42){
          
          namecol[i] <- "#FF8000"
          
        }
        else if (brcol[i] > 1.42 & brcol[i] <= 100){
          
          namecol[i] <- "#FF0000"
          
        }
        else{
          
          namecol[i] <- "#00FF00"
          
        }
        
      }
      
      if (quiet == FALSE){
        
        print(paste("A total of ", length(objects$usnoref),
                    " object(s) identified within the search radius of ", radii, " degree(s) with greater than ", lim, " observations.", sep=""))
        
        par(bg = "black", col = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white")
        
        plot(objects$RA, objects$DEClin, main = paste("Coordinates RA = ", ra, " deg and DEC = ", declin, " deg, radius = ", radii, " deg.", sep=""), xlab = "Right Ascension (degrees)",
             ylab = "Declination (degrees)", pch = 19, cex = ((2.5 * ((rmag - minmag) / (maxmag - minmag)) + 0.5) * max((log10(1 / (2 * radii)) + 1), 0.1)), col = namecol, xlim = c(ra-radra, ra+radra), ylim = c(declin-radii, declin+radii))
        
        if (length(objects$usnoref) < 100){
          
          text(objects$RA, objects$DEClin, labels = objects$usnoref, cex = 0.6, pos = 4, offset = (0.5 * max((log10(1 / (2 * radii)) + 1), 0.01)), col = "#808080")
          
        }
        
        grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
             lwd = par("lwd"), equilogs = TRUE)
        
        par(bg = "white", col = "black", col.axis = "black", col.lab = "black", col.main = "black", col.sub = "black")
        
      }
      
    }
    
  }
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
}





lccounts <- function(data, lim = 100, radii = 0.1) {
  
  start <- Sys.time()
  
  library(RODBC)
  
  print("Analysing all objects...")
  
  channel <- odbcConnect("sktobs", uid = "ross", pwd = "Ccqabw00")
  
  data$Period <- as.numeric(as.vector(data$Period))
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  aper <- as.numeric(as.vector(data$Period))
  
  fullset <- matrix(0, nrow = 0, ncol = 11)
  
  for (k in 1:length(data$Name)){
    
    print(paste("Processing object ", k, " out of ", length(data$Name), ", named '", data$Name[k], "' of class '", data$Type[k], "'.", sep=""))
    
    ra <- RA[k]
    
    declin <- DEC[k]
    
    radra <- radii / abs(cos(DEC[k]*pi/180))
    
    objects <- sqlQuery(channel, paste("SELECT DISTINCT usnoref, entries, RA, DEClin FROM objdat WHERE RA BETWEEN '", 
                                       ra-radra, "' AND '", ra+radra, "' AND DEClin BETWEEN '", declin-radii, "' AND '",
                                       declin+radii, "' AND entries > '", lim, "'", sep=""))
    
    if (length(objects$entries) < 1){
      
      print(paste("No objects detected within search area of ", radii, " degree(s) with greater than ", lim, " entries.", sep=""))
      
    }
    
    else{
      
      geodist <- sqrt((DEC[k] - objects$DEClin)^2.0 + ((RA[k] - objects$RA) / abs(cos(DEC[k]*pi/180)))^2.0)
      
      minmag <- as.numeric(sqlQuery(channel, "SELECT max(r2mag) FROM catdat"))
      
      maxmag <- as.numeric(sqlQuery(channel, "SELECT min(r2mag) FROM catdat"))
      
      brcol <- rep(0, length(objects$usnoref))
      
      rmag <- rep(maxmag, length(objects$usnoref))
      
      namecol <- rep("#00FF00", length(objects$usnoref))
      
      fullfeats <- matrix(0, nrow = length(objects$usnoref), ncol = 11)
      
      for (i in 1:length(objects$usnoref)){
        
        colandmag <- sqlQuery(channel, paste("SELECT BRcolour, r2mag FROM catdat WHERE usnoref = '", objects$usnoref[i], "'", sep=""))
        
        brcol[i] <- colandmag$BRcolour[1]
        
        rmag[i] <- colandmag$r2mag[1]
        
        if (is.na(rmag[i])){
          
          rmag[i] <- 12
          
        }
        
        if (is.na(brcol[i])){
          
          brcol[i] <- 1000
          
        }
        
        if (brcol[i] <= -0.45){
          
          namecol[i] <- "#0000FF"
          
        }
        else if (brcol[i] > -0.45 & brcol[i] <= -0.2){
          
          namecol[i] <- "#00FFFF"
          
        }
        else if (brcol[i] > -0.2 & brcol[i] <= 0.2){
          
          namecol[i] <- "#FFFFFF"
          
        }
        else if (brcol[i] > 0.2 & brcol[i] <= 0.64){
          
          namecol[i] <- "#FFFF80"
          
        }
        else if (brcol[i] > 0.64 & brcol[i] <= 1.06){
          
          namecol[i] <- "#FFFF00"
          
        }
        else if (brcol[i] > 1.06 & brcol[i] <= 1.42){
          
          namecol[i] <- "#FF8000"
          
        }
        else if (brcol[i] > 1.42 & brcol[i] <= 100){
          
          namecol[i] <- "#FF0000"
          
        }
        else{
          
          namecol[i] <- "#00FF00"
          
        }
        
        info <- sqlQuery(channel, paste("SELECT MJD, Rcat, Rcaterr FROM obsdat100split WHERE usnoref = '", 
                                        objects$usnoref[i], "'", sep=""))
        
        ts <- matrix(c(info$MJD, info$Rcat, info$Rcaterr), nrow = length(info$MJD))
        
        prelen <- length(ts[,2])
        
        ts <- unique(ts)
        
        meanmag <- weighted.mean(ts[,2], (1/ts[,3]))
        
        sdmag <- sd(ts[,2])
        
        bright <- median(head(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))
        
        faint <- median(tail(sort(ts[,2]), ceiling(length(ts[,2]) * 0.05)))
        
        fullfeats[i,] <- c(as.character(as.vector(data$Name[k])), as.character(as.vector(data$Type[k])), as.character(as.vector(objects$usnoref[i])), as.numeric(as.vector(objects$RA[i])), as.numeric(as.vector(objects$DEClin[i])), as.numeric(as.vector(objects$entries[i])), sdmag, meanmag, bright, faint, geodist[i]*60*60)
        
      }
      
      fullset <- rbind(fullset, fullfeats)
      
      if (quiet == FALSE){
        
        print(paste("A total of ", length(objects$usnoref),
                    " object(s) identified within the search radius of ", radii, " degree(s) with greater than ", lim, " observations.", sep=""))
        
        par(bg = "black", col = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white")
        
        plot(objects$RA, objects$DEClin, main = paste("Coordinates RA = ", ra, " deg and DEC = ", declin, " deg, radius = ", radii, " deg.", sep=""), xlab = "Right Ascension (degrees)",
             ylab = "Declination (degrees)", pch = 19, cex = ((2.5 * ((rmag - minmag) / (maxmag - minmag)) + 0.5) * max((log10(1 / (2 * radii)) + 1), 0.1)), col = namecol, xlim = c(ra-radra, ra+radra), ylim = c(declin-radii, declin+radii))
        
        if (length(objects$usnoref) < 100){
          
          text(objects$RA, objects$DEClin, labels = objects$usnoref, cex = 0.6, pos = 4, offset = (0.5 * max((log10(1 / (2 * radii)) + 1), 0.01)), col = "#808080")
          
        }
        
        grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
             lwd = par("lwd"), equilogs = TRUE)
        
        par(bg = "white", col = "black", col.axis = "black", col.lab = "black", col.main = "black", col.sub = "black")
        
      }
      
    }
    
  }
  
  colnames(fullset) <- c("Cat.Name", "Cat.Type", "USnoref", "RA", "DEC", "Observations", "SD.Mag", "Mean.Mag", "Bright.Mag", "Faint.Mag", "Dist.from.Cat")
  
  print("Operation Complete.")
  
  close(channel)
  
  print(Sys.time() - start)
  
  fullset <- as.data.frame(fullset)
  
  fullset
  
}