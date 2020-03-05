viewOGLElc <- function(rows, cat = "CCEP", periodogram = FALSE, doplot = FALSE){
  
  start <- Sys.time()
  
  bins <- 100
  
  output <- as.data.frame(matrix(NA, nrow = length(rows), ncol = 100))
  
  if (cat == "CCEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_classical_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "T2CEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_type2_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 13
    
  }
  else if (cat == "M"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_miras.txt", header=FALSE, comment.char="#")
    
    set <- 16
    
  }
  else if (cat == "SR"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_semiregulars.txt", header=FALSE, comment.char="#")
    
    set <- 16
  
  }
  else if (cat == "RR"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_rrlyraes.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else{
    
    stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
    
  }
  
  input <- input[rows,]
  
  input[,1] <- as.character(as.vector(input[,1]))
  
  for (i in 1:nrow(input)){
  
    inputrow <- input[i,]
  
    object <- as.character(inputrow[1])
  
    period <- as.numeric(inputrow[set])
  
    print(paste("Object : '", i, "'. Catalogue period = ", period, " days.", sep = ""))
  
    if (cat == "CCEP"){
  
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_classical_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
  
    }
    else if (cat == "T2CEP"){
    
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_type2_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
    
    }
    else if (cat == "M"){
    
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_miras/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
    
    }
    else if (cat == "SR"){
    
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_semiregulars/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
    
    }
    else if (cat == "RR"){
    
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_rrlyraes/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
    
    }
    else{
    
      stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
    
    }
  
    if (periodogram == TRUE){
    
      library(lomb)
    
      baseline = max(lightcurve[,1]) - min(lightcurve[,1])
    
      freq <- lsp(lightcurve[,1:2], from = 1/baseline, to = 20, ofac = 5, plot = FALSE)$peak.at[1]
      
      perlsp <- lsp(lightcurve[,1:2], from = 0.99*freq, to = 1.01*freq, ofac = 20, plot = FALSE)$peak.at[2]
    
      print(paste("Periodogram period = ", perlsp, " days.", sep = ""))
    
      if (abs(perlsp - period) > 5e-3){
      
        period <- perlsp
      
      }
    
    }
    
    fts <- cbind((((lightcurve[,1] - lightcurve[which.min(lightcurve[,2]),1])/period) + 0.25)%%1 - 0.5, lightcurve[,2], lightcurve[,3])
  
    meanmag <- weighted.mean(fts[,2], 1/(fts[,3]^2.0))
    
    if (doplot == TRUE){
    
      plot(fts[,1:2], pch=19, xlim = c(-0.5, 0.5), ylim = c(max(fts[,2]), min(fts[,2])), xlab = paste("Phase at a period of ", period, " days", sep = ""), ylab = "Apparent Magnitude (mag)")
    
    }
    
    fts[,2] <- fts[,2] - meanmag
    
    co <- harm8fit_lin(fts, 1)
    
    tc <- seq(-0.5, 0.5, length.out = 100)
    
    fit <- co[1] + 
      co[2]*sin(2*pi*tc) + co[3]*cos(2*pi*tc) + 
      co[4]*sin(4*pi*tc) + co[5]*cos(4*pi*tc) + 
      co[6]*sin(6*pi*tc) + co[7]*cos(6*pi*tc) + 
      co[8]*sin(8*pi*tc) + co[9]*cos(8*pi*tc) + 
      co[10]*sin(10*pi*tc) + co[11]*cos(10*pi*tc) + 
      co[12]*sin(12*pi*tc) + co[13]*cos(12*pi*tc) + 
      co[14]*sin(14*pi*tc) + co[15]*cos(14*pi*tc) + 
      co[16]*sin(16*pi*tc) + co[17]*cos(16*pi*tc)
    
    if (doplot == TRUE){
      
      plot(tc, fit, type = "l", lwd = 2)
    
    }
    
    output[i,] <- fit
    
  }
  
  print(Sys.time() - start)
  
  output
  
}





viewOGLElc_new <- function(rows, cat = "CCEP", periodogram = FALSE, doplot = FALSE){
  
  start <- Sys.time()
  
  bins <- 100
  
  output <- as.data.frame(matrix(NA, nrow = length(rows), ncol = 100))
  
  if (cat == "CCEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_classical_cepheids.txt", header=FALSE, comment.char="#")
    
    input <- input[which(input$V4 == "F"),]
    
    set <- 14
    
  }
  else if (cat == "T2CEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_type2_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 13
    
  }
  else if (cat == "ACEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_anomalous_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "DSCT"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_delta_scutis.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "RR"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/blg_rrlyraes.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else{
    
    stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
    
  }
  
  input <- input[rows,]
  
  input[,1] <- as.character(as.vector(input[,1]))
  
  for (i in 1:nrow(input)){
    
    inputrow <- input[i,]
    
    object <- as.character(inputrow[1])
    
    period <- as.numeric(inputrow[set])
    
    print(paste("Object : '", i, "'. Catalogue period = ", period, " days.", sep = ""))
    
    if (cat == "CCEP"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_classical_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "T2CEP"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_type2_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "ACEP"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_anomalous_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "DSCT"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_delta_scutis/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "RR"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/blg_rrlyraes/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else{
      
      stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
      
    }
    
    if (class(lightcurve) == "try-error"){
      
      next
      
    }
    
    if (periodogram == TRUE){
      
      library(lomb)
      
      baseline = max(lightcurve[,1]) - min(lightcurve[,1])
      
      freq <- lsp(lightcurve[,1:2], from = 1/baseline, to = 20, ofac = 5, plot = FALSE)$peak.at[1]
      
      perlsp <- lsp(lightcurve[,1:2], from = 0.99*freq, to = 1.01*freq, ofac = 20, plot = FALSE)$peak.at[2]
      
      print(paste("Periodogram period = ", perlsp, " days.", sep = ""))
      
      if (abs(perlsp - period) > 5e-3){
        
        period <- perlsp
        
      }
      
    }
    
    fts <- cbind((((lightcurve[,1] - lightcurve[which.min(lightcurve[,2]),1])/period) + 0.25)%%1 - 0.5, lightcurve[,2], lightcurve[,3])
    
    meanmag <- weighted.mean(fts[,2], 1/(fts[,3]^2.0))
    
    if (doplot == TRUE){
      
      plot(fts[,1:2], pch=19, xlim = c(-0.5, 0.5), ylim = c(max(fts[,2]), min(fts[,2])), xlab = paste("Phase at a period of ", period, " days", sep = ""), ylab = "Apparent Magnitude (mag)")
      
    }
    
    fts[,2] <- fts[,2] - meanmag
    
    #co <- harm16fit_lin(fts, 1, lambda = 1e-4)
    
    #tc <- seq(-0.5, 0.5, length.out = 100)
    
    #fit <- co[1] + 
      #co[2]*sin(2*pi*tc) + co[3]*cos(2*pi*tc) + 
      #co[4]*sin(4*pi*tc) + co[5]*cos(4*pi*tc) + 
      #co[6]*sin(6*pi*tc) + co[7]*cos(6*pi*tc) + 
      #co[8]*sin(8*pi*tc) + co[9]*cos(8*pi*tc) + 
      #co[10]*sin(10*pi*tc) + co[11]*cos(10*pi*tc) + 
      #co[12]*sin(12*pi*tc) + co[13]*cos(12*pi*tc) + 
      #co[14]*sin(14*pi*tc) + co[15]*cos(14*pi*tc) + 
      #co[16]*sin(16*pi*tc) + co[17]*cos(16*pi*tc)
    
    #if (doplot == TRUE){
      
      #plot(tc, fit, type = "l", lwd = 2)
      
    #}
    
    #output[i,] <- fit
    
    x <- c(fts[,1] - 1, fts[,1], fts[,1] + 1)
    
    y <- c(fts[,2], fts[,2], fts[,2])
    
    output[i,] <- as.numeric(as.vector(approx(x, y, xout = seq(-0.5, 0.5, length.out = 100))$y))
    
  }
  
  print(Sys.time() - start)
  
  output
  
}







viewOGLElc_new2 <- function(rows, cat = "CCEP", periodogram = FALSE, doplot = FALSE){
  
  start <- Sys.time()
  
  bins <- 101
  
  output <- as.data.frame(matrix(NA, nrow = length(rows), ncol = 100))
  
  if (cat == "CCEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_classical_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "T2CEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_type2_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 13
    
  }
  else if (cat == "ACEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_anomalous_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "DSCT"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/all_delta_scutis.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "RR"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/newset/blg_rrlyraes.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else{
    
    stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
    
  }
  
  input <- input[rows,]
  
  input[,1] <- as.character(as.vector(input[,1]))
  
  for (i in 1:nrow(input)){
    
    inputrow <- input[i,]
    
    object <- as.character(inputrow[1])
    
    period <- as.numeric(inputrow[set])
    
    print(paste("Object : '", i, "'. Catalogue period = ", period, " days.", sep = ""))
    
    if (cat == "CCEP"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_classical_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "T2CEP"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_type2_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "ACEP"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_anomalous_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "DSCT"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/all_delta_scutis/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else if (cat == "RR"){
      
      lightcurve <- try(read.table(paste("F:/Documents/BLAP/OGLE_variables/newset/blg_rrlyraes/I/", object, ".dat", sep = ""), quote="\"", comment.char=""), TRUE)
      
    }
    else{
      
      stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
      
    }
    
    if (class(lightcurve) == "try-error"){
      
      next
      
    }
    
    if (periodogram == TRUE){
      
      library(lomb)
      
      baseline = max(lightcurve[,1]) - min(lightcurve[,1])
      
      freq <- lsp(lightcurve[,1:2], from = 1/baseline, to = 20, ofac = 5, plot = FALSE)$peak.at[1]
      
      perlsp <- lsp(lightcurve[,1:2], from = 0.99*freq, to = 1.01*freq, ofac = 20, plot = FALSE)$peak.at[2]
      
      print(paste("Periodogram period = ", perlsp, " days.", sep = ""))
      
      if (abs(perlsp - period) > 5e-3){
        
        period <- perlsp
        
      }
      
    }
    
    fts <- cbind((((lightcurve[,1] - lightcurve[which.min(lightcurve[,2]),1])/period) + 0.25)%%1 - 0.5, lightcurve[,2], lightcurve[,3])
    
    meanmag <- weighted.mean(fts[,2], 1/(fts[,3]^2.0))
    
    if (doplot == TRUE){
      
      plot(fts[,1:2], pch=19, xlim = c(-0.5, 0.5), ylim = c(max(fts[,2]), min(fts[,2])), xlab = paste("Phase at a period of ", period, " days", sep = ""), ylab = "Apparent Magnitude (mag)")
      
    }
    
    fts[,2] <- fts[,2] - meanmag
    
    splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
    
    tery <- rep(0, (bins-1))
    
    tere <- tery
    
    tert <- tery
    
    for (k in 1:(bins-1)){
      
      tery[k] <- weighted.mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2, drop = F], 1/(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F])^2.0)
      
      tere[k] <- weighted.mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F], 1/(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F])^2.0)
      
      tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
      
    }
    
    fts <- cbind(tert, tery, tere)
    
    fts <- fts[order(fts[,1]), , drop = F]
    
    output[i,] <- fts[,2]
    
  }
  
  print(Sys.time() - start)
  
  output
  
}



viewOGLElc_old <- function(rows, cat = "CCEP", periodogram = FALSE){
  
  start <- Sys.time()
  
  bins <- 100
  
  output <- as.data.frame(matrix(NA, nrow = length(rows), ncol = 100))
  
  if (cat == "CCEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_classical_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else if (cat == "T2CEP"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_type2_cepheids.txt", header=FALSE, comment.char="#")
    
    set <- 13
    
  }
  else if (cat == "M"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_miras.txt", header=FALSE, comment.char="#")
    
    set <- 16
    
  }
  else if (cat == "SR"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_semiregulars.txt", header=FALSE, comment.char="#")
    
    set <- 16
    
  }
  else if (cat == "RR"){
    
    input <- read.delim("F:/Documents/BLAP/OGLE_variables/blg_rrlyraes.txt", header=FALSE, comment.char="#")
    
    set <- 14
    
  }
  else{
    
    stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
    
  }
  
  input <- input[rows,]
  
  input[,1] <- as.character(as.vector(input[,1]))
  
  for (i in 1:nrow(input)){
    
    inputrow <- input[i,]
    
    object <- as.character(inputrow[1])
    
    period <- as.numeric(inputrow[set])
    
    print(paste("Catalogue period = ", period, " days.", sep = ""))
    
    if (cat == "CCEP"){
      
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_classical_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
      
    }
    else if (cat == "T2CEP"){
      
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_type2_cepheids/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
      
    }
    else if (cat == "M"){
      
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_miras/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
      
    }
    else if (cat == "SR"){
      
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_semiregulars/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
      
    }
    else if (cat == "RR"){
      
      lightcurve <- read.table(paste("F:/Documents/BLAP/OGLE_variables/blg_rrlyraes/I/", object, ".dat", sep = ""), quote="\"", comment.char="")
      
    }
    else{
      
      stop(paste("'", cat, "' is not a valid category selection.", sep = ""))
      
    }
    
    if (periodogram == TRUE){
      
      library(lomb)
      
      baseline = max(lightcurve[,1]) - min(lightcurve[,1])
      
      freq <- lsp(lightcurve[,1:2], from = 1/baseline, to = 20, ofac = 5, plot = FALSE)$peak.at[1]
      
      perlsp <- lsp(lightcurve[,1:2], from = 0.99*freq, to = 1.01*freq, ofac = 20, plot = FALSE)$peak.at[2]
      
      print(paste("Periodogram period = ", perlsp, " days.", sep = ""))
      
      if (abs(perlsp - period) > 5e-3){
        
        period <- perlsp
        
      }
      
    }
    
    fts <- cbind((((lightcurve[,1] - lightcurve[which.min(lightcurve[,2]),1])/period) + 0.25)%%1 - 0.5, lightcurve[,2], lightcurve[,3])
    
    plot(fts[,1:2], pch=19, xlim = c(-0.5, 0.5), ylim = c(max(fts[,2]), min(fts[,2])), xlab = paste("Phase at a period of ", period, " days", sep = ""), ylab = "Apparent Magnitude (mag)")
    
    splitphase <- seq(from = -0.5, to = 0.5, length.out = bins)
    
    tery <- rep(0, (bins-1))
    
    tere <- tery
    
    tert <- tery
    
    for (k in 1:(bins-1)){
      
      tery[k] <- weighted.mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),2, drop = F], 1/(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F])^2.0)
      
      tere[k] <- weighted.mean(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F], 1/(fts[which((fts[,1] > splitphase[k]) & (fts[,1] <= splitphase[k+1])),3, drop = F])^2.0)
      
      tert[k] <- (splitphase[k] + splitphase[k+1]) / 2.0
      
    }
    
    fts <- cbind(tert, tery, tere)
    
    row.has.na <- apply(fts, 1, function(x){any(is.na(x))})
    
    fts <- fts[!row.has.na, , drop = F]
    
    fts <- fts[order(fts[,1]), , drop = F]
    
    meanmag <- weighted.mean(fts[,2], 1/(fts[,3]^2.0))
    
    x <- c(fts[,1] - 1, fts[,1], fts[,1] + 1)
    
    y <- c(fts[,2], fts[,2], fts[,2])
    
    err <- c(fts[,2], fts[,2], fts[,2])
    
    kern <- Matern32$new(0.15)
    gpk <- GauPro_kernel_model$new(matrix(x, ncol=1), y, kernel=kern, parallel=TRUE)
    if (requireNamespace("MASS", quietly = TRUE)) {
      plot(gpk, xmin = -1.0, xmax = 1.0)
    }
    
    output[i,] <- as.numeric(as.vector(gpk$predict(seq(-0.5, 0.5, length.out = 100)) - meanmag))
    
  }
  
  print(Sys.time() - start)
  
  output
  
}