LCPE <- function(data, radii = 0.001, lim = 100, fit = "lsp", jit = 0, noimod = 1, sigk = 2.5, seascoe, trendrem = FALSE, rangcut = 0.2, redo = 5, seed = 20, pop = 20, ps = 20, nogen = 100, crov = 0.65, mut = 0.03, fdif = 0.6, dfrac = 0.7, clipcut = 0.2, bins = 20, usenyq = FALSE, quiet = TRUE, fulltest = TRUE, tenper = FALSE) {
  
  start <- Sys.time()
  
  data$Period <- as.numeric(as.character(data$Period))
  
  data <- data[!is.na(data$Period),]
  
  print("Initializing light curve period estimation using a hybrid genetic algorithm...")
  
  if (length(crov) <= 1){
    
    print("Setting Crossover rate to be identical across all generations.")
    
    crov <- rep(crov, nogen)
    
  }
  else if (length(crov) > nogen){
    
    print("Warning: Crossover rate vector length larger than the max generations.")
    
    print("All values for invalid generations are being removed...")
    
    crov <- crov[1:nogen]
    
  }
  
  if (length(mut) <= 1){
    
    print("Setting Mutation rate to be identical across all generations.")
    
    mut <- rep(mut, nogen)
    
  }
  else if (length(mut) > nogen){
    
    print("Warning: Mutation rate vector length larger than the max generations.")
    
    print("All values for invalid generations are being removed...")
    
    mut <- mut[1:nogen]
    
  }
  
  if (fit == "lsp"){
    
    print("Using Lomb-Scargle model for fitness function...")
    
  }
  if (fit == "bgls"){
    
    print("Using Bayesian Generalised Lomb-Scargle model for fitness function...")
    
  }
  else if (fit == "ce"){
    
    print("Using Conditional Entropy for fitness function...")
    
  }
  else if (fit == "ml"){
    
    print("Using Machine Learned templates for fitness function...")
    
  }
  else{
    
    print("Invalid fitness function selected. Overriding...")
    
    print("Using Bayesian Generalised Lomb-Scargle model for fitness function...")
    
    fit <- "bgls"
    
  }
  
  print("Using settings:")
  
  print(paste("Repeats = ", redo, ". Seed = ", seed, ".", sep = ""))
  
  print(paste("Population = ", pop, ". Pairups = ", ps, ".", sep = ""))
  
  print(paste("No. of Generations = ", nogen, ". Crossover Rate = ", sep = ""))
  
  print(crov)
  
  print("Mutation Rate = ")
  
  print(mut)
  
  print(paste("Selection Pressure = ", fdif, ". Death Fraction = ", dfrac, ".", sep = ""))
  
  if (fulltest == TRUE){
    
    print(paste("Boolean 'fulltest' has been set to TRUE. Algorithm will determine candidate period from model comparison with the ", redo, " candidate periods. This will require the computation of ", redo*redo, " Vuong-Closeness tests.", sep=""))
    
  }
  else{
    
    print(paste("Boolean 'fulltest' has been set to FALSE. Algorithm will determine candidate period from model comparison with a single constant model. This will require the computation of ", redo, " Vuong-Closeness tests.", sep=""))
    
  }
  
  result <- data.frame()
  
  RA <- as.numeric(as.vector(data$RA))
  
  DEC <- as.numeric(as.vector(data$DEC))
  
  for (i in 1:length(data$Name)){
    
    print(paste("Processing object ", i, " out of ", length(data$Name), ", named '", data$Name[i], "' of class '", data$Type[i], "'.", sep=""))
    
    object <- LCPEAN(RA[i], DEC[i], radii, lim, fit = fit, jit = jit, noimod = noimod, sigk = sigk, seascoe = seascoe, trendrem = trendrem, rangcut = rangcut, redo = redo, seedno = seed, pop = pop, ps = ps, nogen = nogen, crov = crov, mut = mut, fdif = fdif, dfrac = dfrac, clipcut = clipcut, bins = bins, alldb = FALSE, geo = TRUE, quiet = quiet, fulltest = TRUE, tenper = tenper)
    
    object <- as.data.frame(object)
    
    object <- cbind(data$Name[i], data$Type[i], data$Period[i], object[1,])
    
    result <- rbind(result, object[1,])
    
  }
  
  posalias <- abs((1 + as.numeric(data$Period)^-1)^-1)
  
  negalias <- abs((1 - as.numeric(data$Period)^-1)^-1)
  
  posalias2 <- abs((2 + as.numeric(data$Period)^-1)^-1)
  
  negalias2 <- abs((2 - as.numeric(data$Period)^-1)^-1)
  
  result <- cbind(result, posalias, negalias, posalias2, negalias2)
  
  colnames(result) <- c("Name", "Type", "Period", "USnoref", "#.Observations", "Output.Periods", "Fit.Stat", "Genetic.Period", "Vuong.Period", "Vuong.Confidence", "Vuong.Day.Confidence", "Num.Spurious", "Nearest.Spurious", "Period.+1Alias", "Period.-1Alias", "Period.+0.5Alias", "Period.-0.5Alias")
  
  print("Operation Complete.")
  
  print(Sys.time() - start)
  
  result
  
}