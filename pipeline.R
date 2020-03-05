pipeline <- function(varcat, lim = 100, ep = 0.05, seedno = 10, lctrenddata, varinddata = NULL, fullcatdata = NULL, optfn = "F1", clsoset = 0, clsaoset = 0, varoset = 0, onlymatch = FALSE, usecolour = FALSE){
  
  library(RODBC)
  
  set.seed(seedno)
  
  print("Executing Skycam Automated Variability Classification Pipeline.")
  
  # Firstly we want to generate a full set of variability indices for every object in the database
  # with more than lim observations. Only if varinddata doesnt exist.
  
  if (is.null(varinddata)){
    
    objdatset <- sellimobj(lim)
    
    print("Calculating Variability Indices for these objects... **This may take some time**")
    
    print("This can be overridden using the varinddata function argument for already generated data.")
    
    objdata <- featgen(objdatset, radii = 0.042, quiet = T, sigk = 4, nonper = T, lctrenddata = lctrenddata)
    
    objdata <- objdata[complete.cases(objdata[,-3]),]
    
  }
  else{
    
    print("Variability Indices supplied for a set of Skycam-T objects. Using these...")
    
    n <- nrow(varinddata)
    
    objdatset <- as.data.frame(cbind(as.character(as.vector(varinddata$USnoref)), rep("--", n), rep("NEW", n), as.numeric(as.vector(varinddata$Right.Ascension)), as.numeric(as.vector(varinddata$Declination)), rep("", n), rep("--", n)))
    
    colnames(objdatset) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude")
    
    print("Data Imported.")
    
    objdata <- varinddata[complete.cases(varinddata[,-3]),]
    
  }
  
  if (usecolour == TRUE){
  
    if (any(is.na(objdata$B.R.Colour)) | any(objdata$B.R.Colour == 1000)){
    
      print("Missing B-R Colour information. Imputing missing values using missForest...")
    
      objdata <- colimpute(objdata, maxiter = 10, ntree = 100)
    
      print("Imputation Complete.")
    
    }
    else{
    
      print("No Missing data for B-R Colour. Skipping missForest imputation procedure.")
    
    }
    
  }
  else{
    
    objdata <- objdata[,-which(names(objdata) %in% "B.R.Colour")]
    
  }
  
  print("Using catalogued Variable Stars to construct classification calibration set...")
  
  if (is.null(fullcatdata)){
    
    print("Calculating full feature set for catalogue objects... **This may take some time**")
    
    print("This can be overridden using the fullcatdata function argument for already generated data.")
    
    vardata <- featgen(varcat, radii = 0.042, quiet = T, sigk = 4, nonper = F, lctrenddata = lctrenddata)
    
    print(vardata)
    
    vardata <- vardata[complete.cases(vardata),]
    
  }
  else{
    
    print("Classification Caibration set full features supplied for a set of Skycam-T objects. Using these...")
    
    n <- nrow(fullcatdata)
    
    varcat <- as.data.frame(cbind(as.character(as.vector(fullcatdata$Name)), rep("--", n), as.character(as.vector(fullcatdata$Type)), as.numeric(as.vector(fullcatdata$Right.Ascension)), as.numeric(as.vector(fullcatdata$Declination)), as.numeric(as.vector(fullcatdata$AAVSO.Period)), as.numeric(as.vector(fullcatdata$Median.Magnitude))))
    
    colnames(objdatset) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude")
    
    print("Data Imported.")
    
    vardata <- fullcatdata[complete.cases(fullcatdata),]
    
  }
  
  if (usecolour == TRUE){
  
    print("Calibrating the B-R Colours to the Skycam-T variability index data...")
    
    for (i in 1:nrow(vardata)){
      
      vardata$B.R.Colour[i] <- objdata$B.R.Colour[which(objdata$USnoref == vardata$USnoref[i])]
      
    }
    
    print("Colour calibration complete.")
  
  }
  
  print(paste(nrow(vardata), " classification calibration objects present.", sep = ""))
  
  print(paste("These objects cover ", length(levels(as.factor(vardata$Type))), " known variability classes.", sep = ""))
  
  print("Selecting catalogue matching objects using the period feature (allowing for Hits, Multiples and Aliases)...")
  
  n.old <- nrow(vardata)
  
  permod <- periodmodes(vardata, ep = ep, per = T)
      
  traindata <- vardata[which(permod$Hit == 1 | permod$Multiple == 1 | permod$Sub.Multiple == 1),]
  
  trainadata <- vardata[which(permod$Alias.1 == 1 | permod$Alias.0.5 == 1),]
  
  n <- nrow(traindata)

  n2 <- nrow(trainadata)
  
  print(paste("This filter results in ", n, " matching classification training objects with ", n.old-n, " removed objects.", sep = ""))
  
  print(paste("This filter results in ", n2, " aliased classification training objects with ", n.old-n2, " removed objects.", sep = ""))
  
  print("Beginning the training of the variable star classifiers...")
  
  print("Optimising the classifier hyperparameters using a 5-fold cross validation with 2 repeats.")
  
  print(paste("'", optfn, "' selected as the optimisation statistic.", sep = ""))
  
  bigcv <- bigparamsel(traindata[,c(2, 10:ncol(traindata))], optfn = optfn, oset = clsoset)
  
  bigcvstats <- bigparamselmean(bigcv$result)
  
  bigacv <- bigparamsel(trainadata[,c(2, 10:ncol(trainadata))], optfn = optfn, oset = clsaoset)
  
  bigacvstats <- bigparamselmean(bigacv$result)
  
  if (optfn == "BS"){
    
    beststat <- min(bigcvstats$meanvals)
    
    bestans <- which.min(bigcvstats$meanvals)
    
    bestastat <- min(bigacvstats$meanvals)
    
    bestaans <- which.min(bigacvstats$meanvals)
    
  }
  else{
    
    beststat <- max(bigcvstats$meanvals)
    
    bestans <- which.max(bigcvstats$meanvals)
    
    bestastat <- max(bigacvstats$meanvals)
    
    bestaans <- which.max(bigacvstats$meanvals)
    
  }
  
  print("Fine-tuning the hyperparameters around the optimal values...")
  
  bigcv <- bigparamsel(traindata[,c(2, 10:ncol(traindata))], optfn = optfn, oset = clsoset, mtrs = seq(as.numeric(bigcv$config$mtrs[bestans])-2, as.numeric(bigcv$config$mtrs[bestans])+2, by = 1), ntrees = c(250, 500, 1000, 2000), nodesizes = seq(as.numeric(bigcv$config$nodesizes[bestans])-2, as.numeric(bigcv$config$nodesizes[bestans])+2, by = 1))
  
  bigcvstats <- bigparamselmean(bigcv$result)
  
  print(bigcv$cutoffs)
  
  bigacv <- bigparamsel(trainadata[,c(2, 10:ncol(trainadata))], optfn = optfn, oset = clsaoset, mtrs = seq(as.numeric(bigacv$config$mtrs[bestaans])-2, as.numeric(bigacv$config$mtrs[bestaans])+2, by = 1), ntrees = c(250, 500, 1000, 2000), nodesizes = seq(as.numeric(bigacv$config$nodesizes[bestaans])-2, as.numeric(bigacv$config$nodesizes[bestaans])+2, by = 1))
  
  bigacvstats <- bigparamselmean(bigacv$result)
  
  print(bigacv$cutoffs)
  
  print("Training final model using the optimum hyperparameters...")
  
  if (optfn == "BS"){
    
    beststat <- min(bigcvstats$meanvals)
    
    bestans <- which.min(bigcvstats$meanvals)
    
    bestastat <- min(bigacvstats$meanvals)
    
    bestaans <- which.min(bigacvstats$meanvals)
    
  }
  else{
    
    beststat <- max(bigcvstats$meanvals)
    
    bestans <- which.max(bigcvstats$meanvals)
    
    bestastat <- max(bigacvstats$meanvals)
    
    bestaans <- which.max(bigacvstats$meanvals)
    
  }
  
  print(paste("Best matching hyperparameters located with optimisation statistic of ", beststat, ".", sep = ""))
  
  print(paste("These are: mtr = ", bigcv$config$mtrs[bestans], ", ntrees = ", bigcv$config$ntrees[bestans], ", and nodesize = ", bigcv$config$nodesizes[bestans], ".", sep = ""))
  
  print(paste("Best aliased hyperparameters located with optimisation statistic of ", bestastat, ".", sep = ""))
  
  print(paste("These are: mtr = ", bigacv$config$mtrs[bestaans], ", ntrees = ", bigacv$config$ntrees[bestaans], ", and nodesize = ", bigacv$config$nodesizes[bestaans], ".", sep = ""))
  
  clscutoffs <- as.vector(bigcv$cutoffs[[bestans]])
  
  print(clscutoffs)
  
  clsacutoffs <- as.vector(bigacv$cutoffs[[bestaans]])
  
  print(clsacutoffs)
  
  classmodel <- crossval(traindata[,c(2, 10:ncol(traindata))], train = T, mtr = as.numeric(bigcv$config$mtrs[bestans]), ntree = as.numeric(bigcv$config$ntrees[bestans]), nodesize = as.numeric(bigcv$config$nodesizes[bestans]), quiet = T, featsel = T, oset = clsoset)
  
  classamodel <- crossval(trainadata[,c(2, 10:ncol(trainadata))], train = T, mtr = as.numeric(bigacv$config$mtrs[bestaans]), ntree = as.numeric(bigacv$config$ntrees[bestaans]), nodesize = as.numeric(bigacv$config$nodesizes[bestaans]), quiet = T, featsel = T, oset = clsaoset)
  
  print("Training of variable star classification models complete.")
  
  print("Calibrating the probabilities of the matching classification model using the Brier-Score of a 5-fold cross validation on the training set.")
  
  cvpreds <- crossval(traindata[,c(2, 10:ncol(traindata))], seed = seedno, k=5, cvpred = T, mtr = as.numeric(bigcv$config$mtrs[bestans]), ntree = as.numeric(bigcv$config$ntrees[bestans]), nodesize = as.numeric(bigcv$config$nodesizes[bestans]), quiet = T, featsel = T, oset = clsoset)$result
  
  corprobs <- probcorrect(cvpreds, seed = seedno, type = "crossval", levs = levels(as.factor(traindata$Type)))
  
  A <- corprobs$A
  
  B <- corprobs$B
  
  print(paste("Calibration complete. Optimal values: A = ", A, " and B = ", B, sep = ""))
  
  print("Calibrating the probabilities of the aliased classification model using the Brier-Score of a 5-fold cross validation on the training set.")
  
  cvapreds <- crossval(trainadata[,c(2, 10:ncol(trainadata))], seed = seedno, k=5, cvpred = T, mtr = as.numeric(bigacv$config$mtrs[bestaans]), ntree = as.numeric(bigacv$config$ntrees[bestaans]), nodesize = as.numeric(bigacv$config$nodesizes[bestaans]), quiet = T, featsel = T)$result
  
  coraprobs <- probcorrect(cvapreds, seed = seedno, type = "crossval", levs = levels(as.factor(trainadata$Type)))
  
  Aa <- coraprobs$A
  
  Ba <- coraprobs$B
  
  print(paste("Calibration complete. Optimal values: A = ", Aa, " and B = ", Ba, sep = ""))
  
  print("Classifying the current variable star set using the classification model...")
  
  varprobs <- classgen(vardata, classmodel, traindata, showprobs = T)
  
  varaprobs <- classgen(vardata, classamodel, trainadata, showprobs = T)
  
  print("Operation complete.")
  
  print("Calibrating probabilities using the calibration values...")
  
  varprobs <- probcorother(varprobs, A = A, B = B, cutoffs = clscutoffs, type = "class")
  
  varcorcuts <- varprobs$corcutoffs
  
  varprobs <- varprobs$corprobs
  
  varprobs <- selcutset(varprobs, cutoff = varcorcuts, clsnm = levels(as.factor(traindata$Type)), doprobs = FALSE)
  
  varprobs2 <- selcutset(varprobs, cutoff = varcorcuts, clsnm = levels(as.factor(traindata$Type)), doprobs = TRUE)
  
  varaprobs <- probcorother(varaprobs, A = Aa, B = Ba, cutoffs = clsacutoffs, type = "class")
  
  varacorcuts <- varaprobs$corcutoffs
  
  varaprobs <- varaprobs$corprobs
  
  varaprobs <- selcutset(varaprobs, cutoff = varacorcuts, clsnm = levels(as.factor(trainadata$Type)), doprobs = FALSE)
  
  varaprobs2 <- selcutset(varaprobs, cutoff = varacorcuts, clsnm = levels(as.factor(trainadata$Type)), doprobs = TRUE)
  
  print("Operation complete.")
  
  print("Validating the classifications and selecting high-confidence objects to search for additional variability...")
  
  #varcands <- classvalid(varprobs, cut = c(0.8, 5))
  
  #varacands <- classvalid(varaprobs, cut = c(0.8, 5))
  
  varcands <- varprobs[which(varprobs$Name %in% traindata$Name),]
  
  varacands <- varaprobs[which(varaprobs$Name %in% trainadata$Name),]
  
  print("Selecting Variability Indices for these variability candidates...")
  
  varfcands <- unique(rbind(varcands, varacands))
  
  if (onlymatch == TRUE){
    
    varfcands <- varcands
    
  }
  
  varcandsind <- objdata[which(objdata$USnoref %in% varfcands$USnoref),]
  
  varcandsind <- varcandsind[order(varcandsind$Name),]
  
  rownames(varcandsind) <- 1:nrow(varcandsind)
  
  varcandsind$Type <- as.character(as.vector(varcandsind$Type))
  
  varcandsind$Type <- "VAR"
  
  print("Combining these indices with a set of randomly selected objects to act as non-variable candidates...")
  
  nvarcandsind <- objdata[-(which(objdata$USnoref %in% varcands$USnoref)),]
  
  nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]
  
  nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]
  
  rownames(nvarcandsind) <- 1:nrow(nvarcandsind)
  
  nvarcandsind$Type <- as.character(as.vector(nvarcandsind$Type))
  
  nvarcandsind$Type <- "NVAR"
  
  vardetset <- rbind(varcandsind, nvarcandsind)
  
  rownames(vardetset) <- 1:nrow(vardetset)
  
  print("Beginning the training of the variability detection model...")
  
  print("Optimising the model hyperparameters using a 5-fold cross validation with 2 repeats.")
  
  print(paste("'", optfn, "' selected as the optimisation statistic.", sep = ""))
  
  bigvarcv <- bigparamsel(vardetset[,c(2, 10:ncol(vardetset))], optfn = optfn, mtrs = seq(2, 20, by = 2), ntrees = c(1000, 1500, 2000), oset = varoset)
  
  bigvarcvstats <- bigparamselmean(bigvarcv$result)
  
  print("Fine-tuning the hyperparameters around the optimal values...")
  
  if (optfn == "BS"){
    
    bestvarstat <- min(bigvarcvstats$meanvals)
    
    bestvarans <- which.min(bigvarcvstats$meanvals)
    
  }
  else{
    
    bestvarstat <- max(bigvarcvstats$meanvals)
    
    bestvarans <- which.max(bigvarcvstats$meanvals)
    
  }
  
  bigvarcv <- bigparamsel(vardetset[,c(2, 10:ncol(vardetset))], optfn = optfn, oset = varoset, mtrs = seq(as.numeric(bigvarcv$config$mtrs[bestvarans])-2, as.numeric(bigvarcv$config$mtrs[bestvarans])+2, by = 1), ntrees = c(250, 500, 1000, 2000), nodesizes = seq(as.numeric(bigvarcv$config$nodesizes[bestvarans])-2, as.numeric(bigvarcv$config$nodesizes[bestvarans])+2, by = 1))
  
  bigvarcvstats <- bigparamselmean(bigvarcv$result)
  
  print("Training final model using the optimum hyperparameters...")
  
  if (optfn == "BS"){
    
    bestvarstat <- min(bigvarcvstats$meanvals)
    
    bestvarans <- which.min(bigvarcvstats$meanvals)
    
  }
  else{
    
    bestvarstat <- max(bigvarcvstats$meanvals)
    
    bestvarans <- which.max(bigvarcvstats$meanvals)
    
  }
  
  print(paste("Best hyperparameters located with optimisation statistic of ", bestvarstat, ".", sep = ""))
  
  print(paste("These are: mtr = ", bigvarcv$config$mtrs[bestvarans], ", ntrees = ", bigvarcv$config$ntrees[bestvarans], ", and nodesize = ", bigvarcv$config$nodesizes[bestvarans], ".", sep = ""))
  
  varcutoffs <- as.vector(bigvarcv$cutoffs[[bestvarans]])
  
  print(varcutoffs)
  
  varmodel <- crossval(vardetset[,c(2, 10:ncol(vardetset))], train = T, mtr = as.numeric(bigvarcv$config$mtrs[bestvarans]), ntree = as.numeric(bigvarcv$config$ntrees[bestvarans]), nodesize = as.numeric(bigvarcv$config$nodesizes[bestvarans]), quiet = T, featsel = T, oset = varoset)
  
  print("Training of variability detection model complete.")
  
  #print("Calibrating the probabilities of this model using the Brier-Score of a 5-fold cross validation on the training set.")
  
  #cvdetpreds <- crossval(vardetset[,c(2, 10:ncol(vardetset))], seed = seedno, k=5, cvpred = T, mtr = as.numeric(bigvarcv$config$mtrs[bestvarans]), ntree = as.numeric(bigvarcv$config$ntrees[bestvarans]), nodesize = as.numeric(bigvarcv$config$nodesizes[bestvarans]), quiet = T, featsel = T, oset = varoset)
  
  #cordetprobs <- probcorrect(cvdetpreds, seed = seedno, type = "crossval", levs = levels(as.factor(vardetset$Type)))
  
  #Av <- cordetprobs$A
  
  #Bv <- cordetprobs$B
  
  #cordetprobs <- cordetprobs$corprobs
  
  #print(paste("Calibration complete. Optimal values: A = ", Av, " and B = ", Bv, sep = ""))
  
  #print("Using this model to detect additional variable candidates in the dataset...")
  
  vardetprobs <- vargen(objdata, varmodel, vardetset, cutoffs = varcutoffs, showprobs = T)
  
  print("Operation complete.")
  
  #print("Calibrating probabilities using the calibration values...")
  
  #vardetprobs <- probcorother(vardetprobs, A = Av, B = Bv, type = "vardet")
  
  #print("Operation complete.")
  
  print("Identifying new candidate variables from the probability of being a variable...")
  
  newvars <- vardetprobs[which(!(vardetprobs$USnoref %in% vardata$USnoref)),]
  
  newvars <- newvars[order(newvars$VAR.prob, decreasing = T),]
  
  print(paste("There are ", nrow(newvars[which(newvars$VAR.prob >= varcutoffs[2]),]), " variables in this selection.", sep = ""))
  
  newvarcands <- newvars[which(newvars$VAR.prob >= varcutoffs[2]),]
  
  #print("Generating full feature set for these candidate variables...")
  
  #newdata <- featgen(newvarcands, radii = 0.042, quiet = T, gaufil = gaufil, LPVcut = c(0.1, 0.7, 0.7), sigk = 4, nonper = F)
  
  #newdata <- vardata[complete.cases(newdata[,-3]),]
  
  #print("Classifying the candidate variables...")
  
  #newprobs <- classgen(newdata, classmodel, traindata, showprobs = T)
  
  #print("Operation complete.")
  
  #print("Calibrating probabilities using the calibration values...")
  
  #newprobs <- probcorother(newprobs, A = A, B = B, type = "class")
  
  #print("Operation complete.")
  
  #print("Validating the classifications and selecting high-confidence objects to search for additional variability...")
  
  #newcands <- classvalid(newprobs, cut = c(0.8, 5))
  
  out = list(newvarcands = newvarcands, varprobs = varprobs, varaprobs = varaprobs, varcands = varcands, varacands = varacands, varfcands = varfcands, varcandsind = varcandsind, objdata = objdata, traindata = traindata, trainadata = trainadata, classmodel = classmodel, classamodel = classamodel, varmodel = varmodel, A = A, B = B, Aa = Aa, Ba = Ba, clscutoffs = clscutoffs, clsacutoffs = clsacutoffs, varcutoffs = varcutoffs, varcorcuts = varcorcuts, varacorcuts = varacorcuts, varprobs2 = varprobs2, varaprobs2 = varaprobs2, bigcv = bigcv, bigcvstats = bigcvstats, bigacv = bigacv, bigacvstats = bigacvstats, bigvarcv = bigvarcv, bigvarcvstats = bigvarcvstats)
  
  out
  
}






pipeline2 <- function(varcat, lim = 100, ep = 0.05, seedno = 10, varinddata = NULL, fullcatdata = NULL, optfn = "AUC", usecolour = FALSE){
  
  library(RODBC)
  
  set.seed(seedno)
  
  print("Executing Skycam Automated Variability Classification Pipeline.")
  
  gaufil <- c(0, 100)
  
  # Firstly we want to generate a full set of variability indices for every object in the database
  # with more than lim observations. Only if varinddata doesnt exist.
  
  if (is.null(varinddata)){
    
    objdatset <- sellimobj(lim)
    
    print("Calculating Variability Indices for these objects... **This may take some time**")
    
    print("This can be overridden using the varinddata function argument for already generated data.")
    
    objdata <- featgen(objdatset, radii = 0.042, quiet = T, gaufil = gaufil, LPVcut = c(0.1, 0.7, 0.7), sigk = 4, nonper = T)
    
    objdata <- objdata[complete.cases(objdata[,-3]),]
    
  }
  else{
    
    print("Variability Indices supplied for a set of Skycam-T objects. Using these...")
    
    n <- nrow(varinddata)
    
    objdatset <- as.data.frame(cbind(as.character(as.vector(varinddata$USnoref)), rep("--", n), rep("NEW", n), as.numeric(as.vector(varinddata$Right.Ascension)), as.numeric(as.vector(varinddata$Declination)), rep("", n), rep("--", n)))
    
    colnames(objdatset) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude")
    
    print("Data Imported.")
    
    objdata <- varinddata[complete.cases(varinddata[,-3]),]
    
  }
  
  if (usecolour == TRUE){
  
    if (any(is.na(objdata$B.R.Colour)) | any(objdata$B.R.Colour == 1000)){
    
      print("Missing B-R Colour information. Imputing missing values using missForest...")
    
      objdata <- colimpute(objdata, maxiter = 10, ntree = 100)
    
      print("Imputation Complete.")
    
    }
    else{
    
      print("No Missing data for B-R Colour. Skipping missForest imputation procedure.")
    
    }
    
  }
  else{
    
    objdata <- objdata[,-which(names(objdata) %in% "B.R.Colour")]
    
  }
  
  print("Using catalogued Variable Stars to construct classification calibration set...")
  
  if (is.null(fullcatdata)){
    
    print("Calculating full feature set for catalogue objects... **This may take some time**")
    
    print("This can be overridden using the fullcatdata function argument for already generated data.")
    
    vardata <- featgen(varcat, radii = 0.042, quiet = T, gaufil = gaufil, LPVcut = c(0.1, 0.7, 0.7), sigk = 4, nonper = F)
    
    vardata <- vardata[complete.cases(vardata),]
    
  }
  else{
    
    print("Classification Caibration set full features supplied for a set of Skycam-T objects. Using these...")
    
    n <- nrow(fullcatdata)
    
    varcat <- as.data.frame(cbind(as.character(as.vector(fullcatdata$Name)), rep("--", n), as.character(as.vector(fullcatdata$Type)), as.numeric(as.vector(fullcatdata$Right.Ascension)), as.numeric(as.vector(fullcatdata$Declination)), as.numeric(as.vector(fullcatdata$AAVSO.Period)), as.numeric(as.vector(fullcatdata$Median.Magnitude))))
    
    colnames(objdatset) <- c("Name", "AUID", "Type", "RA", "DEC", "Period", "Magnitude")
    
    print("Data Imported.")
    
    vardata <- fullcatdata[complete.cases(fullcatdata),]
    
  }
  
  if (usecolour == TRUE){
  
    print("Calibrating the B-R Colours to the Skycam-T variability index data...")
  
    for (i in 1:nrow(vardata)){
    
      vardata$B.R.Colour[i] <- objdata$B.R.Colour[which(objdata$USnoref == vardata$USnoref[i])]
    
    }
  
    print("Colour calibration complete.")
  
  }
  
  print(paste(nrow(vardata), " classification calibration objects present.", sep = ""))
  
  print(paste("These objects cover ", length(levels(as.factor(vardata$Type))), " known variability classes.", sep = ""))
  
  print("Selecting catalogue matching objects using the period feature (allowing for Hits, Multiples and Aliases)...")
  
  n.old <- nrow(vardata)
  
  permod <- periodmodes(vardata, ep = ep, per = T)
  
  traindata <- vardata[which(permod$Unknown == 0),]
  
  n <- nrow(traindata)
  
  print(paste("This filter results in ", n, " classification training objects with ", n.old-n, " removed objects.", sep = ""))
  
  out = list(objdata = objdata, traindata = traindata)
  
  out
  
}