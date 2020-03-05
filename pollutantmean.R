pollutantmean <- function(directory, pollutant, id = 1:332) {

  files_full <- list.files(directory, full.names=TRUE)
  
  dat <- data.frame()
  for (i in id) {
    dat <- rbind(dat, read.csv(files_full[i]))
  }
  
  test <- dat[,pollutant]
  
  good <- complete.cases(test)
  
  test2 <- test[good]
  
  mean(test2)
}