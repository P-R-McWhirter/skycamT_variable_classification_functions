run_analysis <- function(outputfile) {
  
  ## The Samsung data should be extracted into the working directory using the same structure as it has
  ## in the zipped folder. The folder UCI HAR Dataset should be visible in the working directory.
  
  library(reshape2)
  
  ## First we will grab the testing and training dataset feature vectors vs observations
  ## data from the appropriate folders as if we had extracted the zip file into our working directory.
  
  trainData <- read.table("./UCI HAR Dataset/train/X_train.txt")
  testData <- read.table("./UCI HAR Dataset/test/x_test.txt")
  
  ## Let's also read in the class information and the subject number for each feature vector.
  
  trainClass <- read.table("./UCI HAR Dataset/train/y_train.txt")
  testClass <- read.table("./UCI HAR Dataset/test/y_test.txt")
  trainSub <- read.table("./UCI HAR Dataset/train/subject_train.txt")
  testSub <- read.table("./UCI HAR Dataset/test/subject_test.txt")
  
  ## Next we will create a dataset containing all the observations from both input datasets.
  
  data <- rbind(trainData, testData)
  
  ## For future work, we will also do the same for the class and subject datasets.
  
  class <- rbind(trainClass, testClass)
  subject <- rbind(trainSub, testSub)
  
  ## Enter column names listed in the features text file into a dataframe.
  ## This describes the variables used in the dataset.
  
  columnNames <- read.table("./UCI HAR Dataset/features.txt")
  
  ## We also grab the activity names from the associated text file and load it into a dataframe.
  
  activityNames <- read.table("./UCI HAR Dataset/activity_labels.txt")
  
  ## Now let's add these names to the inputData dataframe.
  
  names(data) <- columnNames[,2]
  
  ## We want all variables relating to the Mean and Standard Deviation of the measurements.
  ## We can search for partial matches to 'mean()' and 'std()' by using grepl and placing the
  ## desired values into a character vector.
  
  desired = c("mean()", "std()")
  
  data <- data[,colnames(data)[grepl(paste(desired,collapse="|"),colnames(data))]]
  
  ## Finally, now that we have the important variables selected, we can add our activity 'class' information
  ## to the 'data' data frame as an activity column. The same is repeated with the subject information.
  
  data$Activity <- class[,1]
  data$Subject <- subject[,1]
  
  ## Next, we define a character vector containing the different activities currently represented by numbers.
  
  activ = as.vector(activityNames[,2])
  
  ## Finally, this for loop will replace the numbers representing the activities with the names
  ## presented in the above character vector.
  
  for (i in 1:length(activ)){
  
  data$Activity[data$Activity == i] <- activ[i]
  
  }
  
  ## Next we want to melt the dataframe by Subject and Activity.
  
  dataMelt <- melt(data, id.vars=c("Subject","Activity"), variable.name = "Measurement_Name", value.name = "Value")
  
  ## Now we cast the molten dataframe to produce a dataframe with the Subjects and their activities against the
  ## aggregated mean value of each of the variables producing a wide tidy dataset.
  
  finalData <- dcast(dataMelt, Subject + Activity ~ ..., fun.aggregate = mean)
  
  ## Finally, we will write this table to a file whose name is the argument of the function.
  
  write.table(finalData, file = outputfile, row.names = FALSE)
  
}