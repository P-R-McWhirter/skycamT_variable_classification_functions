neuralnet <- function(traindata, testdata, hid = 200, batch = 128, epochs = 200, seedno = 20){
  
  library(keras)
  
  set.seed(seedno)
  
  # Data Preparation ---------------------------------------------------
  
  batch_size <- batch
  num_classes <- length(levels(as.factor(traindata[,1])))
  
  #data[,-1] <- scale(data[,-1])
  
  x_train <- as.matrix(traindata[,2:ncol(traindata)])
  y_train <- as.numeric(traindata[,1])-1
  x_test <- as.matrix(testdata[,2:ncol(testdata)])
  y_test <- as.numeric(testdata[,1])-1
  
  y_act <- as.factor(y_test)
  
  cat(nrow(x_train), 'train samples\n')
  cat(nrow(x_test), 'test samples\n')
  
  # Convert class vectors to binary class matrices
  y_train <- to_categorical(y_train, num_classes)
  y_test <- to_categorical(y_test, num_classes)
  
  # Define Model --------------------------------------------------------------
  
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = hid, activation = 'relu', input_shape = length(2:ncol(traindata))) %>% 
    layer_dropout(rate = 0.2) %>% 
    layer_dense(units = num_classes, activation = 'softmax')
  
  #summary(model)
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  # Training & Evaluation ----------------------------------------------------
  
  # Fit model to data
  history <- model %>% fit(
    x_train, y_train,
    batch_size = batch_size,
    epochs = epochs,
    verbose = 0,
    validation_split = 0.2
  )
  
  #plot(history)
  
  score <- model %>% evaluate(
    x_test, y_test,
    verbose = 0
  )
  
  # Output metrics
  cat('Test loss:', score[[1]], '\n')
  cat('Test accuracy:', score[[2]], '\n')
  
  y_test_set <- predict_proba(model, x_test)
  
  y_test_hat <- rep(0, nrow(y_test_set))
  
  for (i in 1:nrow(y_test_set)){
    
    y_test_hat[i] <- which.max(y_test_set[i,])
    
  }
  
  y_test_hat <- as.factor(y_test_hat)
  
  levels(y_act) <- levels(traindata$Type)
  
  levels(y_test_hat) <- levels(traindata$Type)
  
  print(confusionMatrix(y_test_hat, y_act))
  
  fm <- try(as.numeric(as.vector(confusionMatrix(y_test_hat, y_act)$byClass[,7])), TRUE)
  
  if (class(fm) == "try-error"){
    
    fm <- as.numeric(as.vector(confusionMatrix(y_test_hat, y_act)$byClass[7]))
    
  }
  
  fm[which(is.na(fm) | is.nan(fm))] <- 0
  
  out <- list(model = model, x_test = x_test, y_act = y_act, y_test_hat = y_test_hat, fm = fm)
  
}







deepneuralnet <- function(traindata, testdata, hid = c(200, 200, 200), batch = 128, epochs = 40, seedno = 20){
  
  library(keras)
  
  set.seed(seedno)
  
  # Data Preparation ---------------------------------------------------
  
  batch_size <- batch
  num_classes <- length(levels(as.factor(traindata[,1])))
  
  print(num_classes)
  
  #data[,-1] <- scale(data[,-1])
  
  x_train <- as.matrix(traindata[,2:ncol(traindata)])
  y_train <- as.numeric(traindata[,1])-1
  x_test <- as.matrix(testdata[,2:ncol(testdata)])
  y_test <- as.numeric(testdata[,1])-1
  
  y_act <- as.factor(y_test)
  
  cat(nrow(x_train), 'train samples\n')
  cat(nrow(x_test), 'test samples\n')
  
  # Convert class vectors to binary class matrices
  y_train <- to_categorical(y_train, num_classes)
  y_test <- to_categorical(y_test, num_classes)
  
  # Define Model --------------------------------------------------------------
  
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = hid[1], activation = 'relu', input_shape = length(2:ncol(traindata))) %>% 
    layer_dropout(rate = 0.2) %>% 
    layer_dense(units = hid[2], activation = 'relu') %>% 
    layer_dropout(rate = 0.2) %>% 
    layer_dense(units = hid[3], activation = 'relu') %>% 
    layer_dropout(rate = 0.2) %>% 
    layer_dense(units = num_classes, activation = 'softmax')
  
  #summary(model)
  
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = optimizer_rmsprop(),
    metrics = c('accuracy')
  )
  
  # Training & Evaluation ----------------------------------------------------
  
  # Fit model to data
  history <- model %>% fit(
    x_train, y_train,
    batch_size = batch_size,
    epochs = epochs,
    verbose = 0,
    validation_split = 0.2
  )
  
  #plot(history)
  
  score <- model %>% evaluate(
    x_test, y_test,
    verbose = 0
  )
  
  # Output metrics
  cat('Test loss:', score[[1]], '\n')
  cat('Test accuracy:', score[[2]], '\n')
  
  y_test_set <- predict_proba(model, x_test)
  
  y_test_hat <- rep(0, nrow(y_test_set))
  
  for (i in 1:nrow(y_test_set)){
    
    y_test_hat[i] <- which.max(y_test_set[i,])+1
    
  }
  
  y_test_hat <- as.factor(y_test_hat)
  
  levels(y_act) <- levels(traindata$Type)
  
  levels(y_test_hat) <- levels(traindata$Type)
  
  print(confusionMatrix(y_test_hat, y_act))
  
  fm <- try(as.numeric(as.vector(confusionMatrix(y_test_hat, y_act)$byClass[,7])), TRUE)
  
  if (class(fm) == "try-error"){
    
    fm <- as.numeric(as.vector(confusionMatrix(y_test_hat, y_act)$byClass[7]))
    
  }
  
  fm[which(is.na(fm) | is.nan(fm))] <- 0
  
  out <- list(model = model, x_test = x_test, y_act = y_act, y_test_hat = y_test_hat, fm = fm)
  
}