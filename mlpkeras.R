mlpkeras <- function(seedno = 20){
  
  library(keras)
  
  set.seed(seedno)
  
  # Data Preparation ---------------------------------------------------
  
  batch_size <- 128
  num_classes <- 13
  epochs <- 4000
  
  data <- mlptraindata[sample(1:nrow(mlptraindata), nrow(mlptraindata)),c(2, 10:ncol(mlptraindata))]
  
  data[,-1] <- scale(data[,-1])
  
  x_train <- as.matrix(data[1:500,2:ncol(data)])
  y_train <- as.numeric(data[1:500,1])-1
  x_test <- as.matrix(data[501:866,2:ncol(data)])
  y_test <- as.numeric(data[501:866,1])-1
  
  y_act <- as.factor(y_test)
  
  cat(nrow(x_train), 'train samples\n')
  cat(nrow(x_test), 'test samples\n')
  
  # Convert class vectors to binary class matrices
  y_train <- to_categorical(y_train, num_classes)
  y_test <- to_categorical(y_test, num_classes)
  
  # Define Model --------------------------------------------------------------
  
  model <- keras_model_sequential()
  model %>% 
    layer_dense(units = 1024, activation = 'relu', input_shape = length(2:ncol(data))) %>% 
    layer_dropout(rate = 0.4) %>% 
    layer_dense(units = 1024, activation = 'relu') %>%
    layer_dropout(rate = 0.3) %>%
    layer_dense(units = 256, activation = 'relu') %>%
    layer_dropout(rate = 0.2) %>%
    layer_dense(units = 13, activation = 'softmax')
  
  summary(model)
  
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
  
  plot(history)
  
  score <- model %>% evaluate(
    x_test, y_test,
    verbose = 0
  )
  
  # Output metrics
  cat('Test loss:', score[[1]], '\n')
  cat('Test accuracy:', score[[2]], '\n')
  
  y_test_hat <- as.factor(predict_classes(model, x_test))
  
  levels(y_act) <- levels(mlptraindata$Type)
  
  levels(y_test_hat) <- levels(mlptraindata$Type)
  
  print(confusionMatrix(y_test_hat, y_act))
  print(mean(as.numeric(y_act) == as.numeric(y_test_hat)))
  
  fm <- try(as.numeric(as.vector(confusionMatrix(y_test_hat, y_act)$byClass[,7])), TRUE)
  
  if (class(fm) == "try-error"){
    
    fm <- as.numeric(as.vector(confusionMatrix(y_test_hat, y_act)$byClass[7]))
    
  }
  
  fm[which(is.na(fm) | is.nan(fm))] <- 0
  
  print(mean(fm))
  
}