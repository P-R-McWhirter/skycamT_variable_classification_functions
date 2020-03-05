cnntrend <- function(data){
  
  library(keras)
  
  set.seed(20)
  
  # Data Preparation -----------------------------------------------------
  
  batch_size <- 128
  num_classes <- length(levels(data[,1]))
  epochs <- 500
  
  # Input image dimensions
  img_rows <- 1
  img_cols <- 100
  
  # The data, shuffled and split between train and test sets
  
  data <- data[sample(1:nrow(data), nrow(data)),]
  
  trainsize <- floor(0.7*nrow(data))
  
  x_train <- as.matrix(data[1:trainsize,2:101])
  y_train <- as.numeric(data[1:trainsize,1])-1
  x_test <- as.matrix(data[(trainsize+1):nrow(data),2:101])
  y_test <- as.numeric(data[(trainsize+1):nrow(data),1])-1
  
  y_act <- as.factor(y_test)
  
  # Redefine  dimension of train/test inputs
  x_train <- array_reshape(x_train, c(nrow(x_train), img_rows, img_cols, 1))
  x_test <- array_reshape(x_test, c(nrow(x_test), img_rows, img_cols, 1))
  input_shape <- c(img_rows, img_cols, 1)
  
  # Transform RGB values into [0,1] range
  #x_train <- x_train / 255
  #x_test <- x_test / 255
  
  cat('x_train_shape:', dim(x_train), '\n')
  cat(nrow(x_train), 'train samples\n')
  cat(nrow(x_test), 'test samples\n')
  
  # Convert class vectors to binary class matrices
  y_train <- to_categorical(y_train, num_classes)
  y_test <- to_categorical(y_test, num_classes)
  
  
  # Define Model -----------------------------------------------------------
  
  model <- keras_model_sequential()
  model %>%
    layer_conv_2d(filters = 192, kernel_size = c(1,3), activation = 'relu', strides = c(1,1),
                  input_shape = input_shape) %>% 
    #layer_max_pooling_2d(pool_size = c(1, 2)) %>% 
    layer_conv_2d(filters = 512, kernel_size = c(1,2), activation = 'relu', strides = c(1,1)) %>% 
    #layer_max_pooling_2d(pool_size = c(1, 2)) %>% 
    layer_conv_2d(filters = 400, kernel_size = c(1,2), activation = 'relu', strides = c(1,1)) %>%
    layer_dropout(rate = 0.25) %>% 
    layer_flatten(input_shape = input_shape) %>% 
    layer_dense(units = 256, activation = 'relu') %>% 
    layer_dropout(rate = 0.5) %>% 
    layer_dense(units = num_classes, activation = 'softmax')
  
  print(summary(model))
  
  # Compile model
  model %>% compile(
    loss = loss_categorical_crossentropy,
    optimizer = optimizer_adadelta(),
    metrics = c('accuracy')
  )
  
  # Train & Evaluate -------------------------------------------------------
  
  model %>% fit(
    x_train, y_train,
    batch_size = batch_size,
    epochs = epochs,
    verbose = 1,
    validation_data = list(x_test, y_test)
  )
  scores <- model %>% evaluate(
    x_test, y_test, verbose = 0
  )
  
  # Output metrics
  cat('Test loss:', scores[[1]], '\n')
  cat('Test accuracy:', scores[[2]], '\n')
  
  y_test_hat <- as.factor(predict_classes(model, x_test))
  
  print(confusionMatrix(y_test_hat, y_act))
  print(mean(as.numeric(y_act) == as.numeric(y_test_hat)))
  
}




fnntrend <- function(data){
  
  library(keras)
  
  set.seed(20)
  
  # Data Preparation -----------------------------------------------------
  
  batch_size <- 128
  num_classes <- length(levels(data[,1]))
  epochs <- 10000
  
  # Input image dimensions
  img_rows <- 1
  img_cols <- ncol(data)-1
  
  # The data, shuffled and split between train and test sets
  
  data <- data[sample(1:nrow(data), nrow(data)),]
  
  trainsize <- floor(0.7*nrow(data))
  
  x_train <- as.matrix(data[1:trainsize,-1])
  y_train <- as.numeric(data[1:trainsize,1])-1
  x_test <- as.matrix(data[(trainsize+1):nrow(data),-1])
  y_test <- as.numeric(data[(trainsize+1):nrow(data),1])-1
  
  y_act <- as.factor(y_test)
  
  # Redefine  dimension of train/test inputs
  # reshape
  x_train <- array_reshape(x_train, c(nrow(x_train), img_cols))
  x_test <- array_reshape(x_test, c(nrow(x_test), img_cols))
  # rescale
  x_train <- x_train / 255
  x_test <- x_test / 255
  
  # Transform RGB values into [0,1] range
  #x_train <- x_train / 255
  #x_test <- x_test / 255
  
  cat('x_train_shape:', dim(x_train), '\n')
  cat(nrow(x_train), 'train samples\n')
  cat(nrow(x_test), 'test samples\n')
  
  # Convert class vectors to binary class matrices
  y_train <- to_categorical(y_train, num_classes)
  y_test <- to_categorical(y_test, num_classes)
  # Define Model -----------------------------------------------------------
  
  model <- keras_model_sequential() 
  model %>% 
    layer_dense(units = 120, activation = 'relu', input_shape = c(img_cols)) %>% 
    layer_dropout(rate = 0.5) %>% 
    layer_dense(units = 20, activation = 'relu') %>% 
    layer_dropout(rate = 0.25) %>% 
    layer_dense(units = num_classes, activation = 'softmax')
  
  print(summary(model))
  
  # Compile model
  model %>% compile(
    loss = loss_categorical_crossentropy,
    optimizer = optimizer_adadelta(),
    metrics = c('accuracy')
  )
  
  # Train & Evaluate -------------------------------------------------------
  
  model %>% fit(
    x_train, y_train,
    batch_size = batch_size,
    epochs = epochs,
    verbose = 0,
    validation_data = list(x_test, y_test)
  )
  scores <- model %>% evaluate(
    x_test, y_test, verbose = 0
  )
  
  # Output metrics
  cat('Test loss:', scores[[1]], '\n')
  cat('Test accuracy:', scores[[2]], '\n')
  
  y_test_hat <- as.factor(predict_classes(model, x_test))
  
  print(confusionMatrix(y_test_hat, y_act))
  print(mean(as.numeric(y_act) == as.numeric(y_test_hat)))
  
}