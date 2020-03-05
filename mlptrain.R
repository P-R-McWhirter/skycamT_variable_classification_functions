mlptrain <- function(train, test, hidden = 200) {
  
  require(mxnet)
  
  train <- data.matrix(train)
  test <- data.matrix(test)
  
  train[,1] <- as.numeric(train[,1])
  test[,1] <- as.numeric(test[,1])
  
  train.x <- train[,-1]
  train.y <- train[,1]
  
  test.x <- test[,-1]
  test.y <- test[,1]
  
  mx.set.seed(100)
  
  model <- mx.mlp(train.x, train.y, hidden_node=hidden, out_node=5, activation = "tanh", out_activation="softmax",
                  num.round=10000, array.batch.size=1000, learning.rate=0.005, momentum=0.9, initializer=mx.init.uniform(0.07),
                  eval.metric=mx.metric.accuracy, epoch.end.callback = mx.callback.log.train.metric(100))
  
  model
  
}