varcandsind <- pipelinerun$traindata[,1:44]

varcandsind$Type <- "VAR"

set.seed(10)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset <- rbind(varcandsind, nvarcandsind)

vartestmodel <- crossval(vardetset[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(20)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset2 <- rbind(varcandsind, nvarcandsind)

vartestmodel2 <- crossval(vardetset2[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(30)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset3 <- rbind(varcandsind, nvarcandsind)

vartestmodel3 <- crossval(vardetset3[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(40)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset4 <- rbind(varcandsind, nvarcandsind)

vartestmodel4 <- crossval(vardetset4[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(50)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset5 <- rbind(varcandsind, nvarcandsind)

vartestmodel5 <- crossval(vardetset5[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(60)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset6 <- rbind(varcandsind, nvarcandsind)

vartestmodel6 <- crossval(vardetset6[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(70)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset7 <- rbind(varcandsind, nvarcandsind)

vartestmodel7 <- crossval(vardetset7[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(80)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset8 <- rbind(varcandsind, nvarcandsind)

vartestmodel8 <- crossval(vardetset8[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(90)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset9 <- rbind(varcandsind, nvarcandsind)

vartestmodel9 <- crossval(vardetset9[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(100)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset10 <- rbind(varcandsind, nvarcandsind)

vartestmodel10 <- crossval(vardetset10[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(110)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset11 <- rbind(varcandsind, nvarcandsind)

vartestmodel11 <- crossval(vardetset11[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(120)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset12 <- rbind(varcandsind, nvarcandsind)

vartestmodel12 <- crossval(vardetset12[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(130)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset13 <- rbind(varcandsind, nvarcandsind)

vartestmodel13 <- crossval(vardetset13[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(140)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset14 <- rbind(varcandsind, nvarcandsind)

vartestmodel14 <- crossval(vardetset14[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(150)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset15 <- rbind(varcandsind, nvarcandsind)

vartestmodel15 <- crossval(vardetset15[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(160)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset16 <- rbind(varcandsind, nvarcandsind)

vartestmodel16 <- crossval(vardetset16[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(170)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset17 <- rbind(varcandsind, nvarcandsind)

vartestmodel17 <- crossval(vardetset17[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(180)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset18 <- rbind(varcandsind, nvarcandsind)

vartestmodel18 <- crossval(vardetset18[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(190)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset19 <- rbind(varcandsind, nvarcandsind)

vartestmodel19 <- crossval(vardetset19[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

set.seed(200)

nvarcandsind <- obj2dat[-(which(obj2dat$USnoref %in% varcandsind$USnoref)),]

nvarcandsind <- nvarcandsind[sample(1:nrow(nvarcandsind), min(c(nrow(varcandsind), nrow(nvarcandsind)))),]

nvarcandsind <- nvarcandsind[order(nvarcandsind$Name),]

nvarcandsind$Type <- "NVAR"

vardetset20 <- rbind(varcandsind, nvarcandsind)

vartestmodel20 <- crossval(vardetset20[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 200, oset = 0, scaleit = F, train = T)

vardetimp <- cbind(vartestmodel$importance[,4], vartestmodel2$importance[,4], vartestmodel3$importance[,4], vartestmodel4$importance[,4], vartestmodel5$importance[,4], vartestmodel6$importance[,4], vartestmodel7$importance[,4], vartestmodel8$importance[,4], 
                   vartestmodel9$importance[,4], vartestmodel10$importance[,4], vartestmodel11$importance[,4], vartestmodel12$importance[,4], vartestmodel13$importance[,4], vartestmodel14$importance[,4], vartestmodel15$importance[,4], vartestmodel16$importance[,4], 
                   vartestmodel17$importance[,4], vartestmodel18$importance[,4], vartestmodel19$importance[,4], vartestmodel20$importance[,4])

vardetmean <- rowMeans(vardetimp)

vardetsd <- rep(0, 35)

for (i in 1:35){
  
  vardetsd[i] <- sd(vardetimp[i,])
  
}

vardetstats <- cbind(vardetmean, vardetsd)