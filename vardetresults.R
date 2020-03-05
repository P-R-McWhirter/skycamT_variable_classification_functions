vardetresults$f1 <- crossval(vardetset[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f2 <- crossval(vardetset2[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f3 <- crossval(vardetset3[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f4 <- crossval(vardetset4[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f5 <- crossval(vardetset5[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f6 <- crossval(vardetset6[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f7 <- crossval(vardetset7[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f8 <- crossval(vardetset8[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f9 <- crossval(vardetset9[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f10 <- crossval(vardetset10[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f11 <- crossval(vardetset11[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f12 <- crossval(vardetset12[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f13 <- crossval(vardetset13[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f14 <- crossval(vardetset14[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f15 <- crossval(vardetset15[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f16 <- crossval(vardetset16[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f17 <- crossval(vardetset17[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f18 <- crossval(vardetset18[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f19 <- crossval(vardetset19[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

vardetresults$f20 <- crossval(vardetset20[,c(2, 10:44)], k = 5, r = 2, seed = c(100, 1000), alg = "RF", mtr = 4, ntree = 1500, nodesize = 25, quiet = F, featsel = T, output = "F1", hid = 100, oset = 0, scaleit = F, train = F)$result

