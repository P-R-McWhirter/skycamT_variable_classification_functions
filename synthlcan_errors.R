synthlcan_errors <- function(){
  
  out <- list()
  
  out$regcad_sin_pergm_snr2_20_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$regcad_saw_pergm_snr2_20_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$regcad_ecl_pergm_snr2_20_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$regcad_eeb_pergm_snr2_20_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$sktcad_sin_pergm_snr2_20_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$sktcad_saw_pergm_snr2_20_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$sktcad_ecl_pergm_snr2_20_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$sktcad_eeb_pergm_snr2_20_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$regcad_sin_pergm_snr2_30_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$regcad_saw_pergm_snr2_30_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$regcad_ecl_pergm_snr2_30_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$regcad_eeb_pergm_snr2_30_1 <- synthlcan_pergm(times_t[1:25,], testpers[1:25], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$sktcad_sin_pergm_snr2_30_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$sktcad_saw_pergm_snr2_30_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$sktcad_ecl_pergm_snr2_30_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$sktcad_eeb_pergm_snr2_30_1 <- synthlcan_pergm(times_s[1:25,], testpers[1:25], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$regcad_sin_pergm_snr2_20_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$regcad_saw_pergm_snr2_20_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$regcad_ecl_pergm_snr2_20_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$regcad_eeb_pergm_snr2_20_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$sktcad_sin_pergm_snr2_20_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$sktcad_saw_pergm_snr2_20_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$sktcad_ecl_pergm_snr2_20_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$sktcad_eeb_pergm_snr2_20_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$regcad_sin_pergm_snr2_30_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$regcad_saw_pergm_snr2_30_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$regcad_ecl_pergm_snr2_30_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$regcad_eeb_pergm_snr2_30_2 <- synthlcan_pergm(times_t[251:275,], testpers[251:275], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$sktcad_sin_pergm_snr2_30_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$sktcad_saw_pergm_snr2_30_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$sktcad_ecl_pergm_snr2_30_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$sktcad_eeb_pergm_snr2_30_2 <- synthlcan_pergm(times_s[251:275,], testpers[251:275], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$regcad_sin_pergm_snr2_20_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$regcad_saw_pergm_snr2_20_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$regcad_ecl_pergm_snr2_20_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$regcad_eeb_pergm_snr2_20_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$sktcad_sin_pergm_snr2_20_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$sktcad_saw_pergm_snr2_20_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$sktcad_ecl_pergm_snr2_20_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$sktcad_eeb_pergm_snr2_20_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$regcad_sin_pergm_snr2_30_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$regcad_saw_pergm_snr2_30_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$regcad_ecl_pergm_snr2_30_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$regcad_eeb_pergm_snr2_30_3 <- synthlcan_pergm(times_t[501:525,], testpers[501:525], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$sktcad_sin_pergm_snr2_30_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$sktcad_saw_pergm_snr2_30_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$sktcad_ecl_pergm_snr2_30_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$sktcad_eeb_pergm_snr2_30_3 <- synthlcan_pergm(times_s[501:525,], testpers[501:525], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$regcad_sin_pergm_snr2_20_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$regcad_saw_pergm_snr2_20_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$regcad_ecl_pergm_snr2_20_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$regcad_eeb_pergm_snr2_20_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$sktcad_sin_pergm_snr2_20_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "sin", redo = 5, lcseed = 20)
  
  out$sktcad_saw_pergm_snr2_20_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "saw", redo = 5, lcseed = 20)
  
  out$sktcad_ecl_pergm_snr2_20_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "ecl", redo = 5, lcseed = 20)
  
  out$sktcad_eeb_pergm_snr2_20_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "eeb", redo = 5, lcseed = 20)
  
  out$regcad_sin_pergm_snr2_30_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$regcad_saw_pergm_snr2_30_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$regcad_ecl_pergm_snr2_30_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$regcad_eeb_pergm_snr2_30_4 <- synthlcan_pergm(times_t[751:775,], testpers[751:775], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out$sktcad_sin_pergm_snr2_30_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "sin", redo = 5, lcseed = 30)
  
  out$sktcad_saw_pergm_snr2_30_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "saw", redo = 5, lcseed = 30)
  
  out$sktcad_ecl_pergm_snr2_30_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "ecl", redo = 5, lcseed = 30)
  
  out$sktcad_eeb_pergm_snr2_30_4 <- synthlcan_pergm(times_s[751:775,], testpers[751:775], snr = 2, type = "eeb", redo = 5, lcseed = 30)
  
  out
  
}




synthlcan_errors_grape <- function(){
  
  out <- list()
  
  #out$regcad_sin_pergm_snr2_30_1 <- synthlcan_grape(times_t[1:25,], testpers[1:25], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_saw_pergm_snr2_30_1 <- synthlcan_grape(times_t[1:25,], testpers[1:25], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_ecl_pergm_snr2_30_1 <- synthlcan_grape(times_t[1:25,], testpers[1:25], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_eeb_pergm_snr2_30_1 <- synthlcan_grape(times_t[1:25,], testpers[1:25], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_sin_pergm_snr2_30_1 <- synthlcan_grape(times_s[1:25,], testpers[1:25], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_saw_pergm_snr2_30_1 <- synthlcan_grape(times_s[1:25,], testpers[1:25], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_ecl_pergm_snr2_30_1 <- synthlcan_grape(times_s[1:25,], testpers[1:25], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_eeb_pergm_snr2_30_1 <- synthlcan_grape(times_s[1:25,], testpers[1:25], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_sin_pergm_snr2_30_2 <- synthlcan_grape(times_t[251:275,], testpers[251:275], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_saw_pergm_snr2_30_2 <- synthlcan_grape(times_t[251:275,], testpers[251:275], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_ecl_pergm_snr2_30_2 <- synthlcan_grape(times_t[251:275,], testpers[251:275], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_eeb_pergm_snr2_30_2 <- synthlcan_grape(times_t[251:275,], testpers[251:275], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_sin_pergm_snr2_30_2 <- synthlcan_grape(times_s[251:275,], testpers[251:275], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_saw_pergm_snr2_30_2 <- synthlcan_grape(times_s[251:275,], testpers[251:275], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_ecl_pergm_snr2_30_2 <- synthlcan_grape(times_s[251:275,], testpers[251:275], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_eeb_pergm_snr2_30_2 <- synthlcan_grape(times_s[251:275,], testpers[251:275], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_sin_pergm_snr2_30_3 <- synthlcan_grape(times_t[501:525,], testpers[501:525], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_saw_pergm_snr2_30_3 <- synthlcan_grape(times_t[501:525,], testpers[501:525], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_ecl_pergm_snr2_30_3 <- synthlcan_grape(times_t[501:525,], testpers[501:525], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_eeb_pergm_snr2_30_3 <- synthlcan_grape(times_t[501:525,], testpers[501:525], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_sin_pergm_snr2_30_3 <- synthlcan_grape(times_s[501:525,], testpers[501:525], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_saw_pergm_snr2_30_3 <- synthlcan_grape(times_s[501:525,], testpers[501:525], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_ecl_pergm_snr2_30_3 <- synthlcan_grape(times_s[501:525,], testpers[501:525], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_eeb_pergm_snr2_30_3 <- synthlcan_grape(times_s[501:525,], testpers[501:525], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_sin_pergm_snr2_30_4 <- synthlcan_grape(times_t[751:775,], testpers[751:775], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_saw_pergm_snr2_30_4 <- synthlcan_grape(times_t[751:775,], testpers[751:775], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_ecl_pergm_snr2_30_4 <- synthlcan_grape(times_t[751:775,], testpers[751:775], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  #out$regcad_eeb_pergm_snr2_30_4 <- synthlcan_grape(times_t[751:775,], testpers[751:775], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_sin_pergm_snr2_30_4 <- synthlcan_grape(times_s[751:775,], testpers[751:775], snr = 2, type = "sin", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_saw_pergm_snr2_30_4 <- synthlcan_grape(times_s[751:775,], testpers[751:775], snr = 2, type = "saw", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_ecl_pergm_snr2_30_4 <- synthlcan_grape(times_s[751:775,], testpers[751:775], snr = 2, type = "ecl", redo = 5, lcseed = 30, seedno = 1)
  
  out$sktcad_eeb_pergm_snr2_30_4 <- synthlcan_grape(times_s[751:775,], testpers[751:775], snr = 2, type = "eeb", redo = 5, lcseed = 30, seedno = 1)
  
  out
  
}