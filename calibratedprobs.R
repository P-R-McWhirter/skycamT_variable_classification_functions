plot(precaldata$probvortestrfprobs, precaldata$probvortesttrueprobs, type = "n", pch = 19, lwd = 3, main = "Calibration of RF Probabilities on aliased data", ylab = "True Class Probability", xlab = "RF Estimated Probability", ylim = c(0, 1), xlim = c(0, 1))

abline(0, 1, lty = 1, lwd = 2, col = "darkgray")

lines(precaldata$probvortestrfprobs, precaldata$probvortesttrueprobs, type = "o", pch = 19, col = "black", lwd = 3, lty = 1)

lines(precaldata$probvortestrfprobs, probvortesttrueprobs, type = "o", pch = 19, col = "blue", lwd = 3, lty = 2)

segments(probcortestsplits[1:13], probvortesttrueprobs, probcortestsplits[2:14], probvortesttrueprobs, col = "blue")

segments(probcortestsplits[1:13], precaldata$probvortesttrueprobs, precaldata$probcortestsplits[2:14], precaldata$probvortesttrueprobs)

segments(precaldata$probvortestrfprobs,  precaldata$probvortesttrueprobs+probcortesttrueerror2, precaldata$probvortestrfprobs, precaldata$probvortesttrueprobs-probcortesttrueerror2)

segments(precaldata$probvortestrfprobs,  probvortesttrueprobs+probcortesttrueerror, precaldata$probvortestrfprobs, probvortesttrueprobs-probcortesttrueerror, col = "blue")

legend("bottomright", inset = .05, cex = 1.5, title = "Probabilities", c("Original", "Adjusted"), lty = c(1,2), lwd = c(4,4), col = c("black", "blue"), pch = 19, bg = "grey96")