geometricMean = function(x) {
  exp(mean(log(x), na.rm = T))
}
