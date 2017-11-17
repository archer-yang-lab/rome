huberreg <- function(r, delta) {
  0.5 * r^2 * (abs(r) <= delta) + delta * (abs(r) - delta/2) * (abs(r) > delta)
} 
