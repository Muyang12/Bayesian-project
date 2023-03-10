########## R script: addcircle ###########

# Adds circles for simulated datasets.

# Last changed: 28 JAN 2020

# Input
# c1,c2 circle centre
# r     radius
# mu    intensity
# Output
# X     image matrix

addcircle <- function(c1,c2,r,mu,X){
  
  m1 <- nrow(X); m2 <- ncol(X)
  
  for (i1 in 1:m1){
    for (i2 in 1:m2){
      d <- (i1-c1)^2 + (i2-c2)^2 - r^2
      if(d<0)X[i1,i2] <- mu
    }
  }
  
  return(X)
}

############ End of addcircle ############
