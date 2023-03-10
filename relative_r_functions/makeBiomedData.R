########## R script: makeBiomedData ###########

# To produce a simuluated biomedical dataset.
# The output also contains the true vector x.

# Last changed: 06 FEB 2020

makeBiomedData <- function(control_list){

if(exists("n1",       where=control_list)){      n1 = control_list[["n1"]]}      else{n1=40}
if(exists("n2",       where=control_list)){      n2 = control_list[["n2"]]}      else{n2=40}
if(exists("K",        where=control_list)){       K = control_list[["K"]]}       else{cat("Error in call to makeArchaeoData: K missing.\n")}
if(exists("noise",    where=control_list)){   noise = control_list[["noise"]]}   else{noise="gaussian"}
if(exists("sigmaeps", where=control_list)){sigmaeps = control_list[["sigmaeps"]]}else{sigmaeps=1}
if(exists("sc",       where=control_list)){      sc = control_list[["sc"]]}      else{sc=1}

  # Set dimensions of X:
  m1 <- n1*sc; m2 <- n2*sc
  
  # Add background intensity, or eta=T uses pre-set value:
  eta <- 4
  
  # Set theta value:
  theta <- scan(text = "1099 14 40 7 20 14 3 8 20 1.6 22 26 1.1 0.57 2.3")
  # theta[14] is now obsolete, was delta
  
  #addcircle <- function(c1,c2,r,mu,X){
    
    #m1 <- nrow(X); m2 <- ncol(X)
    
    #for (i1 in 1:m1){
      #for (i2 in 1:m2){
        #d <- (i1-c1)^2 + (i2-c2)^2 - r^2
       # if(d<0)X[i1,i2] <- mu
      #}
    #}
   # return(X)
  #}
  
  
  if(eta == T){eta <- theta[15]} else{theta[15] <- eta}
  
  if( (min(theta[-length(theta)]) > 0) & #each elements are positive
      (max(theta[c( 2, 5, 8,11)]) < m1) & # set the center location x-axis
      (max(theta[c( 3, 6, 9,12)]) < m2) & # set the center location y-axis
      (max(theta[c( 4, 7,10,13)]) < min(m1,m2)/2 )) # radius
  {
    hot.Ci <- theta[ 1]
    eta <- theta[15]
  } else {stop("Error in makeBiomedData: Invalid theta values.")}
  
  X <- matrix(0, m1, m2)
  X <- addcircle(theta[2]*sc,theta[3]*sc,theta[4]*sc,hot.Ci,X)
  X <- addcircle(theta[5]*sc,theta[6]*sc,theta[7]*sc,hot.Ci,X)
  X <- addcircle(theta[8]*sc,theta[9]*sc,theta[10]*sc,hot.Ci,X)
  X <- addcircle(theta[11]*sc,theta[12]*sc,theta[13]*sc,hot.Ci,X)
    
  # Calcule forward projection:
  MU <- matrix(data_control$K %*% as.vector(X),29,58) + 5
  
  # Generate data:
  Y <- matrix(0,29, 58)
  
  if(noise == "poisson"){
    Y <- matrix(rpois(n1*n2,MU),n1,n2)
  } else if (noise == "gaussian"){
    Y <- MU + sigmaeps*matrix(rnorm(n1*n2),n1,n2)
  } else if (noise == "none"){
    Y <- MU
  }  else {stop("Error in makeBiomedData: Noise option not recognized.")}
  
  return(list(X=X,MU=MU,Y=Y))
  
}
  
############ End of makeBiomedData ############
