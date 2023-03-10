
######################### load observation data 
load("4_cylinders_data.RData")
data_control<-list()
data_control[["n1"]]=nrow(Y)
data_control[["n2"]]=ncol(Y)
data_control[["ntrunc"]]=10
n1<-nrow(Y)
n2<-ncol(Y)
sc<-1  
m1 <- n1*sc
m2 <- n2*sc

#######################################getKB function
getKB <- function(control_list){
  
  if(exists("n1",     where=control_list)){n1     =control_list[["n1"]]}else{n1=40}       
  if(exists("n2",     where=control_list)){n2     =control_list[["n2"]]}else{n2=40}       
  if(exists("ntrunc", where=control_list)){ntrunc =control_list[["ntrunc"]]}else{ntrunc=F}    
  if(exists("delta",  where=control_list)){delta  =control_list[["delta"]]}else{delta=0.5}       
  if(exists("sc",     where=control_list)){sc     =control_list[["sc"]]}else{sc=1}        
  
  # Set dimensions of X:
  m1 <- n1*sc; m2 <- n2*sc
  
  if(is.logical(ntrunc)){
    if(!ntrunc){cat("Warning in getKB: Kernel truncation not being used, the calculation may be slow...\n")
      ntrunc <- ceiling(sqrt(n1^2 + n2^2))} 
    else {ntrunc <- 10}}
  
  K <- matrix(0, n1*n2, m1*m2)
  
  for (i1 in 1:m1){
    for (i2 in 1:m2){
      
      ii <- (i2 - 1)*m1 + i1
      
      j1.min <- max(1, round((i1-1/2)/sc+1/2 - ntrunc, 0))
      j1.max <- min(n1,round((i1-1/2)/sc+1/2 + ntrunc, 0))
      j2.min <- max(1, round((i2-1/2)/sc+1/2 - ntrunc, 0))
      j2.max <- min(n2,round((i2-1/2)/sc+1/2 + ntrunc, 0))
      
      for (j1 in j1.min:j1.max){
        for (j2 in j2.min:j2.max){
          
          jj <- (j2 - 1)*n1 + j1
          
          k1 <- (j1 - 1/2) - (i1 - 1/2)/sc
          k2 <- (j2 - 1/2) - (i2 - 1/2)/sc
          
          K[jj,ii] <- dnorm(k1,0,delta)*dnorm(k2,0,delta)/sc^2
          
        }}     # j1/j2
    }}# i1//i2
  
  return(K)
}



###############################################################add circle function
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


##################################3# write optimization function 


optFunction <- function(theta){
  
  center_x1<-theta[1]
  center_x2<-theta[2]
  center_x3<-theta[3]
  center_x4<-theta[4]
  center_y1<-theta[5]
  center_y2<-theta[6]
  center_y3<-theta[7]
  center_y4<-theta[8]
  radius_1<-theta[9]
  radius_2<-theta[10]
  radius_3<-theta[11]
  radius_4<-theta[12]
  mu<-theta[13]
  eta<-theta[14]
  data_control[["delta"]]=theta[15]
  
  # Set dimensions of X:
  
  X <- matrix(0, m1, m2)
  X <- addcircle(center_x1,center_y1,radius_1,mu,X)
  X <- addcircle(center_x2,center_y2,radius_2,mu,X)
  X<- addcircle(center_x3,center_y3,radius_3,mu,X)
  X <- addcircle(center_x4,center_y4,radius_4,mu,X)
  
  # Calculate forward projection/transformation matrix:
  K<-getKB(data_control)
  MU <- matrix(K %*% as.vector(X),n1,n2)
  
  fit<-mean((Y-MU-eta)^2)
  return(fit)
}
############################
#install.packages("DEoptim")
library("DEoptim")
#global optimum
globalest<-DEoptim(fn=optFunction,lower = c(21,21,8,13,16,26,20,41,1,1,1,1,1000,1,0.5),upper=c(24,24,11,16,19,29,23,44,10,10,10,10,1300,10,1.5))
View(Y)
#local
# see message 
optim(par=c(21,21,8,13,16,26,20,41,1,1,1,1,1000,0.5),fn=optFunction,method = "L-BFGS-B" ,NULL,lower = c(21,21,8,13,16,26,20,41,1,1,1,1,1000,0.5),upper=c(29,29,29,29,58,58,58,58,14,14,14,14,1300,Inf,1.5))

# only return initial value 
optim(par=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1),fn=optFunction,method = "BFGS")

#######################put the estimated parameter into makeBiomeData function 
makeBiomedData <- function(control_list){
  
  if(exists("n1",       where=control_list)){      n1 = control_list[["n1"]]}      else{n1=40}
  if(exists("n2",       where=control_list)){      n2 = control_list[["n2"]]}      else{n2=40}
  if(exists("K",        where=control_list)){       K = control_list[["K"]]}       else{cat("Error in call to makeArchaeoData: K missing.\n")}
  if(exists("noise",    where=control_list)){   noise = control_list[["noise"]]}   else{noise="gaussian"}
  if(exists("sigmaeps", where=control_list)){sigmaeps = control_list[["sigmaeps"]]}else{sigmaeps=1}
  if(exists("sc",       where=control_list)){      sc = control_list[["sc"]]}      else{sc=1}
  
  # Set theta value: (using code to see the optimized value: globalest$optim)
  globalest$optim
  theta <- scan(text = "13 13 21 14 40 4 14 41 6 1 2.5 6.6 1100 27")
  
  if( (min(theta[-length(theta)]) > 0) & 
      (max(theta[c( 1, 2, 3,4)]) < m1) & 
      (max(theta[c( 5, 6, 7,8)]) < m2) & 
      (max(theta[c( 9, 10,11,12)]) < min(m1,m2)/2 ))
  {
    hot.Ci <- theta[13]
    eta <- theta[14]
  } else {stop("Error in makeBiomedData: Invalid theta values.")}
  
  
  X <- matrix(0, m1, m2)
  X <- addcircle(theta[1]*sc,theta[5]*sc,theta[9]*sc,hot.Ci,X)
  X <- addcircle(theta[2]*sc,theta[6]*sc,theta[10]*sc,hot.Ci,X)
  X <- addcircle(theta[3]*sc,theta[7]*sc,theta[11]*sc,hot.Ci,X)
  X <- addcircle(theta[4]*sc,theta[8]*sc,theta[12]*sc,hot.Ci,X)
  
  MU <- matrix(K %*% as.vector(X),n1,n2) + eta
  
  if(noise == "poisson"){
    Y <- matrix(rpois(n1*n2,MU),n1,n2)
  } else if (noise == "gaussian"){
    Y <- MU + sigmaeps*matrix(rnorm(n1*n2),n1,n2)
  } else if (noise == "none"){
    Y <- MU
  }  else {stop("Error in makeBiomedData: Noise option not recognized.")}
  
  return(list(X=X,MU=MU,Y=Y))
  
}

#create simulation data with Poisson noise

data_control[["noise"]]<-"poisson"
data_control[["K"]] <- K
res <- makeBiomedData(data_control)
MU <- res$MU ; X  <- res$X ; Y_sim  <- res$Y

# Draw some pictures:

par(mfrow=c(1,3),mar=c(0.5, 1, 1, 1), mgp=c(1.7, 0.65, 0))

image(X,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(a)",bg="white",adj=1,text.width=0)

image(MU, xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(b)",bg="white",adj=1,text.width=0)

image(Y_sim,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(c)",bg="white",adj=1,text.width=0)

dev.off()
#global optimal
21.138683   22.446608    8.389241   14.171953   16.000423   26.644610   20.401289   41.000044    3.011304    1.032866    1.550189    7.086662 1089.794486    8.406445    0.888268
