#The function for creation Transformation matrix
#..........................................
getKB <- function(data_control){
  
  if(exists("n1",     where=data_control)){n1     =data_control[["n1"]]}else{n1=40}       
  if(exists("n2",     where=data_control)){n2     =data_control[["n2"]]}else{n2=40}       
  if(exists("ntrunc", where=data_control)){ntrunc =data_control[["ntrunc"]]}else{ntrunc=F}    
  if(exists("delta",  where=data_control)){delta  =data_control[["delta"]]}else{delta=0.5}       
  if(exists("sc",     where=data_control)){sc     =data_control[["sc"]]}else{sc=1}        
  
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
########################........................................
