mcmc_hyper<-function(control_list){
  #if(exists("X",        where=data_control)){ X = data_control[["X"]]}       else{X=F}
  #if(exists("Y",        where=data_control)){ Y = data_control[["Y"]]}       else{cat("Error in call to mcmc: Y missing.\n")}
  #if(exists("sc",       where=data_control)){n1 = control_list[["sc"]]}      else{sc=F}
  if(exists("K",        where=control_list)){ K = control_list[["K"]]}       else{cat("Error in call to mcmc: K missing.\n")}
  if(exists("nburn",    where=control_list)){nburn = control_list[["nburn"]]}       else{nburn=100}
  if(exists("iteration",    where=control_list)){iteration = control_list[["iteration"]]}       else{iteration=100}
  if(exists("nthin",    where=control_list)){nthin = control_list[["nthin"]]}       else{nthin=1}
  # if(exists("ntrunc", where=control_list)){ntrunc =control_list[["ntrunc"]]}  else{ntrunc=F}
  if(exists("X",        where=control_list)){ X = control_list[["X"]]}       else{X=F}
  if(exists("Y",        where=control_list)){ Y = control_list[["Y"]]}       else{cat("Error in call to mcmc: Y missing.\n")}
  if(exists("n1",       where=control_list)){n1 = control_list[["n1"]]}      else{n1=40}
  if(exists("n2",       where=control_list)){n2 = control_list[["n2"]]}      else{n2=40}

  sc=1 ; m1 <- n1*sc ;m2 <- n2*sc;
  nprior<-2*4+4*(m1-2)*(m2-2)+6*(m1-2)+6*(m2-2)#Defining the number of prior
  var.zeta.sample<-0.1 #Defining the variance of random sampling tau
  zeta.int<-matrix(10,m1,m2)#Defining the initial value of tau
  x_int<-matrix(10,m1,m2) #Defining the initial value of X
  x_trans_int<-matrix(2,m1,m2) #Defining the initial value of X
  x_int<-exp(x_trans_int)
  var.x_trans.sample<-matrix(10,m1,m2)#Defining the variance of random sampling x
  nprop.x_trans <- matrix(0, m1, m2)#initial accepted times
  raccept.x_trans.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  x_trans.store<-array(0,dim =c(m1,m2,(nburn+iteration)))  
  zeta.int<-matrix(10,m1,m2)#Defining the initial value of tau
  var.zeta.sample<-matrix(0.5,m1,m2)#Defining the variance of random sampling tau
  tau.int=zeta.int
  ##################### Accept ratio................
  ####The number of accept times and the accept ratio for each unknown parameters
  
  nprop.x <-  matrix(0, m1, m2)#initial accepted times
  naccept.x<-matrix(0,m1,m2)# accept ratio for X
  
  #.................................
  nprop.zeta.x<- matrix(0,m1,m2)# initial accepted times
  naccept.zeta.x<-matrix(0,m1,m2)# accept ratio for tau
  
  raccept.x.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  raccept.zeta.x.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  ##########################Create the space to store the sampling results
  #........create the array,vector to store the update the values...............
  xstore<-array(0,dim =c(m1,m2,(nburn+iteration)))
  #cof.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  mustore<-array(0,dim=c(n1,n2,(nburn+iteration)))
  zeta.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  tau.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  x.post.dif<-array(0,dim=c(m1,m2,(nburn+iteration)))
  differ1<-matrix(0,m1,m2)
  differ2<-matrix(0,m1,m2)
  #............
  zeta.post.dif<-array(0,dim=c(m1,m2,(nburn+iteration)))
  mu<-matrix(K%*%as.vector(x_int),m1,m2)
  
  Rcpp::sourceCpp(file="R_functions/energy_cpp_laplace_non_homogenous.cpp")
  
  
  
  for (k in 1:(nburn+iteration)){
    
    for (ithin in 1:nthin){
      
      for(i in 1:m1){
        for (j in 1:m2) {
          x_trans.new = x_trans_int
          x_trans.new[i,j] = x_trans.new[i,j] + var.x_trans.sample[i,j]*rnorm(1)
          x_int<-exp(x_trans_int)
          x.new<-exp(x_trans.new)
          nprop.x[i,j] = nprop.x[i,j] + 1
            jj<-((j-1)*m1)+i
            logxprior<-(-m1)*m2*log(2)-sum(log(zeta.int))-sum(energy_non_homogenous(x_int)/(zeta.int))+sum(x_trans_int)
            loglikelihood<-sum (Y*log(mu))-sum(mu)##remove the factorial part log(Y!)
            logxprior.new<-(-m1)*m2*log(2)-sum(log(zeta.int))-sum(energy_non_homogenous(x.new)/(zeta.int))+sum(x_trans.new)           
            update_value<-matrix((x.new[i,j]-x_int[i,j])*K[,jj],n1,n2)
            mu.new<-mu+update_value
            loglikelihood.new<-sum (Y*log(mu.new))-sum(mu.new)##remove the factorial part log(Y!)
            differ1[i,j]=(logxprior.new+loglikelihood.new)-(logxprior+loglikelihood)
            if(differ1[i,j]> log(runif(1))) {
              naccept.x[i,j]=naccept.x[i,j]+1
              x_int[i,j]=x.new[i,j]
              x_trans_int[i,j]=x_trans.new[i,j]
              mu=mu.new
              logxprior=logxprior.new
              loglikelihood=loglikelihood.new
            }
          }
          ##then do sigma.x
          zeta.new=zeta.int
          zeta.new[i,j]= zeta.new[i,j] + var.zeta.sample[i,j]*rnorm(1)
          nprop.zeta.x[i,j] = nprop.zeta.x[i,j] + 1
          if( zeta.new[i,j] >=0) {
            loghyperprior<-m1*m2*log(lambda)-lambda*sum(zeta.int)
            loghyperprior.new<-m1*m2*log(lambda)-lambda*sum(zeta.new)
            logxprior.new<-(-m1)*m2*log(2)-sum(log(zeta.new))-sum(energy_non_homogenous(x_int)/(zeta.new))+sum(x_trans_int)
            differ2[i,j]<-(logxprior.new+loghyperprior.new)-(loghyperprior+logxprior)
            if(differ2[i,j]>log(runif(1))){
              naccept.zeta.x[i,j] = naccept.zeta.x[i,j] + 1
              zeta.int[i,j]=zeta.new[i,j]
              loghyperprior=loghyperprior.new
              logxprior=logxprior.new
            }
          }
      }
    }##end thining, run every 10 iteration then only get the 10th, 20th, 30th... iterations
    
    if((k <= nburn )&(k == 10*round(k/10))){
      raccept.x=naccept.x/nprop.x
      var.x_trans.sample = var.x_trans.sample*(1 + raccept.x/0.234)/2 
      tmp = array(0,dim = c(2,m1,m2))
      tmp[1,,]=matrix(10,m1,m2)
      tmp[2,,]=var.x_trans.sample
      var.x_trans.sample = apply(tmp,c(2,3),max)
      tmp = array(0,dim = c(2,m1,m2))
     # tmp[1,,]=matrix(10,m1,m2)
     # tmp[2,,]=var.x.sample
     # var.x.sample = apply(tmp,c(2,3),max)
      raccept.zeta.x = naccept.zeta.x/nprop.zeta.x
      raccept.zeta.x.store[,,k]=raccept.zeta.x
      naccept.zeta.x=naccept.zeta.x=matrix(0,m1,m2)
      var.zeta.sample = var.zeta.sample*(1 + raccept.zeta.x/0.234)/2
      tmp1<-array(0,dim=c(2,m1,m2))
      tmp1[1,,]=matrix(1,m1,m2)
      tmp1[2,,]=var.zeta.sample
      var.zeta.sample = apply(tmp1,c(2,3),max)
    }
    
  ##store results 
  xstore[,,k]=x_int
  x.mean<-apply(xstore[,,nburn:nburn+iteration],c(1,2),mean)
  x_trans.store[,,k]=x_trans_int
  #Ratio 
  x.post.dif[,,k]<-differ1
  zeta.post.dif[,,k]<-differ2
  
  #sigma.x
  zeta.store[,,k]<-zeta.int
  
  zeta.mean<-mean(zeta.store[nburn:nburn+iteration])
  #mu
  mustore[,,k]=mu
    
  } #end
  zeta.mean<-apply(zeta.store[,,nburn:nburn+iteration],c(1,2),mean)
  x.mean<-apply(xstore[,,nburn:nburn+iteration],c(1,2),mean)
  return(list(x_trans.store=x_trans.store,data_control=data_control,Y=Y,x.mean=x.mean,xstore=xstore,x.post.dif=x.post.dif,zeta.post.dif=zeta.post.dif,zeta.store=zeta.store,mustore=mustore,raccept.x.store=raccept.x.store,raccept.zeta.x.store=raccept.zeta.x.store))
}


