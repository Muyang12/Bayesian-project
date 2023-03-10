distribution_external<-function(data_control){
  if(exists("X",        where=data_control)){ X = data_control[["X"]]}       else{X=F}
  if(exists("Y",        where=data_control)){ Y = data_control[["Y"]]}       else{cat("Error in call to mcmc: Y missing.\n")}
  if(exists("n1",       where=data_control)){n1 = data_control[["n1"]]}      else{n1=40}
  if(exists("n2",       where=data_control)){n2 = data_control[["n2"]]}      else{n2=40}
  if(exists("K",        where=data_control)){ K = data_control[["K"]]}       else{cat("Error in call to mcmc: K missing.\n")}
  if(exists("nburn",    where=data_control)){nburn = data_control[["nburn"]]}       else{nburn=1}
  if(exists("iteration",    where=data_control)){iteration = data_control[["iteration"]]}       else{iteration=10}
  if(exists("nthin",    where=data_control)){nthin = data_control[["nthin"]]}       else{nthin=10}
  if(exists("ntrunc", where=data_control)){ntrunc =data_control[["ntrunc"]]}  else{ntrunc=F}
  if(exists("lambda",where=data_control)){lambda=data_control[["lambda"]]}else{cat("Error in call to mcmc: lambda missing.\n")}
  if(exists("factor1",where=data_control)){factor1=data_control[["factor1"]]}else{factor1=1}
  if(exists("factor2",where=data_control)){factor2=data_control[["factor2"]]}else{factor2=1}
  if(exists("opt",where=data_control)){opt=data_control[["opt"]]}else{opt=0.5}
  if(exists("opt_par",where=data_control)){opt_par=data_control[["opt_par"]]}else{opt_par=0.5}
  # Set dimensions of X (here these are fixed):
  sc=1 ; m1 <- n1*sc ; m2 <- n2*sc; 
  # ntrunc.thresh = 10
  lambda_z<-(lambda-mean(lambda))/sd(lambda)
  lambda_p<-factor1*pnorm(factor2*lambda_z)
  lambda_p[which(lambda_p==1)]=0.99999
  # ntrunc.thresh = 10
  
  nprior<-2*4+4*(m1-2)*(m2-2)+6*(m1-2)+6*(m2-2)#Defining the number of prior
  x_int<-matrix(10,m1,m2) #Defining the initial value of X
  theta_matrix.new<-theta_matrix_int<-matrix(0,m1,m2)
  var.x.sample<-matrix(10,m1,m2)#Defining the variance of random sampling x
  zeta.l<-matrix(100,m1,m2)#Defining the initial value of tau  
  zeta.g<-matrix(10,m1,m2)
  #Define the variance of sigma.x
  
  ##################### Accept ratio................
  ####The number of accept times and the accept ratio for each unknown parameters
  raccept.x<-matrix(0,m1,m2)
  nprop.x <-  matrix(0, m1, m2)#initial accepted times
  naccept.x<-matrix(0,m1,m2)# accept ratio for X
  raccept.x.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  #.................................
  nprop.theta<-matrix(0,m1,m2)
  naccept.theta<-matrix(0,m1,m2)
  raccept.theta.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  K_bay<-matrix(0,m1,m2)
  K_store<-array(0,dim=c(m1,m2,(nburn+iteration)))
  
  ##########################Create the space to store the sampling results
  #........create the array,as.vector to store the update the values...............
  xstore<-array(0,dim =c(m1,m2,(nburn+iteration)))
  mustore<-array(0,dim=c(n1,n2,(nburn+iteration)))
  p_store<-array(0,dim=c(n1,n2,(nburn+iteration)))
  x.post.dif<-array(0,dim=c(m1,m2,(nburn+iteration)))
  
  theta_store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  theta.post.dif<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  differ1<-matrix(0,m1,m2)
  differ2<-matrix(0,m1,m2)
  logxprior_matrix<-list()
  #............
  
  mu<-matrix(K%*%as.vector(x_int),m1,m2)
  #data_control[["energy_type"]]="laplace"
  Rcpp::sourceCpp(file="/nobackup/mmmz/Final _year/final_chapter/distribution/R_functions/energy_cpp_gaussian_non_homogenous.cpp")
  Rcpp::sourceCpp(file= "/nobackup/mmmz/Final _year/final_chapter/distribution/R_functions/energy_cpp_laplace_non_homogenous.cpp")
  theta_matrix_int=matrix(0,29,58)
  
  ###################################
  p_selection=matrix(opt,29,58)
  p_matrix_int=matrix(opt_par,29,58)
  #...............................MCMC
  
  for (k in 1:(nburn+iteration)){
    
    for (ithin in 1:nthin){
      #############Update x, for the f(y/x)f(x/theta)
      for(i in 1:m1){
        for (j in 1:m2) {
          x.new = x_int
          x.new[i,j] = x.new[i,j] + var.x.sample[i,j]*rnorm(1)
          nprop.x[i,j] = nprop.x[i,j] + 1
          if(x.new[i,j]>=0) {
            jj<-((j-1)*m1)+i
            #not right should only one contains 
            logxprior_matrix<-list(t_l=(-log(2*(zeta.l))-(energy_non_homogenous_laplace(x_int)/zeta.l)),t_g=(-log(sqrt(2*pi)*zeta.g)-(energy_non_homogenous_gaussian(x_int)/(2*zeta.g^2))))
            logxprior<-sum(theta_matrix_int*logxprior_matrix$t_l+(1-theta_matrix_int)*logxprior_matrix$t_g)
            loglikelihood<-sum (Y*log(mu))-sum(mu)##remove the factorial part log(Y!)
            logxprior_matrix.new<-list(t_l=(-log(2*(zeta.l))-(energy_non_homogenous_laplace(x.new)/zeta.l)),t_g=(-log(sqrt(2*pi)*zeta.g)-(energy_non_homogenous_gaussian(x.new)/(2*zeta.g^2))))
            logxprior.new<-sum(theta_matrix_int*logxprior_matrix.new$t_l+(1-theta_matrix_int)*logxprior_matrix.new$t_g)
            update_value<-matrix((x.new[i,j]-x_int[i,j])*K[,jj],n1,n2)
            mu.new<-mu+update_value
            loglikelihood.new<-sum (Y*log(mu.new))-sum(mu.new)##remove the factorial part log(Y!)
            differ1[i,j]=(logxprior.new+loglikelihood.new)-(logxprior+loglikelihood)
            if(differ1[i,j]> log(runif(1))) {
              naccept.x[i,j]=naccept.x[i,j]+1
              x_int[i,j]=x.new[i,j]
              mu=mu.new
              logxprior=logxprior.new
              loglikelihood=loglikelihood.new
            }
          }
          ############ Update theta
          
          theta_matrix.new=theta_matrix_int
          theta_matrix.new[i,j]=1-theta_matrix_int[i,j]
          nprop.theta[i,j]=nprop.theta[i,j]+1#
          loghyper.theta<-sum(theta_matrix_int*log(p_selection))+sum((1-theta_matrix_int)*log(1-p_selection))
          loghyper.theta.new<-sum(theta_matrix.new*log(p_selection))+sum((1-theta_matrix.new)*log(1-p_selection))
          logxprior<-sum(theta_matrix_int*logxprior_matrix$t_l+(1-theta_matrix_int)*logxprior_matrix$t_g)
          logxprior.new<-sum(theta_matrix.new*logxprior_matrix$t_l+(1-theta_matrix.new)*logxprior_matrix$t_g)
          differ2[i,j]=(logxprior.new+loghyper.theta.new)-(logxprior+loghyper.theta)
          if (differ2[i,j]>log(runif(1))){
            naccept.theta[i,j]=naccept.theta[i,j]+1
            theta_matrix_int=theta_matrix.new
            logxprior=logxprior.new
            loghyper.theta=loghyper.theta.new
          }
          }
          }
          
          
          
          #Bayesian factor
          p_matrix_new=p_matrix_int
          p_matrix_new=lambda_p
          null_1=(theta_matrix_int*log(p_matrix_int))+((1-theta_matrix_int)*log(1-p_matrix_int))
          null_2=(theta_matrix_int*log(p_matrix_new))+((1-theta_matrix_int)*log(1-p_matrix_new))
          K_bay<-null_1-null_2
          lab_neg<-ifelse(K_bay<0,1,0)
          p_selection=lab_neg*lambda_p+(1-lab_neg)*p_matrix_int       
    }
    
      
    
      
      
      ##end thining, run every 10 iteration then only get the 10th, 20th, 30th... iterations
      
      if((k <= nburn )&(k == 10*round(k/10))){
        raccept.x = naccept.x/nprop.x
        nprop.x = naccept.x = matrix(0, m1, m2)
        var.x.sample = var.x.sample*(1 + raccept.x/0.234)/2
        tmp = array(0,dim = c(2,m1,m2))
        tmp[1,,]=matrix(0.3,m1,m2)
        tmp[2,,]=var.x.sample
        var.x.sample = apply(tmp,c(2,3),max)
      }
      
      
      ##store results
      xstore[,,k]=x_int
      theta_store[,,k]=theta_matrix_int
      K_store[,,k]=K_bay
      p_store[,,k]=p_selection
      
      #Ratio
      x.post.dif[,,k]<-differ1
      theta.post.dif[,,k]<-differ2
      
      #mu
      mustore[,,k]=mu
      raccept.x.store[,,k]=raccept.x
    }  
    
  
  #end the iteration
  x.mean<-apply(xstore[,,nburn:nburn+iteration],c(1,2),mean)
  theta.mean<-apply(theta_store[,,nburn:nburn+iteration],c(1,2),mean)
  p_store.mean<-apply(p_store[,,nburn:nburn+iteration],c(1,2),mean)
  return(list(opt=opt,opt_par=opt_par,p_store=p_store,p_store.mean=p_store.mean,K_store=K_store,lambda_p=lambda_p,theta.mean=theta.mean,x.mean=x.mean,raccept.x=raccept.x,theta_store=theta_store,naccept.x=naccept.x,data_control=data_control,Y=Y,x.mean=x.mean,xstore=xstore,x.post.dif=x.post.dif,mustore=mustore))
}




