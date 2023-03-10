switch<-function(data_control){
  if(exists("X",        where=data_control)){ X = data_control[["X"]]}       else{X=F}
  if(exists("Y",        where=data_control)){ Y = data_control[["Y"]]}       else{cat("Error in call to mcmc: Y missing.\n")}
  if(exists("n1",       where=data_control)){n1 = data_control[["n1"]]}      else{n1=40}
  if(exists("n2",       where=data_control)){n2 = data_control[["n2"]]}      else{n2=40}
  if(exists("K",        where=data_control)){ K = data_control[["K"]]}       else{cat("Error in call to mcmc: K missing.\n")}
  if(exists("nburn",    where=data_control)){nburn = data_control[["nburn"]]}       else{nburn=1}
  if(exists("iteration",    where=data_control)){iteration = data_control[["iteration"]]}       else{iteration=1}
  if(exists("nthin",    where=data_control)){nthin = data_control[["nthin"]]}       else{nthin=10}
  if(exists("ntrunc", where=data_control)){ntrunc =data_control[["ntrunc"]]}  else{ntrunc=F}
  if(exists("alpha", where=data_control)){alpha=data_control[["alpha"]]} else{alpha=F}
  if(exists("beta", where=data_control)){beta=data_control[["beta"]]} else{beta=F}
  if(exists("v_l", where=data_control)){v_l=data_control[["v_l"]]} else{v_l=100}
  if(exists("v_g", where=data_control)){v_g=data_control[["v_g"]]} else{v_g=10}
  #............
  
  sc=1 ; m1 <- n1*sc ; m2 <- n2*sc;
  # ntrunc.thresh = 10
  #v_l=1000
  #v_g=100
  nprior<-2*4+4*(m1-2)*(m2-2)+6*(m1-2)+6*(m2-2)#Defining the number of prior
  x_int<-matrix(10,m1,m2) #Defining the initial value of X
  theta_matrix_int<-matrix(0,m1,m2)
  var.x.sample<-matrix(10,m1,m2)#Defining the variance of random sampling x
  var.p.sample<-0.5 #Defining the variance of random sampling p
  zeta.l<-matrix(v_l,m1,m2)#Defining the initial value of tau  
  zeta.g<-matrix(v_g,m1,m2)
  p_int<-0.5
  #Define the variance of sigma.x
  ##################### Accept ratio................
  ####The number of accept times and the accept ratio for each unknown parameters
  raccept.x<-matrix(0,m1,m2)
  nprop.x <-  matrix(0,m1,m2)#initial accepted times
  naccept.x<-matrix(0,m1,m2)# accept ratio for X
  raccept.x.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  #.................................
  nprop.theta<-matrix(0,m1,m2)
  naccept.theta<-matrix(0,m1,m2)
  raccept.theta.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  raccept.p=0
  naccept.p=0
  nprop.p<-0
  raccept.p.store<-as.vector(rep(0,(nburn+iteration)))
  

  ##########################Create the space to store the sampling results
  #........create the array,as.vector to store the update the values...............
  xstore<-array(0,dim =c(m1,m2,(nburn+iteration)))
  mustore<-array(0,dim=c(n1,n2,(nburn+iteration)))
  x.post.dif<-array(0,dim=c(m1,m2,(nburn+iteration)))
  
  theta_store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  theta.post.dif<-array(0,dim =c(m1,m2,(nburn+iteration)))
  p.store<-as.vector(rep(0,(nburn+iteration)))
  p.post.dif<-as.vector(rep(0,(nburn+iteration)))
  differ1<-matrix(0,m1,m2)
  differ2<-matrix(0,m1,m2)
  differ3<-0
  #............
  
  mu<-matrix(K%*%as.vector(x_int),m1,m2)
  #data_control[["energy_type"]]="laplace"
  Rcpp::sourceCpp(file="R_functions/energy_cpp_gaussian_non_homogenous.cpp")
  Rcpp::sourceCpp(file="R_functions/energy_cpp_laplace_non_homogenous.cpp")
 logxprior_matrix<-list(t_l=(-log(2*(zeta.l))-(energy_non_homogenous_laplace(x_int)/zeta.l)),t_g=(-log(sqrt(2*pi)*zeta.g)-(energy_non_homogenous_gaussian(x_int)/(2*zeta.g^2))))
  logxprior<-sum(theta_matrix_int*logxprior_matrix$t_l+(1-theta_matrix_int)*logxprior_matrix$t_g)
  loglikelihood<-sum (Y*log(mu))-sum(mu)##remove the factorial part log(Y!)
  ###################################
  
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
          nprop.theta[i,j]=nprop.theta[i,j]+1
          loghyper.theta<-sum(theta_matrix_int)*log(p_int)+sum(1-theta_matrix_int)*log(1-p_int)
          loghyper.theta.new<-sum(theta_matrix.new)*log(p_int)+sum(1-theta_matrix.new)*log(1-p_int)
          logxprior<-sum(theta_matrix_int*logxprior_matrix$t_l+(1-theta_matrix_int)*logxprior_matrix$t_g)
          logxprior.new<-sum(theta_matrix.new*logxprior_matrix$t_l+(1-theta_matrix.new)*logxprior_matrix$t_g)
          differ2[i,j]=(logxprior.new+loghyper.theta.new)-(logxprior+loghyper.theta)
          if (differ2[i,j]>log(runif(1))){
            naccept.theta[i,j]=naccept.theta[i,j]+1
            theta_matrix_int=theta_matrix.new
            logxprior=logxprior.new
            loghyper.theta=loghyper.theta.new
          }
      
      ##################### Updata P
      p.new<-p_int
      p.new<-p_int+var.p.sample*rnorm(1)
      nprop.p=nprop.p+1
      if (p.new<1 & p.new>0) {
        loghyper.theta<-sum(theta_matrix_int)*log(p_int)+sum(1-theta_matrix_int)*log(1-p_int)
        loghyper.theta.new<-sum(theta_matrix_int)*log(p.new)+sum(1-theta_matrix_int)*log(1-p.new)
        loghyper.p<--log(p_int)
        loghyper.p.new<--log(p.new)
        differ3<-(loghyper.theta.new+loghyper.p.new)-(loghyper.theta+loghyper.p)
        if (differ3 >log(runif(1))){
          naccept.p=naccept.p+1
          p_int<-p.new
          loghyper.p=loghyper.p.new
          loghyper.theta=loghyper.theta.new
        }
      }
    }
    }
    }
    
    ##end thining, run every 10 iteration then only get the 10th, 20th, 30th... iterations
    
    if((k <= nburn )&(k == 10*round(k/10))){
      raccept.x = naccept.x/nprop.x
      raccept.x.store[,,k]=raccept.x
      nprop.x = naccept.x = matrix(0, m1, m2)
      var.x.sample = var.x.sample*(1 + raccept.x/0.234)/2
      tmp = array(0,dim = c(2,m1,m2))
      tmp[1,,]=matrix(0.3,m1,m2)
      tmp[2,,]=var.x.sample
      var.x.sample = apply(tmp,c(2,3),max)
      raccept.p=naccept.p/nprop.p
      raccept.p.store[k]=raccept.p
      naccept.p=nprop.p=0
      var.p.sample=var.p.sample*(1+raccept.p/0.234)/2
      var.p.sample=max(0.1,var.p.sample)
    }
    
    
    ##store results
    xstore[,,k]=x_int
    theta_store[,,k]=theta_matrix_int
    p.store[k]=p_int
    
    
    #Ratio
    x.post.dif[,,k]<-differ1
    theta.post.dif[,,k]<-differ2
    p.post.dif[k]<-differ3
    
    #mu
    mustore[,,k]=mu
    raccept.x.store[,,k]=raccept.x
    raccept.p.store[k]=raccept.p
  } 
  
  
  #end the iteration
  x.mean<-apply(xstore[,,nburn:nburn+iteration],c(1,2),mean)
  theta.mean<-apply(theta_store[,,nburn:nburn+iteration],c(1,2),mean)
  
  return(list(theta.mean=theta.mean,x.mean=x.mean,raccept.x=raccept.x,raccept.p=raccept.p,alpha=alpha,beta=beta,theta_store=theta_store,p.store=p.store,naccept.x=naccept.x,data_control=data_control,Y=Y,x.mean=x.mean,xstore=xstore,x.post.dif=x.post.dif,mustore=mustore,raccept.p.store=raccept.p.store))
  }
