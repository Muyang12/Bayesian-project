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
  x_int<-matrix(10,m1,m2) #Defining the initial value of X
  var.x.sample<-matrix(10,m1,m2)#Defining the variance of random sampling x
  tau.int<-matrix(10,m1,m2)#Defining the initial value of tau
  var.tau.sample<-matrix(0.5,m1,m2)#Defining the variance of random sampling tau
  
  ##################### Accept ratio................
  ####The number of accept times and the accept ratio for each unknown parameters
  
  nprop.x <-  matrix(0, m1, m2)#initial accepted times
  naccept.x<-matrix(0,m1,m2)# accept ratio for X
  
  #.................................
  nprop.tau.x<- matrix(0,m1,m2)# initial accepted times
  naccept.tau.x<-matrix(0,m1,m2)# accept ratio for tau
  
  raccept.x.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  raccept.tau.x.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  ##########################Create the space to store the sampling results
  #........create the array,vector to store the update the values...............
  xstore<-array(0,dim =c(m1,m2,(nburn+iteration)))
  mustore<-array(0,dim=c(n1,n2,(nburn+iteration)))
  tau.store<-array(0,dim =c(m1,m2,(nburn+iteration)))
  
  
  x.post.dif<-array(0,dim=c(m1,m2,(nburn+iteration)))
  differ1<-matrix(0,m1,m2)
  differ2<-matrix(0,m1,m2)
  #............
  tau.post.dif<-array(0,dim=c(m1,m2,(nburn+iteration)))
  mu<-matrix(K%*%as.vector(x_int),m1,m2)
  
  Rcpp::sourceCpp(file="R_functions/energy_cpp_laplace_non_homogenous.cpp")
  ###################################
                       lambda_int<-0.001
                       lambda.var=10
                       lambda.post.dif<-vector()
                       lambda.store<-vector()
                       raccept.lambda<-0
                       naccept.lambda<-0
                       nprop.lambda<-0
                       
                       for (k in 1:(nburn+iteration)){
                         
                         for (ithin in 1:nthin){
                           
                           for(i in 1:m1){
                             for (j in 1:m2) {
                               x.new = x_int
                               x.new[i,j] = x.new[i,j] + var.x.sample[i,j]*rnorm(1)
                               nprop.x[i,j] = nprop.x[i,j] + 1
                               if(x.new[i,j]>=0) {
                                 jj<-((j-1)*m1)+i
                                 logxprior<--sum(log(2*tau.int))-sum(energy_non_homogenous(x_int)/(tau.int))
                                 loglikelihood<-sum (Y*log(mu))-sum(mu)##remove the factorial part log(Y!)
                                 logxprior.new<--sum(log(2*tau.int))-sum(energy_non_homogenous(x.new)/(tau.int))
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
                               ##then do sigma.x
                               tau.new=tau.int
                               tau.new[i,j]= tau.new[i,j] + var.tau.sample[i,j]*rnorm(1)
                               nprop.tau.x[i,j] = nprop.tau.x[i,j] + 1
                               if( tau.new[i,j] >=0) {
                                 loghyperprior<--sum(log(tau.int))
                                 loghyperprior.new<--sum(log(tau.new))
                                 logxprior.new<--sum(log(2*tau.new))-sum(energy_non_homogenous(x_int)/(tau.new))
                                 differ2[i,j]<-(logxprior.new+loghyperprior.new)-(loghyperprior+logxprior)
                                 if(differ2[i,j]>log(runif(1))){
                                   naccept.tau.x[i,j] = naccept.tau.x[i,j] + 1
                                   tau.int[i,j]=tau.new[i,j]
                                   loghyperprior=loghyperprior.new
                                   logxprior=logxprior.new
                                 }
                               }
                             }
                           }
                               
                               ########then do lambda_int
                               lambda_new<-lambda_int
                               lambda_new<-lambda_int+lambda.var*rnorm(1)
                               nprop.lambda=nprop.lambda+1
                               if (lambda_new>=0){
                                 loghyperprior<-m1*m2*log(lambda_int)-sum(lambda_int*(tau.int))
                                 loghyperprior.new<-m1*m2*log(lambda_new)-sum(lambda_new*(tau.int))
                                 loglambdaprior<--log(lambda_int)
                                 loglambdaprior.new<--log(lambda_new)
                                 differ3<-(loghyperprior.new+loglambdaprior.new)-(loghyperprior+loglambdaprior)
                                 if (differ3>log(runif(1))){
                                   naccpet.lambda=naccept.lambda+1
                                   lambda_int=lambda_new
                                   loglambdaprior=loglambdaprior.new
                                   loghyperprior=loghyperprior.new
                                 }
                               }
                             
                           
                         }##end thining, run every 10 iteration then only get the 10th, 20th, 30th... iterations
                         
                         if((k <= nburn )&(k == 10*round(k/10))){
                           raccept.x = naccept.x/nprop.x
                           raccept.x.store[,,k]=raccept.x
                           nprop.x = naccept.x = matrix(0, m1, m2)
                           var.x.sample = var.x.sample*(1 + raccept.x/0.234)/2 
                           tmp = array(0,dim = c(2,m1,m2))
                           tmp[1,,]=matrix(10,m1,m2)
                           tmp[2,,]=var.x.sample
                           var.x.sample = apply(tmp,c(2,3),max)
                           raccept.tau.x = naccept.tau.x/nprop.tau.x
                           raccept.tau.x.store[,,k]=raccept.tau.x
                           naccept.tau.x=naccept.tau.x=matrix(0,m1,m2)
                           var.tau.sample = var.tau.sample*(1 + raccept.tau.x/0.234)/2
                           tmp1<-array(0,dim=c(2,m1,m2))
                           tmp1[1,,]=matrix(1,m1,m2)
                           tmp1[2,,]=var.tau.sample
                           var.tau.sample = apply(tmp1,c(2,3),max)
                                      raccept.lambda<-naccept.lambda/nprop.lambda
                                      naccept.lambda=nprop.lambda=0
                                      lambda.var1=lambda.var*(1+raccept.lambda/0.234)/2
                                      lambda.var=max(0.01,lambda.var1)
                         }
                         
                         ##store results 
                         xstore[,,k]=x_int
                         lambda.store[k]=lambda_int
                         tau.store[,,k]=tau.int
                         # cof.store[,,k]=cof.int
                         
                         
                         #Ratio 
                         x.post.dif[,,k]<-differ1
                         tau.post.dif[,,k]<-differ2
                         lambda.post.dif[k]<-differ3
                         
                         
                         #mu
                         mustore[,,k]=mu
                         
                         
                       }
                       #end
                       
                       tau.mean<-apply(tau.store[,,nburn:nburn+iteration],c(1,2),mean)
                       x.mean<-apply(xstore[,,nburn:nburn+iteration],c(1,2),mean)
                       lambda.mean<-mean(lambda.store[nburn:nburn+iteration])
                     return(list(raccept.lambda=raccept.lambda,lambda.store=lambda.store,lambda.mean=lambda.mean,naccept.lambda=naccept.lambda,naccept.tau.x=naccept.tau.x,data_control=data_control,Y=Y,x.mean=x.mean,xstore=xstore,x.post.dif=x.post.dif,tau.post.dif=tau.post.dif,tau.store=tau.store,mustore=mustore,raccept.x=raccept.x,raccept.tau.x.store=raccept.tau.x.store,naccept.x=naccept.x))

}


