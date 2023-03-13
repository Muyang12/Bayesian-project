##########Importing a list of parameters
listsig<-c(0.1,0.2,0.3,0.4,0.5)
###########Importing simulated dataset
load("C:/Users/mmmz/OneDrive - University of Leeds/third year/REAL_DATA/smoothing.RData")
#Setting the iterations and burning-in
nburn=100
iteration=1000
#######Now importing the estimations along with different parameters' setting
for (i in 1:length(listsig)) {
  sigma.x<-listsig[i]
  load(paste0("C:/Users/mmmz/Downloads/dis_ex_fixed_p/2 _ ",sigma.x," _ .RData")) #type name of document
  attach(result)
  x.mse.l<-apply(sweep(xstore,c(1,2),X,"-")^2,c(1,2),mean)
  x.mean.l<-apply(xstore[,,(nburn+1):(nburn+iteration)], c(1,2), mean)
  x.sd.l<-apply(xstore[,,(nburn+1):(nburn+iteration)], c(1,2), sd)
  x.lci.l    = apply(xstore[,,(nburn+1):(nburn+iteration)], c(1,2), quantile, 0.025) #credit interval
  x.uci.l  = apply(xstore[,,(nburn+1):(nburn+iteration)], c(1,2), quantile, 0.975)
  x.ciw.l= x.uci.l - x.lci.l
  # mu.mean.l<-apply(mustore[,,(nburn+1):(nburn+iteration)],c(1,2),mean)
  assign(paste0("result.l",sigma.x),result)
  assign(paste0（"x.mse.","g",sigma.x),x.mse.l)
  assign(paste0("x.sd.","g",sigma.x),x.sd.l)
  assign(paste0("x.mean","g",sigma.x),x.mean.l)
  assign(paste0("x.lci","g",sigma.x),x.lci.l)
  assign(paste0("x.uci","g",sigma.x),x.uci.l)
  assign(paste0("x.ciw","g",sigma.x),x.ciw.l)
  assign(paste0("xstore.","g",sigma.x),xstore)
  # assign(paste0("mustore.",sigma.x),mustore)
  detach(result)
}


