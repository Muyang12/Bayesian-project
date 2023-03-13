
######create the null vector
z<-vector()
########define a function for stroing information
for (i in 1:length(listsig)){
  sigma.x<-listsig[i]
  t<-as.vector(get(paste0("x.mse.g",sigma.x)))
  z<-cbind(z,t)
}

####estimation MSE
colMeans(z)


########If each element is vector
z<-do.call(cbind,(lapply(paste0("x.mse.g",listsig,sep=""),get)))
