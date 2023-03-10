########## R script: RUN_makeBiomedData ##################
setwd("C:/Users/mmmz/OneDrive - University of Leeds/First year/RGA_functions")
setwd("C:/Users/mmmz/OneDrive - University of Leeds/Second year/R code")
source("C:/Users/mmmz/OneDrive - University of Leeds/Second year/R code/getKB.R")
source("makeBiomedData.R")

setwd("C:/Users/mmmz/OneDrive - University of Leeds/fourth year/")
# Set output folder and example name:
  output.folder <- "SIMULATED_DATA_mouse/"
  example <- "gaussian"

# Set how many to replicate (0 = no replicates):
  nrep = 3; noffset = 0
  
# Create control list:
  data_control = list()

# Add variables for Kernel matrix:
  data_control[["n1"]] <- 29
  data_control[["n2"]] <- 58
  data_control[["ntrunc"]] <- 10
  data_control[["type"]] <- "B"
  data_control[["noise"]]<-"gaussian"
  data_control[["sigmaeps"]] <- 0
# Add variables for make_*_data:
  data_control[["eta"]]<-0.5
   data_control[["delta"]]<-1
   data_control[["K"]] <- getKB(data_control)
  #identify matrix
  #data_control[["K"]]<-diag(29*58)
   arch<-RUN_makeBiomeData()
  res <- makeBiomedData(data_control)
  MU <- res$MU ; X  <- res$X ; Y_sim2 <- res$Y
  
  
  
  
  
  Y <- matrix(0,29, 58)
  
  if(noise == "poisson"){
    Y <- matrix(rpois(n1*n2,MU*X_new),n1,n2)
  } else if (noise == "gaussian"){
    Y <- MU + sigmaeps*matrix(rnorm(n1*n2),n1,n2)
  } else if (noise == "none"){
    Y <- MU
  }  else {stop("Error in makeBiomedData: Noise option not recognized.")}
  
  return(list(X=X,MU=MU,Y=Y))
  
  
  MU <- matrix(data_control$K %*% as.vector(X_new),29,58) + 5
  
  Y <- MU + matrix(rnorm(n1*n2),n1,n2)
  n1=29
  n2=58
# Draw some pictures:

  par(mfrow=c(1,3),mar=c(0.5, 1, 1, 1), mgp=c(1.7, 0.65, 0))
  
  image(Y,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
  legend("topright","(a)",bg="white",adj=1,text.width=0)
  
  image(MU, xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
  legend("topright","(b)",bg="white",adj=1,text.width=0)
  
  image(Y_sim2,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
  legend("topright","(c)",bg="white",adj=1,text.width=0)
  

# Save output:
if(nrep==0){
    save(data_control,X,MU,Y,
       file=paste(output.folder,example,
                  "_sd=",data_control[["sigmaeps"]],
                  "_noiseN",data_control[["noise"]],
                  "_data.Rdata",sep=""))
} else {
# Create and save replicate files:
    sigmaeps = data_control[["sigmaeps"]]
    n1 <- data_control[["n1"]]
    n2 <- data_control[["n2"]]
    
    for (irep in noffset+1:nrep){
      Y <- matrix(rpois(n1*n2,MU),n1,n2)
    # Save output:
    save(data_control,X,MU,Y,
         file=paste(output.folder,example,
                    "_eta_",2.3,
                    "_sd=",data_control[["sigmaeps"]],
                    "_noiseN",data_control[["noise"]],
                    "_irep=",irep,
                    "_data.RData",sep=""))      
    }
}
  
MUpoisson<-MU
Xpoisson<-X
Ypoisson<-Y

MUpoissonN<-MU
XpoissonN<-X
YpoissonN<-Y


MUpoisson0<-MU
Xpoisson0<-X
Ypoisson0<-Y

par(mfrow=c(1,3),mar=c(0.5, 1, 1, 1), mgp=c(1.7, 0.65, 0))

image(Xpoisson-XpoissonN,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(a)",bg="white",adj=1,text.width=0)

image(MUpoisson-MUpoissonN, xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(b)",bg="white",adj=1,text.width=0)

image(Ypoisson-YpoissonN,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(c)",bg="white",adj=1,text.width=0)


par(mfrow=c(1,3),mar=c(0.5, 1, 1, 1), mgp=c(1.7, 0.65, 0))

image(Xpoisson-Xpoisson0,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(a)",bg="white",adj=1,text.width=0)

image(MUpoisson-MUpoisson0, xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(b)",bg="white",adj=1,text.width=0)

image(Ypoisson-Ypoisson0,  xlab="",ylab="",col=grey(32:0/32),axes=F);box() 
legend("topright","(c)",bg="white",adj=1,text.width=0)







  # End RUN_makeBiomedData #################################

