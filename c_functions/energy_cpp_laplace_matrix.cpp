#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

NumericMatrix energy_non_laplace_homogenous (NumericMatrix& x) {
  int j1=x.nrow();
  int j2=x.ncol();
  int ixr, ixc;
  NumericMatrix sumenergy(j1,j2);
  
  for ( ixr= 0; ixr <j1-1; ixr++){
    for (ixc=0; ixc<j2-1; ixc++){
      sumenergy(ixr,ixc)= sumenergy(ixr, ixc)+abs(x(ixr,ixc)-x(ixr+1,ixc));
      sumenergy(ixr,ixc)=sumenergy(ixr, ixc)+ abs(x(ixr,ixc)-x(ixr,ixc+1));
      sumenergy(j1-1-ixr,j2-1-ixc)= sumenergy(j1-1-ixr,j2-1-ixc)+abs(x(j1-1-ixr,j2-1-ixc)-x(j1-1-ixr-1,j2-1-ixc));
      sumenergy(j1-1-ixr,j2-1-ixc)=sumenergy(j1-1-ixr,j2-1-ixc)+ abs(x(j1-1-ixr,j2-1-ixc)-x(j1-1-ixr,j2-1-ixc-1));
    }
  }
  
  
  for (ixc = 0; ixc< j2-1; ixc++){
    sumenergy(j1-1, ixc) = sumenergy(j1-1, ixc) + abs(x(j1-1,ixc)-x(j1-1,ixc+1));
    sumenergy(0, j2-1-ixc) = sumenergy(0, j2-1-ixc) + abs(x(0,j2-1-ixc)-x(0,j2-1-ixc-1));
  }
  
  for (ixr=0; ixr<j1-1; ixr++){
    sumenergy(ixr,j2-1) = sumenergy(ixr, j2-1) + abs(x(ixr,j2-1)-x(ixr+1,j2-1));
    sumenergy(j1-1-ixr,0) = sumenergy(j1-1-ixr,0) + abs(x(j1-1-ixr,0)-x(j1-1-ixr-1,0));
  }
  
  
  
  
  
  
  return sumenergy;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
