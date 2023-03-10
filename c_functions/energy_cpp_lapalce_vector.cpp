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

NumericVector energy_diverse (NumericMatrix& x) {
  int j1=x.nrow();
  int j2=x.ncol();
  int ixr, ixc;
  NumericVector sumenergy((j1-1)*(j2)+j1*(j2-1));
  
  for ( ixr= 0; ixr <j1-1; ixr++){
    for (ixc=0;ixc<j2;ixc++){
      sumenergy(ixc*(j1-1)+ixr)= abs(x(ixr,ixc)-x(ixr+1,ixc));
    }
  }
  
 
 for ( ixc= 0; ixc <j2-1; ixc++){
   for (ixr=0;ixr<j1;ixr++){
     sumenergy((j1-1)*j2+(j2-1)*ixr+ixc)= abs(x(ixr,ixc)-x(ixr,ixc+1));
   }
 }
 

  
  return sumenergy;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


