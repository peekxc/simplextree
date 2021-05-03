#include <Rcpp.h>
using namespace Rcpp;

#include "simplextree.h"

template < typename label_t = size_t >
struct csd {
  
  size_t nv;
  
  struct node {
    label_t eps;              // filtration value
    label_t mul;              // multiplicity
    vector< node* > neighbors; // star tree neighbors
  };
   
  vector< vector< node > > A; 
  
  csd(size_t n) : nv(n){
    
  }
  
  template < typename Iter >
  void insert(Iter b){
    
  }
};
  
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
