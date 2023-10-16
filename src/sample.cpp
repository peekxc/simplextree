#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]  
#include "simplextree.h"

#include <chrono>
#include <random>

// [[Rcpp::export]]
void expand_f_bernoulli(SEXP stx, const size_t k, const double p){
  SimplexTree& st = *(Rcpp::XPtr< SimplexTree >(stx));  
  
  // Random number generator
  std::random_device random_device;
  std::mt19937 random_engine(random_device());
  std::uniform_real_distribution< double > bernoulli(0.0, 1.0);
  
  // Perform Bernoulli trials for given k, with success probability p
  st.expansion_f(k, [&](node_ptr parent, idx_t depth, idx_t label){
    double q = bernoulli(random_engine);
    if (p == 1.0 | q < p){ // if successful trial
      std::array< idx_t, 1 > int_label = { label };
      st.insert_it(begin(int_label), end(int_label), parent, depth);
    }
  });
}
