#include <Rcpp.h>
using namespace Rcpp;

#include "simplextree.h"

// [[Rcpp::export]]
void cech_complex(SEXP stx, const size_t k, Function seb_f) {
  SimplexTree& st = *(Rcpp::XPtr< SimplexTree >(stx));  
  
  // Then perform the conditional k-expansion
  st.expansion_f(k, [&](node_ptr parent, idx_t depth, idx_t label){
    
    // Collect simplex to test
    auto k_simplex = st.full_simplex(parent, depth);
    k_simplex.push_back(label);
    
    // Insert new face if criteria met 
    bool miniball_valid = as< bool >(seb_f(k_simplex));
    if (miniball_valid){
      std::array< idx_t, 1 > int_label = { label };
      st.insert_it(begin(int_label), end(int_label), parent, depth);
    }
  });
}

