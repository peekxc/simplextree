#include <Rcpp.h>
using namespace Rcpp;

#include "simplextree.h"

// [[Rcpp::export]]
IntegerVector bench_preorder(SEXP stree) {
  XPtr< SimplexTree > st_ptr(stree);
  SimplexTree& st = *st_ptr;
  auto l_tr = st::preorder< true >(&st, st.root.get());
  vector< size_t > adds; 
  adds.reserve(std::accumulate(st.n_simplexes.begin(), st.n_simplexes.end(), 0));
  st::traverse(l_tr, [&adds](node_ptr cn, idx_t depth, simplex_t sigma){
    adds.push_back(std::accumulate(sigma.begin(), sigma.end(), 0));
    return true; 
  });
  return(wrap(adds));
}

// [[Rcpp::export]]
IntegerVector bench_levelorder(SEXP stree) {
  XPtr< SimplexTree > st_ptr(stree);
  SimplexTree& st = *st_ptr;
  auto l_tr = st::level_order< true >(&st, st.root.get());
  vector< size_t > adds; 
  adds.reserve(std::accumulate(st.n_simplexes.begin(), st.n_simplexes.end(), 0));
  st::traverse(l_tr, [&adds](node_ptr cn, idx_t depth, simplex_t sigma){
    adds.push_back(std::accumulate(sigma.begin(), sigma.end(), 0));
    return true; 
  });
  return(wrap(adds));
}
