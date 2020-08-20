//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.cpp
//----------------------------------------------------------------------
// Copyright (c) 2018 Matt Piekenbrock. All Rights Reserved.
//
// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure
#include <Rcpp.h>
using namespace Rcpp; 

#include "UnionFind.h"


SEXP as_XPtr(UnionFind* uf){
  Rcpp::XPtr< UnionFind > p(uf, false);
  return(p);
}
                   
void printCC(UnionFind* uf){
  for (size_t j = 0; j < uf->size; ++j){ Rcout << uf->Find(j) << " "; }
  Rcout << std::endl;
}

// Export as an Rcpp module
RCPP_MODULE(union_find_module) {
  Rcpp::class_<UnionFind>("UnionFind")
  .constructor< std::size_t >()
  .field_readonly( "size", &UnionFind::size)
  .field_readonly( "parent", &UnionFind::parent)
  .field_readonly( "rank", &UnionFind::rank)
  .method("as_XPtr", &as_XPtr)
  .method("print", &printCC )
  .method("connected_components", &UnionFind::ConnectedComponents)
  .method("find", &UnionFind::Find)
  .method("find_all", &UnionFind::FindAll)
  .method("union", &UnionFind::Union)
  .method("union_all", &UnionFind::UnionAll)
  .method("add_sets", &UnionFind::AddSets)
  ;
}
