//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.cpp
//----------------------------------------------------------------------
// Copyright (c) 2018 Matt Piekenbrock. All Rights Reserved.
//
// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#include "UnionFind.h"

// Export as an Rcpp module
RCPP_MODULE(union_find_module) {
  Rcpp::class_<UnionFind>("UnionFind")
  .constructor< std::size_t >()
  .field_readonly( "size", &UnionFind::size)
  .field_readonly( "parent", &UnionFind::parent)
  .field_readonly( "rank", &UnionFind::rank)
  .method("as_XPtr", &UnionFind::as_XPtr)
  .method("print", &UnionFind::printCC )
  .method("connected_components", &UnionFind::getCC)
  .method("find", &UnionFind::Find)
  .method("find_all", &UnionFind::FindAll)
  .method("union", &UnionFind::Union)
  .method("union_all", &UnionFind::UnionAll)
  .method("add_sets", &UnionFind::AddSets)
  ;
}
