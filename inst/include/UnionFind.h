//----------------------------------------------------------------------
//                        Disjoint-set data structure 
// File:                        union_find.h
//----------------------------------------------------------------------
// Copyright (c) 2018 Matt Piekenbrock. All Rights Reserved.
//
// Class definition based off of data-structure described here:  
// https://en.wikipedia.org/wiki/Disjoint-set_data_structure

#ifndef UNIONFIND_H_
#define UNIONFIND_H_

#include <Rcpp.h>
using namespace Rcpp;

using std::size_t;
using std::vector;

struct UnionFind {
  size_t size; 
  vector< size_t > parent, rank;
  
  // Ctor and Dtor + Rcpp XPtr handle
  UnionFind(const size_t _size);
  ~UnionFind();
  SEXP as_XPtr();
  
  // Universal operations
  const size_t Find(const size_t x); 
  void AddSets(const size_t n_sets);
  void Union(const size_t x, const size_t y); 
  
  // Overloads for internal convenience functions
  void UnionAll(const vector< size_t >& idx);
  vector< size_t > FindAll(const vector< size_t >& idx);
  
  // Auxilliary 
  IntegerVector getCC();
  void printCC();
}; // class UnionFind

#include "UnionFind/UnionFind.hpp" // implementation

#endif 
