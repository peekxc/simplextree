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

#include <cstddef>		// size_t
#include <vector> 		// vector 
#include <numeric> 		// iota
#include <algorithm>	// transform

struct UnionFind {
	using idx_t = std::size_t; 
	using idx_v = std::vector< size_t >;
  idx_t size; 
  idx_v parent, rank;
  
  UnionFind(const idx_t _size) : size(_size), parent(_size), rank(_size){
		std::iota(parent.begin(), parent.end(), 0);
	}
  
  // Main operations
	void Union(const idx_t x, const idx_t y){
		if (x >= size || y >= size){ return; }
		const idx_t xRoot = Find(x), yRoot = Find(y);
		if (xRoot == yRoot){ return; }
		else if (rank[xRoot] > rank[yRoot]) { parent[yRoot] = xRoot; }
		else if (rank[xRoot] < rank[yRoot]) { parent[xRoot] = yRoot; }
		else if (rank[xRoot] == rank[yRoot]) {
			parent[yRoot] = parent[xRoot];
			rank[xRoot] = rank[xRoot] + 1;
		}
	} 
  const idx_t Find(const idx_t x){
		if (x >= size || parent[x] == x){ return x; }
		else {
			parent[x] = Find(parent[x]);
			return parent[x];
		}
	}
  void AddSets(const idx_t n_sets){
		parent.resize(size + n_sets);
		std::iota(parent.begin() + size, parent.end(), size); // parent initialized incrementally
		rank.resize(size + n_sets, 0); // rank all 0 
		size += n_sets; 
	}

	// Convenience functions
  void UnionAll(const idx_v& idx){
		if (idx.size() <= 1){ return; }
		const idx_t n_pairs = idx.size()-1;
		for (idx_t i = 0; i < n_pairs; ++i){ Union(idx[i], idx[i+1]); }
	}
  idx_v FindAll(const idx_v& idx){
		if (idx.size() == 0){ return idx_v(); }
		idx_v cc = idx_v(idx.size());
		std::transform(idx.begin(), idx.end(), cc.begin(), [this](const size_t i){
			return(Find(i));
		});
		return(cc);
	}
  idx_v ConnectedComponents(){
		idx_v cc = idx_v(size);
		for (size_t i = 0; i < size; ++i){ cc[i] = Find(i); }
		return(cc);
	}
}; // class UnionFind

#endif 
