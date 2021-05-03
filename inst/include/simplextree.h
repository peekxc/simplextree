// simplextree.h
// Author: Matt Piekenbrock 
// License: MIT 
// This package provides a simple implementation of the Simplex tree data structure using Rcpp + STL
// The simplex tree was originally introduced in the following paper: 
// Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#ifndef SIMPLEXTREE_H_
#define SIMPLEXTREE_H_

#include <unordered_map>
#include <queue>
#include <vector>
#include <deque>
#include <memory>
#include <tuple>
#include <array>
#include <set>
#include <map>
#include <string>
#include "UnionFind.h"
#include "utility/discrete.h"
#include "utility/combinations.h"
#include "utility/delegate.hpp"
#include "utility/set_utilities.h"

// #include <iostream>

// Type for simplex labels
typedef std::size_t idx_t;

// Maximum expected size of simplices
static constexpr size_t array_threshold = 9; 

// Buffer type
using splex_t = std::vector< idx_t, short_alloc< idx_t, 16, alignof(idx_t)> >;
using splex_alloc_t = typename splex_t::allocator_type::arena_type;

// Aliases
using std::array;
using std::tuple; 
using std::vector;
using std::map;
using std::size_t;
using std::begin; 
using std::end;
using std::set;
using std::find; 
using std::get;

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args){
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// Simplex tree data structure.
// The Simplex Tree is a normal trie structure, with additional restrictions on the storage order, 
// and a auxiliary map that is used to map 'cousin' simplexes at varying depths.
struct SimplexTree {
	struct node; 
	using node_ptr = node*;
	using node_uptr = std::unique_ptr< node >;
	
	struct less_np_label {
    bool operator()(const node_ptr& lhs, const node_uptr& rhs) {
      return lhs->label < rhs->label;
    }
    bool operator()(const node_uptr& lhs, const node_ptr& rhs) {
      return lhs->label < rhs->label;
    }
  };
	
	struct less_ptr {
    bool operator() (const node_uptr& lhs, const node_uptr& rhs) const { return (*lhs) < (*rhs); }
  };
	using node_set_t = set< node_uptr, less_ptr >;
  using simplex_t = vector< idx_t >; 
  using cousin_map_t = std::map< idx_t, vector< node_ptr > >;
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using value_type = node_ptr;
  
  // Fields 
  node_uptr root; 															// empty face; initialized to id = 0, parent = nullptr
  vector< cousin_map_t > level_map;	            // adjacency map between cousins
  std::array< size_t, 32 > n_simplexes = { { 0 } }; 
  // vector< size_t > n_simplexes; 								// tracks the number of simplices if each order
  size_t tree_max_depth; 												// maximum tree depth; largest path from any given leaf to the root. The depth of the root is 0.
  size_t max_id;										 						// maximum vertex id used so far. Only needed by the id generator. 
  size_t id_policy;  														// policy type to generate new ids

	// Node structure stored by the simplex tree. Contains the following fields:
	//  label := integer index type representing the id of simplex it represents
	//  parent := (shared) node pointer to its parent in the trie
	//  children := connected simplexes whose labels > the current simplex's label
	struct node {
		idx_t label;
		node* parent;
		node_set_t children;
		node(idx_t id, node_ptr c_parent) : label(id), parent(c_parent){ }
		node(const node& other) = delete; 
		bool operator== (const node& rhs) const noexcept { // equivalence by pointer address 
			return (this == &rhs); 
		} 
		bool operator< (const node& rhs) const { return (label < rhs.label); } // order by label
		
    // 		struct iterator {
    // 	  	using difference_type = std::ptrdiff_t;
    // 			using value_type = node_ptr; 
    // 			using pointer = node_ptr*; 
    // 			using reference = node_ptr&;
    // 			using iterator_category = std::forward_iterator_tag;
    //   		std::reference_wrapper< node_ptr > cn; 
    //   		iterator(node_ptr current) : cn(std::ref(current)) { };
    // 			bool operator==(const iterator& t) const { return cn == t.cn; }
    // 			bool operator!=(const iterator& t) const { return cn != t.cn; }
    // 			auto operator*() { return cn.get(); };
    // 			auto operator++() { 
    // 			  cn = std::ref(cn.get()->parent);
    // 			  return(*this);
    // 			};
    // 		};
    // 		auto begin() { return iterator(this); };
    //   	auto end() { return iterator(nullptr); };
	};
	
	SimplexTree(const SimplexTree&);
	SimplexTree& operator=(const SimplexTree&);
	
	// Adds node ptr to cousin map
	void add_cousin(node_ptr cn,  const idx_t depth){
	  if (depth_index(depth) >= level_map.size()){
	    level_map.resize(depth_index(depth) + 1);
	  }
	  auto& label_map = level_map[depth_index(depth)][cn->label];
	  auto it = std::find(begin(label_map), end(label_map), cn);
	  if (it == end(label_map)){ label_map.push_back(cn); } // insert 
	}
	
	// Removes node ptr from cousin map
	void remove_cousin(node_ptr cn, const idx_t depth){
	  if (depth_index(depth) >= level_map.size()){ return; }
	  auto& depth_map = level_map[depth_index(depth)];
	  auto cousin_it = depth_map.find(cn->label);
	  if (cousin_it != end(depth_map)){
	    auto& v = cousin_it->second; 
	    v.erase(std::remove(v.begin(), v.end(), cn), v.end());
	  }
	}
	
	template < typename Iter > 
	auto append_node(Iter pos, node_ptr cn, idx_t label, size_t depth) -> node_set_t::iterator;
	
	// Checks if cousins exist 
	bool cousins_exist(const idx_t label, const idx_t depth) const noexcept {
	  if (depth_index(depth) >= level_map.size()){ return false; }
	  return level_map[depth_index(depth)].find(label) != end(level_map[depth_index(depth)]);
	}
	
	const cousin_map_t::mapped_type& cousins(const idx_t label, const idx_t depth) const {
	  return level_map[depth_index(depth)].at(label);
	}
	
	template < typename Lambda >
	void traverse_cousins(const idx_t label, const idx_t depth, Lambda f) const {
	  if (depth_index(depth) >= level_map.size()){ return; }
	  if (cousins_exist(label, depth)){
	    const auto& c_cousins = level_map[depth_index(depth)].at(label);
	    std::for_each(begin(c_cousins), end(c_cousins), f);
	  }
	};
	
  // Constructor
  SimplexTree() : root(new node(-1, nullptr)), tree_max_depth(0), max_id(0), id_policy(0) { };
  
  // Generates a new set of vertex ids, according to the given rule.
  auto generate_ids(size_t) -> vector< size_t >;
  auto degree(idx_t) const -> size_t;
  auto adjacent_vertices(const idx_t) const -> simplex_t;
  auto record_new_simplexes(const idx_t k, const int n) -> void;// record keeping
  auto dimension() const -> idx_t { return tree_max_depth == 0 ? 0 : tree_max_depth - 1; }
  
  template< typename Iterable > 
  auto insert(Iterable v) -> void; 
    
  template< bool use_lex = false, typename Iter >
  auto insert_it(Iter, Iter, node_ptr, const idx_t) -> void;
  
  template< typename Iterable > 
  auto find(Iterable v) const -> node_ptr; 
  
  template< typename Iter >
  auto find_it(Iter, Iter, node_ptr cn) const -> node_ptr; 
  
  auto find_by_id(const node_set_t&, idx_t) const -> node_ptr;
  
  auto remove(node_ptr cn) -> void;
  auto remove_leaf(node_ptr, idx_t) -> void;
  auto remove_subtree(node_ptr parent) -> void;
  
  template < typename Lambda >
  void traverse_facets(node_ptr, Lambda) const;
  
  // Utility 
  auto is_face(simplex_t, simplex_t) const -> bool;
  auto depth(node_ptr cn) const -> size_t;
  auto max_depth(node_ptr cn) const -> size_t;
  auto connected_components() const -> vector< idx_t >;
  auto reindex(vector< idx_t >) -> void;
  auto is_tree() const -> bool;
  auto get_vertices() const -> vector< idx_t >;
  auto clear() -> void;
  
  // Modifying the complex w/ higher order operations
  auto collapse(node_ptr, node_ptr) -> bool;
  auto vertex_collapse(idx_t, idx_t, idx_t) -> bool;
  auto vertex_collapse(node_ptr, node_ptr, node_ptr) -> bool;
  auto contract(simplex_t) -> void;

  auto expansion(const idx_t k) -> void;
  
  template < typename Lambda >
  auto expansion_f(const idx_t, Lambda&&) -> void;
  
  template < typename Lambda >
  auto expand_f(node_set_t&, const idx_t, size_t, Lambda&&) -> void;
  
  template < typename Lambda > 
	void traverse_up(node_ptr, const size_t, Lambda&&) const noexcept;
	
	template< typename OutputIt > 
	void full_simplex_out(node_ptr, const idx_t, OutputIt) const noexcept;
	auto full_simplex(node_ptr cn, const idx_t depth = 0) const noexcept -> simplex_t;
	
	template < typename T > // Assumes T is pointer type
	static constexpr auto node_label(T& cn) -> idx_t { return cn->label; }
	template < typename T > // Assumes T is pointer type
	static constexpr auto node_children(T& cn) -> node_set_t& { return cn->children; }
	static constexpr auto depth_index(const idx_t depth) noexcept -> idx_t { return(depth - 2); }
  
  // Printing 
  template < typename OutputStream > void print_tree(OutputStream&) const;
  template < typename OutputStream > void print_cousins(OutputStream&) const;
  template < typename OutputStream > void print_level(OutputStream&, node_ptr, idx_t) const;
  template < typename OutputStream > void print_subtree(OutputStream&, node_ptr) const;
  template < typename OutputStream > void print_simplex(OutputStream&, node_ptr, bool newline = true) const;


  // Policy for generating ids
  std::string get_id_policy() const;
  void set_id_policy(std::string);
  
};

#include "simplextree/st_iterators.hpp"
#include "simplextree/st_filtration.hpp"
#include "simplextree/st.hpp"


#endif


