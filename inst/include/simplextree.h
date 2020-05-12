// simplextree.h
// Author: Matt Piekenbrock 
// This package provides a simple implementation of the Simplex tree data structure using Rcpp + STL
// The simplex tree was originally introduced in the following paper: 
// Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
// MIT licensed

#ifndef SIMPLEXTREE_H_
#define SIMPLEXTREE_H_

#include <unordered_map>
#include <queue>
#include <vector>
#include <memory>
#include <tuple>
#include <array>
#include <set>
#include <map>
#include <string>
#include "UnionFind.h"
#include "utility/delegate.hpp"
#include "utility/short_alloc.h"
#include "simplextree/st_traits.hpp"

// Type for simplex labels
typedef std::size_t idx_t;

// Maximum expected size of simplices
static constexpr size_t array_threshold = 9; 

// Buffer type
#include "utility/short_alloc.h"
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

// Simplex tree data structure.
// The Simplex Tree is a normal trie structure, with additional restrictions on the storage order, 
// and a auxiliary map that is used to map 'cousin' simplexes at varying depths.
struct SimplexTree {
  //static constexpr node_ptr empty_face = node_ptr(new node(-1, nullptr));
	
	struct node; 
	using node_ptr = node*;
	using node_uptr = u_ptr< node >; 
	using node_set_t = set< node_uptr, less_ptr< node_uptr > >;
  using simplex_t = vector< idx_t >; 
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using value_type = node_ptr;
  
  // Fields 
  node_uptr root; 															// empty face; initialized to id = 0, parent = nullptr
  map< size_t, vector< node_ptr > > level_map; 	// adjacency map between cousins
  vector< size_t > n_simplexes; 								// tracks the number of simplices if each order
  size_t tree_max_depth; 												// maximum tree depth; largest path from any given leaf to the root. The depth of the root is 0.
  size_t max_id;										 						// maximum vertex id used so far. Only needed by the id generator. 
  size_t id_policy;  														// policy type to generate new ids

  // Filtration-specific fields 
	struct indexed_simplex;   
  vector< bool > included;
  vector< indexed_simplex > fc; 
  

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
		// node(const node& other) : label(other.label), parent(other.parent), children(other.children) { }
		bool operator== (const node& rhs) const noexcept { // equivalence by pointer address 
			return (this == &rhs); 
		} 
		bool operator< (const node& rhs) const { return (label < rhs.label); } // order by label
	};

	// template< size_t I = 0 >
	// constexpr idx_t node_label(node_ptr cn){
	// 	if constexpr (I == 0){
	// 		return(cn->label);
	// 	} else {
	// 		return(node_label< I - 1 >(cn->parent));
	// 	}
	// 	//return I == 0 ? cn->label : node_label< I - 1 >(cn->parent);
	// }

	// template < size_t I, typename Func>
	// void node_labels(const node_ptr cn, Func&& f) {
	// 	auto dispatcher = make_index_dispatcher< I >();
	// 	dispatcher([&](auto idx) { f(node_label< (I - idx - 1) >(cn)); });
	// }

	// TODO: fix noexcept for vector case
	// template < size_t I >
	// v_t< I > full_simplex(const node_ptr cn) noexcept {
	// 	if constexpr(I <= 16){
	// 		v_t< I > result;
	// 		auto dispatcher = make_index_dispatcher< I >();
	// 		dispatcher([&result, &cn](auto idx) { result[idx] = node_label< (I - idx - 1) >(cn); });
	// 		return result;
	// 	} else {
	// 		v_t< I > result;
	// 		if (cn == nullptr){ return; }
	// 		else {
	// 			v_t< I - 1 > r_result = full_simplex3< I - 1 >(cn->parent);
	// 			result.push_back(cn->label);
	// 			result.insert(result.end(), r_result.begin(), r_result.end());
	// 		}
	// 	}
	// }
	template< size_t I = 0 >
	constexpr static idx_t node_label_r(node_ptr cn){
		if constexpr (I == 0){
			return(cn->label);
		} else {
			return(node_label_r< I - 1 >(cn->parent));
		}
	}
	
	// Applies function 'f' to 
	template < size_t I, typename Func>
	static void apply_node_labels(const node_ptr cn, Func&& f) {
		auto dispatcher = make_index_dispatcher< I >();
		dispatcher([&](auto idx) { f(node_label_r< (I - idx - 1) >(cn)); });
	}

	template < typename T , typename U = std::decay_t<T> > 
	static constexpr idx_t node_label(const T& cn) { 
		if constexpr (is_ptr_v< U >){ return cn->label; } 
		else { return cn.label; }
	}
	template < typename T, typename U = std::decay_t<T> > 
	static constexpr node_set_t& node_children(const T& cn) { 
		if constexpr (is_ptr_v< U >){ return cn->children; } 
		else { return cn.children; }
	}
	
	// Encoding to convert a label + depth into a key
	static constexpr auto encode_node(const idx_t label, const idx_t depth) noexcept {
		return szudzik_pair< idx_t, size_t >(label, depth);
	}

	auto node_cousins(node_ptr cn, const idx_t depth) const {
		auto ni = level_map.find(encode_node(node_label(cn), depth)); 
		return ni != level_map.end() ? (*ni).second : vector< node_ptr >();
	}

	// Generates a node id equality predicate 
	// inline std::function<bool(const node_ptr)> SimplexTree::eq_node_id(const idx_t label) const{ 
	// 	return [label](const node_ptr cn)->bool{ return(cn->label == label); };
	// }
	struct eq_node_id {
		const idx_t label; 
		eq_node_id(idx_t new_label) : label(new_label) { };
		
		template < typename T, typename U = std::decay_t<T> >
		bool operator()(const T& cn) const noexcept { 
			if constexpr (is_ptr_v< U >){ return cn->label == label; } 
			else { return cn.label == label; }
			// return(cn->label == label); 
		}
	};
	// constexpr auto eq_node_id(const idx_t label) const noexcept {
	// 	return [label](const node_ptr cn)->bool{ return(cn->label == label); };
	// } 

	
  // Constructor + Destructor 
  SimplexTree() : root(new node(-1, nullptr)), tree_max_depth(0), max_id(0), id_policy(0) { };
  // ~SimplexTree(){
  //   std::vector<idx_t>().swap(n_simplexes);
  // };
  
  // Copy constructor
  SimplexTree(const SimplexTree& st) : root(new node(-1, nullptr)), tree_max_depth(0), max_id(0), id_policy(0) {
    deserialize(st.serialize());
    id_policy = st.id_policy;
    included = st.included;
    fc.reserve(st.fc.size()); 
    std::copy(begin(st.fc), end(st.fc), back_inserter(fc));
  };
  
  // Assignment operator
  SimplexTree& operator=(const SimplexTree& st) {
    deserialize(st.serialize());
    id_policy = st.id_policy;
    included = st.included;
    fc.reserve(st.fc.size()); 
    std::copy(begin(st.fc), end(st.fc), back_inserter(fc));
    return *this;
  };
  void record_new_simplexes(const idx_t k, const idx_t n);// record keeping
  
  constexpr idx_t dimension() const {
    return tree_max_depth == 0 ? 0 : tree_max_depth - 1; 
  }
  
  // Generates a new set of vertex ids, according to the given rule.
  vector< size_t > generate_ids(size_t);
  
  size_t vertex_index(const idx_t) const;
  size_t degree(idx_t) const;
  vector< size_t > degree(vector< idx_t >) const;
  simplex_t adjacent_vertices(const idx_t) const;
  
  // Simplex utilities
  void insert_simplex(simplex_t);
  void remove_simplex(simplex_t);
  bool find_simplex(simplex_t) const;
  
  // Extensions for multiple 
  void insert_simplices(vector< vector< idx_t > >);
  void remove_simplices(vector< vector< idx_t > >);
  vector< bool > find_simplices(vector< simplex_t >) const;
  
  // Recursive helper functions
  node_ptr insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth);
  // node_ptr remove(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  void remove_leaf(node_ptr, idx_t);
  void add_children(node_ptr, const simplex_t&, idx_t);
  void insert(idx_t*, const size_t, const size_t, node_ptr, const idx_t);
  node_ptr find_by_id(const node_set_t&, idx_t) const;
  node_ptr find_node(simplex_t) const;
  node_ptr top_node(node_ptr) const;
  node_ptr find_vertex(const idx_t) const;
  
  bool static is_face(simplex_t, simplex_t);
  
  template < typename Lambda >
  void traverse_facets(node_ptr, Lambda) const;
    
  // Node id equality predicate
  std::function<bool(const node_ptr)> eq_node_id(const idx_t) const;
  
  // Utility 
  simplex_t get_labels(const node_set_t& level, idx_t offset = 0) const;
  size_t depth(node_ptr cn) const;
  size_t max_depth(node_ptr cn) const;
  void remove_subtree(node_ptr parent);
  vector< idx_t > connected_components() const;
  void reindex(vector< idx_t >);
  
  bool is_tree() const; // tests if is fully connected and is a tree
  // bool is_cycle(vector< idx_t > v); // Tests if vertex sequence has a cycle.
    
  // Modifying the complex w/ higher order operations
  bool collapse(node_ptr, node_ptr);
  bool collapseR(simplex_t, simplex_t);
  bool vertex_collapseR(idx_t, idx_t, idx_t);
  bool vertex_collapse(node_ptr, node_ptr, node_ptr);
  void contract(simplex_t);
  void expansion(const idx_t k);
  void expand(node_set_t&, const idx_t);
  void get_cousins() const;

  // Constructs the full simplex from a given node, recursively
  template < size_t I = 0 >
	simplex_t full_simplex_r(node_ptr) const noexcept;
  simplex_t full_simplex(node_ptr, const idx_t depth = 0) const noexcept;
  
	void remove_subtree2(simplex_t);

	
	vector< idx_t > get_vertices() const;
  
  // Serialization/unserialization + saving/loading the complex
  vector< simplex_t > serialize() const;
  void deserialize(vector< simplex_t >);  // void deserializeR(List simplices)
  void save(std::string) const;
  void load(std::string);
  
  // Printing 
  void print_tree() const;
  void print_cousins() const;
  void print_level(node_ptr, idx_t) const;
  void print_subtree(node_ptr) const;
  void print_simplex(node_ptr, bool newline = true) const;

	void clear();

  // Policy for generating ids
  std::string get_id_policy() const;
  void set_id_policy(std::string);
  
  // Indexed simplex 
  struct indexed_simplex {
    size_t parent_idx;  // index of its parent simplex in tree
    idx_t label;        // last(sigma) 
    double index;       // diameter/weight of the simplex 
  };
  
  // Filtration
  vector< size_t > simplex_idx(const size_t) const;
  simplex_t expand_simplex(const vector< size_t >) const;
  void rips(vector< double >, const size_t);
  void threshold_function(double);
  void threshold_index(size_t);
  vector< vector< idx_t > > rips_simplices() const;
  vector< double > rips_weights() const;
  size_t rips_index() const;
  double rips_epsilon() const;
};

#include "simplextree/st_iterators.hpp"
#include "simplextree/st.hpp"


#endif


