// simplextree.h
// Author: Matt Piekenbrock 
// This package provides a simple implementation of the Simplex tree data structure using Rcpp + STL
// The simplex tree was originally introduced in the following paper: 
// Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.
// MIT licensed

#ifndef SIMPLEXTREE_H_
#define SIMPLEXTREE_H_

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// // Progress bar
// #include <RProgress.h>

#include <unordered_map>
#include <queue>
#include <vector>
#include <memory>
#include "UnionFind.h"
#include "utilities.h"

typedef std::size_t idx_t;
using std::set;

// struct node; // forward declaration

template <typename T>
struct ptr_comp {
  bool operator() (const s_ptr<T>& lhs, const s_ptr<T>& rhs) const
  { return (*lhs) < (*rhs); }
};

// Node structure stored by the simplex tree. Contains the following fields:
//  label := integer index type representing the id of simplex it represents
//  parent := (shared) node pointer to its parent in the trie
//  children := connected simplexes whose labels > the current simplex's label
struct node {
  using node_ptr = s_ptr< node >; 
  idx_t label;
  node_ptr parent;
  set< node_ptr, ptr_comp< node > > children;
  node(idx_t id, node_ptr c_parent) : label(id), parent(c_parent){ }
  node(const node& other) : label(other.label), parent(other.parent), children(other.children) { }
  // bool operator== (const node& rhs) const { return (this == &rhs); } // equivalence by pointer address
  bool operator< (const node& rhs) const { return (label < rhs.label); } // order by label
};
using node_ptr = s_ptr< node >; 

// Forward declarations 
struct dfs_iter;
struct bfs_iter;

// Simplex tree data structure.
// The Simplex Tree is a normal trie structure, with additional restrictions on the storage order, 
// and a auxiliary map that is used to map 'cousin' simplexes at varying depths.
struct SimplexTree {
  // Aliases
  using node_ptr = s_ptr< node >; 
  using node_set_t = set< node_ptr, ptr_comp< node > >;
  using simplex_t = vector< idx_t >; 
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using value_type = node_ptr;
  
  // Fields 
  node_ptr root; // empty face; initialized to id = 0, parent = nullptr
  std::unordered_map< std::string, vector< node_ptr > > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers
  vector< size_t > n_simplexes; // tracks the number of simplices if each order
  size_t tree_max_depth; // maximum tree depth; largest path from any given leaf to the root. The depth of the root is 0.
  size_t max_id; // maximum vertex id used so far. Only needed by the id generator. 
  size_t id_policy;  // policy type to generate new ids
    
  // Attach friend class iterators  
  friend struct dfs_iter;
  friend struct bfs_iter;
    
  // Constructor + Destructor 
  SimplexTree(){
    root = node_ptr(new node(-1, nullptr));
    level_map = std::unordered_map< std::string, vector< node_ptr > >();
    n_simplexes = vector<idx_t>(); 
    tree_max_depth = 0; 
    max_id = 0; 
    id_policy = 0; 
  };
  ~SimplexTree(){
    std::vector<idx_t>().swap(n_simplexes);
  };
  // Read-only property 
  SEXP dimension(){ return(tree_max_depth == 0 ? R_NilValue : wrap(tree_max_depth-1)); }
  
  // Utilities
  void clear();
  SEXP as_XPtr();
  void record_new_simplexes(const idx_t k, const idx_t n);// record keeping
  
  
  // Generates a new set of vertex ids, according to the given rule.
  vector< size_t > generate_ids(size_t);
  
  // Read-only "properties" that may be useful for the user
  IntegerMatrix get_k_simplices(const size_t);
  vector< idx_t > get_vertices();
  IntegerMatrix get_edges();
  IntegerMatrix get_triangles();
  IntegerMatrix get_quads();
  
  size_t vertex_index(const idx_t);
  vector< size_t > degree(vector< idx_t >);
  size_t degree(idx_t);
  simplex_t adjacent_vertices(const idx_t);
  
  // Simplex utilities
  void insert_simplex(simplex_t);
  void remove_simplex(simplex_t);
  bool find_simplex(simplex_t);
  
  // Extensions for multiple 
  void insert_simplices(vector< vector< idx_t > >);
  void remove_simplices(vector< vector< idx_t > >);
  vector< bool > find_simplices(vector< simplex_t >);
  
  // R-wrappers
  void insert(SEXP sigma);
  void remove(SEXP sigma);
  LogicalVector find(SEXP sigma);

  // Export utilities
  IntegerMatrix as_adjacency_matrix(); // Exports the 1-skeleton as an adjacency matrix 
  IntegerMatrix as_edge_list(); // Exports the 1-skeleton as an edgelist 
  List as_adjacency_list(); // Exports the 1-skeleton as an adjacency matrix 
  List as_list(); // Exports entire simplex to a list
  
  // Recursive helper functions
  node_ptr insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth);
  node_ptr remove(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  void remove_leaf(node_ptr, idx_t);
  void add_children(node_ptr, const simplex_t&, idx_t);
  void insert(idx_t*, const size_t, const size_t, node_ptr, const idx_t);
  node_ptr find_by_id(const node_set_t&, idx_t);
  node_ptr find_node(simplex_t);
  bool is_face(simplex_t, simplex_t);
  node_ptr top_node(node_ptr);

  // Find a vertex 
  node_ptr find_vertex(const idx_t);
  
  // Node id equality predicate
  std::function<bool(const node_ptr)> eq_node_id(const idx_t);
  
  // Generic way of traversing complex
  template <typename Lambda> 
  void trav_switch(node_ptr, Lambda, std::string, List);
  List traverse_int(SEXP, Function, std::string, Rcpp::Nullable<List>, bool);
  
  // Overloads for regular traversals 
  void traverse(Function, std::string);
  void traverse(SEXP, Function, std::string);
  void traverse(SEXP, Function, std::string, Rcpp::Nullable<List>);
  
  // Overloads for list traversals 
  List ltraverse(Function f, std::string);
  List ltraverse(SEXP simp, Function f, std::string);
  List ltraverse(SEXP, Function, std::string, Rcpp::Nullable<List>);
  
  // Overloads for vector traversals 
  SEXP straverse(Function f, std::string);
  SEXP straverse(SEXP simp, Function f, std::string);
  SEXP straverse(SEXP, Function, std::string, Rcpp::Nullable<List>);

  // Traversal specializations
  template <typename Lambda> void traverse_dfs(node_ptr, Lambda f);
  template <typename Lambda, typename P1, typename P2> void traverse_dfs_if(node_ptr, Lambda, P1, P2);
  template <typename Lambda> void traverse_bfs(node_ptr, Lambda f);
  template <typename Lambda> void traverse_cofaces(node_ptr, Lambda f);
  template <typename Lambda> void traverse_link(node_ptr, Lambda f);
  template <typename Lambda> void traverse_skeleton(node_ptr, Lambda f, size_t k);
  template <typename Lambda> void traverse_max_skeleton(node_ptr, Lambda f, size_t k);
  
  // Utility 
  simplex_t get_labels(const node_set_t& level, idx_t offset = 0);
  size_t depth(node_ptr cn);
  size_t max_depth(node_ptr cn);
  void remove_subtree(node_ptr parent);
  vector< idx_t > connected_components();
  void reindex(SEXP);
  void reindex(vector< idx_t >);
  
  bool is_tree(); // tests if is fully connected and is a tree
  // bool is_cycle(vector< idx_t > v); // Tests if vertex sequence has a cycle.
    
  bool collapse(node_ptr, node_ptr);
  bool collapseR(simplex_t, simplex_t);
  bool vertex_collapseR(idx_t, idx_t, idx_t);
  bool vertex_collapse(node_ptr, node_ptr, node_ptr);
  void contract(simplex_t);
  void expansion(const idx_t k);
  void expand(set< node_ptr, ptr_comp< node > >&, const idx_t);
  
  void get_cousins();

  // Constructs the full simplex from a given node, recursively
  void full_simplex_r(node_ptr, simplex_t&);
  simplex_t full_simplex(node_ptr);
  void remove_subtree2(simplex_t);
  
  // Locate cofaces
  vector< node_ptr > locate_cofaces(node_ptr cn);
  vector< node_ptr > expand_subtree(node_ptr sigma); // unpacks a subtree
  vector< node_ptr > expand_subtrees(vector< node_ptr >); // unpacks multiple subtrees
  
  // Locate links  
  vector< node_ptr > link(node_ptr);
    
  // Depth-first search iterator utility
  dfs_iter begin_dfs(node_ptr), end_dfs(); 
  bfs_iter begin_bfs(node_ptr), end_bfs(); 
  
  // Serialization/unserialization + saving/loading the complex
  vector< simplex_t > serialize();
  void deserialize(vector< simplex_t >);  // void deserializeR(List simplices)
  void save(std::string);
  void load(std::string);
  
  // Printing 
  void print_tree();
  void print_level(node_ptr, idx_t);
  void print_subtree(node_ptr);
  void print_simplex(node_ptr);
  
  // Policy for generating ids
  std::string get_id_policy();
  void set_id_policy(std::string);
  
};

#include <stack>
struct dfs_iter {
  using node_ptr = s_ptr< node >; 
  typedef dfs_iter iterator;
  typedef std::ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef node_ptr value_type;
  typedef node_ptr * pointer;
  typedef node_ptr & reference;
  
  // Members
  node_ptr current; 
  std::stack< node_ptr > node_stack; 
  // TODO: Track depth as well as current simplex in expanded form
  // vector< size_t > depth_stack; 
  // vector< idx_t > c_simplex; 
  // size_t c_depth; 
  
  // Constructor allows DFS searching at any arbitrary subtree
  dfs_iter(node_ptr trie_node) : current(trie_node){
    node_stack = std::stack< node_ptr >();
    // TODO: Track depth as well as current simplex in expanded form
    // c_depth = 0;
    // depth_stack = vector< size_t >();
  }  
  bool operator==(const dfs_iter& other){ return(this->current == other.current); } const
  bool operator!=(const dfs_iter& other){ return(!(this->current == other.current)); } const 
  node_ptr& operator*(){ return(current); } const 
  dfs_iter& operator++(){ // prefix
    if (current != nullptr && !current->children.empty()){
      for (auto nv = current->children.rbegin(); nv != current->children.rend(); ++nv){ 
        node_stack.emplace(*nv); 
      }
      // TODO: Track depth as well as current simplex in expanded form
      // c_depth++;
      // if (c_depth >= depth_stack.size()){
      //   depth_stack.resize(c_depth+1);
      //   c_simplex.resize(c_depth+1);
      // }
      // depth_stack.at(c_depth) += current->children.size()-1;
    }
    if (node_stack.empty()){
      current = nullptr;
    } else { 
      current = node_stack.top();
      node_stack.pop(); 
      // TODO: Track depth as well as current simplex in expanded form
      // c_simplex.at(c_depth) = current->label;
      // if (depth_stack.at(c_depth) <= 0){ --c_depth; } 
      // else { depth_stack.at(c_depth) -= 1; }
    }
    return *this; 
  } 
  // dfs_iter operator++ (int) { // postfix operator
  //   dfs_iter clone(*this);
  //   ++offset;
  //   return clone;
  // }

};

struct bfs_iter {
  using node_ptr = s_ptr< node >; 
  typedef bfs_iter iterator;
  typedef std::ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef node_ptr value_type;
  typedef node_ptr * pointer;
  typedef node_ptr & reference;
  // Members
  node_ptr current; 
  std::queue< node_ptr > node_queue; 
  
  // Constructor allows DFS searching at any arbitrary subtree
  bfs_iter(node_ptr trie_node) : current(trie_node){
    node_queue = std::queue< node_ptr >();
  }  
  bool operator==(const bfs_iter& other){ return(this->current == other.current); } const
  bool operator!=(const bfs_iter& other){ return(!(this->current == other.current)); } const 
  node_ptr& operator*(){ return(current); } const 
  bfs_iter& operator++(){ // prefix
    if (current != nullptr && !current->children.empty()){
      for (auto nv = current->children.begin(); nv != current->children.end(); ++nv){ 
        node_queue.emplace(*nv); 
      }
    }
    if (node_queue.empty()){
      current = nullptr;
    } else { 
      current = node_queue.front();
      node_queue.pop(); 
    }
    return *this; 
  } 
};

#include "simplextree/simplextree.hpp" // implementation
#endif
