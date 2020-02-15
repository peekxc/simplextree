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
#include <tuple>
#include "UnionFind.h"

#include <Rcpp/Benchmark/Timer.h>

typedef std::size_t idx_t;

// Aliases
using std::vector;
using std::map;
using std::size_t;
using std::begin; 
using std::end;
using std::set;
using std::find; 
using std::get;
template <typename T> using s_ptr = std::shared_ptr<T>; // Shared pointer
template <typename T> using u_ptr = std::unique_ptr<T>; // Unique pointer
template <typename T> using w_ptr = std::weak_ptr<T>;   // Weak pointer// struct node; // forward declaration

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
  //static constexpr node_ptr empty_face = node_ptr(new node(-1, nullptr));
  
  // Aliases
  using node_ptr = s_ptr< node >; 
  using node_set_t = set< node_ptr, ptr_comp< node > >;
  using simplex_t = vector< idx_t >; 
  using difference_type = std::ptrdiff_t;
  using size_type = std::size_t;
  using value_type = node_ptr;
  
  // Fields 
  node_ptr root; // empty face; initialized to id = 0, parent = nullptr
  // std::unordered_map< std::string, vector< node_ptr > > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers
  map< size_t, vector< node_ptr > > level_map; 
  vector< size_t > n_simplexes; // tracks the number of simplices if each order
  size_t tree_max_depth; // maximum tree depth; largest path from any given leaf to the root. The depth of the root is 0.
  size_t max_id; // maximum vertex id used so far. Only needed by the id generator. 
  size_t id_policy;  // policy type to generate new ids
  
  // Forward declaration indexed simplex 
  struct indexed_simplex;   
  
  // Filtration-specific fields 
  vector< bool > included;
  vector< indexed_simplex > fc; 
    
  // Attach friend class iterators  
  friend struct dfs_iter;
  friend struct bfs_iter;
    
  // Constructor + Destructor 
  SimplexTree() : root(new node(-1, nullptr)), tree_max_depth(0), max_id(0), id_policy(0) {
    // level_map = std::unordered_map< std::string, vector< node_ptr > >();
    // n_simplexes = vector<idx_t>(); 
  };
  ~SimplexTree(){
    std::vector<idx_t>().swap(n_simplexes);
  };
  
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
  
  // Read-only property 
  SEXP dimension(){ return(tree_max_depth == 0 ? R_NilValue : wrap(tree_max_depth-1)); }
  
  // Utilities
  void clear();
  SEXP as_XPtr();
  void record_new_simplexes(const idx_t k, const idx_t n);// record keeping
  
  // Generates a new set of vertex ids, according to the given rule.
  vector< size_t > generate_ids(size_t);
  
  // Read-only "properties" that may be useful for the user
  IntegerMatrix get_k_simplices(const size_t) const;
  vector< idx_t > get_vertices() const;
  IntegerMatrix get_edges() const;
  IntegerMatrix get_triangles() const;
  IntegerMatrix get_quads() const;
  
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
  
  // R-wrappers
  void insert(SEXP sigma);
  void remove(SEXP sigma);
  LogicalVector find(SEXP sigma) const;

  // Export utilities
  IntegerMatrix as_adjacency_matrix() const; // Exports the 1-skeleton as an adjacency matrix 
  IntegerMatrix as_edge_list() const; // Exports the 1-skeleton as an edgelist 
  List as_adjacency_list() const; // Exports the 1-skeleton as an adjacency matrix 
  List as_list() const; // Exports entire simplex to a list
  
  // Recursive helper functions
  node_ptr insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth);
  node_ptr remove(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  void remove_leaf(node_ptr, idx_t);
  void add_children(node_ptr, const simplex_t&, idx_t);
  void insert(idx_t*, const size_t, const size_t, node_ptr, const idx_t);
  node_ptr find_by_id(const node_set_t&, idx_t) const;
  node_ptr find_node(simplex_t) const;
  bool is_face(simplex_t, simplex_t) const;
  node_ptr top_node(node_ptr) const;
  node_ptr find_vertex(const idx_t) const;
  
  // Node id equality predicate
  std::function<bool(const node_ptr)> eq_node_id(const idx_t) const;
  
  // Generic way of traversing complex
  template <typename Lambda> 
  void trav_switch(node_ptr, Lambda, std::string, List) const;
  
  List traverse_int(SEXP, Function, std::string, Rcpp::Nullable<List>, bool) const;
  
  // Overloads for regular traversals 
  void traverse(Function, std::string) const;
  void traverse(SEXP, Function, std::string) const;
  void traverse(SEXP, Function, std::string, Rcpp::Nullable<List>) const;
  
  // // Overloads for list traversals 
  // List ltraverse(Function f, std::string) const;
  // List ltraverse(SEXP simp, Function f, std::string) const;
  // List ltraverse(SEXP, Function, std::string, Rcpp::Nullable<List>) const;
  // 
  // // Overloads for vector traversals 
  // SEXP straverse(Function f, std::string) const;
  // SEXP straverse(SEXP simp, Function f, std::string) const;
  // SEXP straverse(SEXP, Function, std::string, Rcpp::Nullable<List>) const;

  // Traversal specializations
  template <typename Lambda> void traverse_dfs(node_ptr, Lambda f) const;
  template <typename Lambda, typename P1, typename P2> void traverse_dfs_if(node_ptr, Lambda, P1, P2) const;
  template <typename Lambda> void traverse_bfs(node_ptr, Lambda f) const;
  template <typename Lambda, typename P1, typename P2> void traverse_bfs_if(node_ptr, Lambda, P1, P2) const;
  template <typename Lambda> void traverse_cofaces(node_ptr, Lambda f) const;
  template <typename Lambda> void traverse_link(node_ptr, Lambda f) const;
  template <typename Lambda> void traverse_skeleton(node_ptr, Lambda f, size_t k) const;
  template <typename Lambda> void traverse_max_skeleton(node_ptr, Lambda f, size_t k) const;
  template <typename Lambda> void traverse_facets(node_ptr, Lambda) const;
  
  // Simplex versions
  template < typename Lambda > void traverse_faces_s(node_ptr, Lambda) const;
  
  // Utility 
  simplex_t get_labels(const node_set_t& level, idx_t offset = 0) const;
  size_t depth(node_ptr cn) const;
  size_t max_depth(node_ptr cn) const;
  void remove_subtree(node_ptr parent);
  vector< idx_t > connected_components() const;
  void reindex(SEXP);
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
  void expand(set< node_ptr, ptr_comp< node > >&, const idx_t);
  
  void get_cousins() const;

  // Constructs the full simplex from a given node, recursively
  void full_simplex_r(node_ptr, simplex_t&) const;
  simplex_t full_simplex(node_ptr) const;
  void remove_subtree2(simplex_t);
  
  // Locate cofaces
  vector< node_ptr > locate_cofaces(node_ptr cn) const;
  vector< node_ptr > expand_subtree(node_ptr sigma) const; // unpacks a subtree
  vector< node_ptr > expand_subtrees(vector< node_ptr >) const; // unpacks multiple subtrees
  
  // Locate links  
  vector< node_ptr > link(node_ptr) const;
    
  // Depth-first search iterator utility
  dfs_iter begin_dfs(node_ptr) const, end_dfs() const; 
  bfs_iter begin_bfs(node_ptr) const, end_bfs() const; 
  
  // Serialization/unserialization + saving/loading the complex
  vector< simplex_t > serialize() const;
  void deserialize(vector< simplex_t >);  // void deserializeR(List simplices)
  void save(std::string) const;
  void load(std::string);
  
  // Printing 
  void print_tree() const;
  void print_level(node_ptr, idx_t) const;
  void print_subtree(node_ptr) const;
  void print_simplex(node_ptr) const;
  
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
