// Simple but limited implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <map>
#include <unordered_map>
#include <queue>
#include <vector>
#include <memory>
#include "utility_rcpp.h"

typedef std::size_t idx_t;
using std::vector;
using std::size_t;

// Node structure stored by the simplex tree. Contains the following fields:
//  label := integer index type representing the id of simplex it represents
//  parent := (shared) node pointer to its parent in the trie
//  children := connected simplexes whose labels > the current simplex's label
struct node {
  using node_ptr = s_ptr< node >; 
  idx_t label;
  s_ptr< node > parent;
  std::map< idx_t, node_ptr > children;
  node(idx_t id, node_ptr c_parent) : label(id), parent(c_parent){ 
    children = std::map< idx_t, node_ptr >(); 
  }
  node(const node& other) : label(other.label), parent(other.parent), children(other.children) { }
};
using node_ptr = s_ptr< node >; 

// Forward declarations 
struct dfs_iter;
struct bfs_iter;

// Simplex tree data structure.
// The Simplex Tree is a normal trie structure, with additional restrictions on the storage order, 
// and a auxiliary map that is used to map 'cousin' simplexes at varying depths.
struct SimplexTree {
  using node_ptr = s_ptr< node >; 
  node_ptr root; // empty face; initialized to id = 0, parent = nullptr
  std::unordered_map< std::string, std::vector< node_ptr > > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers
  vector<idx_t> n_simplexes; // tracks the number of simplices if each order
  size_t tree_max_depth; // maximum tree depth; largest path from any given leaf to the root. The depth of the root is 0.
    
  // Attach friend class iterators  
  friend struct dfs_iter;
  friend struct bfs_iter;
    
  // Constructor + Destructor 
  SimplexTree();
  ~SimplexTree();
  
  // Utilities
  SEXP as_XPtr();
  void record_new_simplexes(const idx_t k, const idx_t n);// record keeping
  
  // User-facing API 
  size_t vertex_index(const idx_t v_id);
  // vector< idx_t > vertex_available(idx_t n_vertices);
  vector< idx_t > adjacent_vertices(const int v);
  void add_vertices(const idx_t v_i);
  void remove_vertices(vector< idx_t > vertex_ids);
  // void remove_vertex_cofaces(const int v);
  void insert_simplex(std::vector<idx_t> labels);
  bool find_simplex(vector< idx_t > simplex);
  // void remove_simplex(const IntegerVector& simplex);
  void expansion(const idx_t k);
  IntegerVector get_vertices();

  // Export utilities
  IntegerMatrix as_adjacency_matrix(); // Exports the 1-skeleton as an adjacency matrix 
  List as_adjacency_list(); // Exports the 1-skeleton as an adjacency matrix 
  IntegerMatrix as_edge_list(); // Exports the 1-skeleton as an edgelist 
  
  // Recursive helper functions
  node_ptr insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth);
  node_ptr remove(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  void remove_leaf(node_ptr c_parent, idx_t child_label);
  void add_children(node_ptr c_parent, const std::vector<idx_t>& new_children, idx_t depth);
  void insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  node_ptr find (idx_t label);
  node_ptr find (vector<idx_t> simplex);
  bool is_face(vector<idx_t> tau, vector<idx_t> sigma);
  node_ptr top_node(node_ptr cn);
  void remove_simplex(vector< idx_t > simplex);
  
  // Generic way of traversing complex
  void apply(SEXP, Function, std::string);
  // void as_list_helper(List& res, node_ptr c_node, idx_t depth, idx_t* simplex, const size_t n_keys);
  
  // Utility 
  // idx_t get_dfs_height(node_ptr cnode, idx_t c_height);
  void print_level(node_ptr cnode, idx_t level);
  std::vector<idx_t> get_labels(const std::map<idx_t, node_ptr>& level, const idx_t offset = 0);
  // idx_t intersection_size(std::vector<idx_t> v1, std::vector<idx_t> v2);
  // std::vector<idx_t> intersection(std::vector<idx_t> v1, std::vector<idx_t> v2);
  void expand(std::map<idx_t, node_ptr>& v, const idx_t k, idx_t depth, idx_t* simplex);
  size_t depth(node_ptr cn);
  size_t max_depth(node_ptr cn);
  void remove_subtree(node_ptr parent);
  void collapse(node_ptr tau, node_ptr sigma);
  void collapseR(vector< idx_t >, vector< idx_t >);
  void contract(vector< idx_t > edge);
    
  // Constructs the full simplex from a given node, recursively
  void full_simplex_r(node_ptr node, vector<idx_t>& res);
  vector<idx_t> full_simplex(node_ptr node);
  void remove_subtree2(vector< idx_t > simplex);
  
  // Locate cofaces
  vector<node_ptr> locate_cofaces(node_ptr cn);
  vector<node_ptr> subtree_to_vec(node_ptr sigma); // unpacks a subtree
  
  // Locate links  
  vector< node_ptr > link(node_ptr);
    
  // Depth-first search iterator utility
  dfs_iter begin_dfs(node_ptr), end_dfs(); 
  bfs_iter begin_bfs(node_ptr), end_bfs(); 
  
  // Serialization/unserialization
  void serialize(std::string);
  void unserialize(std::string);
  
  List as_list();
  
  void print_tree();
  void print_subtree(node_ptr cn);
  void print_cousins();
  void print_simplex(node_ptr cn);
  void print_cofaces(vector<idx_t> simplex);
  void test(); // debugging 
    
  // Standrad typedefs for iterating   
  typedef std::ptrdiff_t difference_type;
  typedef std::size_t size_type;
  typedef node_ptr value_type;
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
        node_stack.emplace(nv->second); 
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
        node_queue.emplace(nv->second); 
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

