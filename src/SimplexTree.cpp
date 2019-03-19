// SimplexTree.cpp
// Simple but limited implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include "SimplexTree.h"

// Constructs the simplex tree
SimplexTree::SimplexTree() {  
  root = node_ptr(new node(-1, nullptr));
  level_map = std::unordered_map< std::string, vector< node_ptr > >();
  n_simplexes = vector<idx_t>(); 
  tree_max_depth = 0; 
}

SimplexTree::~SimplexTree() { 
  std::vector<idx_t>().swap(n_simplexes);
}


SEXP SimplexTree::as_XPtr(){
  Rcpp::XPtr< SimplexTree> p(this, false);
  return(p);
}

// Search the level map (cousins) to quickly get the adjacency relations. 
// The set of adjcency relations are the 0-simplexes connected to a given vertex v. 
vector< idx_t > SimplexTree::adjacent_vertices(const int v){

  // Resulting vector to return
  vector< idx_t > res = vector< idx_t >(); 
  
  // First extract the vertices which labels > v
  std::string key = std::to_string(v) + "-2";
  std::unordered_map< std::string, vector< node_ptr > >::iterator it = level_map.find(key);
  if (it != level_map.end()){
    // Get the cousins parents 
    vector< node_ptr > cousins = (*it).second;
    std::for_each(cousins.begin(), cousins.end(), [&res](const node_ptr& c_node ){
      res.push_back((*c_node).parent->label);
    });
  }
  
  // Then get the vertices with labels < v
  std::map< idx_t, node_ptr > top_vertices = root->children;
  if (top_vertices.find(v) != top_vertices.end()){
    node_ptr vertex = top_vertices.at(v);
    std::map< idx_t, node_ptr > v_children = vertex->children;
    for (auto& kv: v_children){
      res.push_back(kv.first);
    }
  }
  
  // Return 
  vector< idx_t >::iterator tmp = std::unique(res.begin(), res.end());
  res.resize( std::distance(res.begin(), tmp) );
  return(res);
}

// Modifies the number of simplices at dimension k by +/- n
void SimplexTree::record_new_simplexes(const idx_t k, const idx_t n){
  if (n_simplexes.size() < k+1){ n_simplexes.resize(k+1); }
  n_simplexes.at(k) += n;
  while(n_simplexes.back() == 0 && n_simplexes.size() > 0){ n_simplexes.resize(n_simplexes.size() - 1); }
}

// Finds the 0-based index of the vertex in the top nodes, or -1 otherwise.
size_t SimplexTree::vertex_index(const idx_t v_id){
  std::map< idx_t, node_ptr > top_vertices = root->children;
  auto it = top_vertices.find(v_id);
  if (it != top_vertices.end()){
    return(std::distance(top_vertices.begin(), it));
  }
  return(-1);
}

// Returns a 0-based integer vector of the first ids available to allocate new vertices
// vector< idx_t > SimplexTree::vertex_available(idx_t n_vertices){
//   std::map< idx_t, node_ptr > top_vertices = root->children;
//   idx_t max = top_vertices.size() + n_vertices;
//   vector< idx_t > new_idx = vector< idx_t >();
//   for (idx_t cc = 0; cc < max; ++cc){
//     if (top_vertices.find(cc) == top_vertices.end()){
//       new_idx.push_back(cc);
//       if (new_idx.size() == n_vertices){ break; }
//     }
//   }
//   return(new_idx);
// }
  
// Emplace v_i new 0-simplices. Ids are automatically assigned. 
// void SimplexTree::add_vertices(const idx_t v_i){
//   const size_t n = root->children.size();
//   for (size_t i = 0; i < v_i; ++i){ 
//     size_t v_id = n+i+1;
//     root->children.emplace(v_id, node_ptr(new node(v_id, root))); 
//   }
//   record_new_simplexes(0, v_i);
// }

// Removes a vertex from the simplex tree, including all edges connected to it.
// void SimplexTree::remove_vertices(vector< idx_t > vertex_ids){
//   if (vertex_ids.size() == 0){ return; }
//   std::map< idx_t, s_ptr<node> > top_vertices = root->children;
//   vector< idx_t > edge_to_remove = { 0, 0 };
//   for (vector< idx_t >::const_iterator vid = vertex_ids.begin(); vid != vertex_ids.end(); ++vid){
//     
//     // First remove any of its cofaces
//     remove_vertex_cofaces(*vid);
//     
//     // Then remove the edges with labels > the query vertex
//     for (auto& top_vertex: top_vertices) {
//       if (top_vertex.first >= *vid){ continue; }
//       edge_to_remove[0] = top_vertex.first;
//       edge_to_remove[1] = *vid;
//       remove_edge(edge_to_remove);
//     }
//     
//     // Remove the vertex itself and its edges with label < the query vertex  
//     if (root->children.find(*vid) != root->children.end()){
//       node_ptr c_node = root->children.at(*vid);
//       int n_children = c_node->children.size(); 
//       c_node->children.clear();
//       record_new_simplexes(1, -n_children);
//       remove_leaf(root, *vid);
//     }
//     
//     // Remove any reference in the level map
//     std::string key = std::to_string(*vid) + "-1";
//     level_map.erase(key);
//   }
// }

// Removes all cofaces containing vertex v
// void SimplexTree::remove_vertex_cofaces(const int v){
//   using simplices = std::map< idx_t, s_ptr<node> >;
//   simplices top_nodes = root->children;
//   simplices::iterator it = top_nodes.find(v);
//   if (it != top_nodes.end()){
//     
//     // First, find the query vertex's connected edges with label > than the vertex 
//     simplices v_children = (*it).second->children;
//     for (auto& child: v_children){
//       
//       // Then get those vertices cousins
//       std::string key = std::to_string(child.first)+"-1";
//       if (level_map.find(key) != level_map.end()){
//         std::vector< node_ptr >& adj_vertices = level_map[key];
//         
//         // If the query vertex is the parent of any such vertices, remove it 
//         std::vector< node_ptr >::iterator c_node; 
//         c_node = std::find_if(adj_vertices.begin(), adj_vertices.end(), [&v](node_ptr vi){
//           return(vi->parent->label == v);
//         });
//         if (c_node != adj_vertices.end()){ adj_vertices.erase(c_node); }
//       }
//     }
//   }
// }


// void SimplexTree::remove_edge( vector<idx_t> edge ){
//   if (edge.size() != 2){ stop("Invalid query. 'remove_edge' takes as input a vector of length 2."); }
//   if (edge[0] > edge[1]){ int tmp = edge[0]; edge[0] = edge[1]; edge[1] = tmp; }
//   bool edge_exists = find_simplex(edge);
//   if (!edge_exists){ return; }
//   else {
//     s_ptr<node> v_ptr = root->children.at(edge[0]);
//     remove_leaf(v_ptr, edge[1]);
//   }
// }

// Recursive helper 
// node_ptr remove(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth){
//   if (i >= n_keys || labels == nullptr || c_node == nullptr){ return nullptr; } // base case + safety checks
//   // Base case: If arrived at the largest k-simplex to remove
//   if (depth == n_keys){
//     // std::string key = std::to_string(labels[j]) + "-" + std::to_string(depth);
//     // size_t has_children = c_node->children.size();
//     // if (!bool(j_exists)){
//     //   // Rcout << "Creating new node " << labels[j] << ", parent: " << c_node->label << std::endl;
//     //   node_ptr new_node = node_ptr(new node(labels[j], c_node));
//     //   insert_child(c_node, new_node, depth);
//     //   if (depth > 0){ // keep track of nodes which share ids at the same depth
//     //     std::string key = std::to_string(labels[j]) + "-" + std::to_string(depth);
//     //     if level_map[key] push_back(new_node);
//     //   }
//     // }
//   }
//   
//   // for (int j = i; j < n_keys; ++j){
//   //   // Rcout << "Recursing with child node: " << labels[j] <<  " of parent: " << c_node->label << " and grandparent: " << (c_node->parent == nullptr ? 0 : c_node->parent->label)  << std::endl;
//   //   node_ptr child_node = c_node->children.at(labels[j]);
//   //   remove(labels, j + 1, n_keys, child_node, depth + 1);
//   // }
// }  



// Remove a child node directly from the parent, if it exists
// This will check that the child is a leaf, and if so, then it will: 
// 1) Remove the child from the parents children map
// 2) Remove the child from the level map 
void SimplexTree::remove_leaf(node_ptr parent, idx_t child_label){
  if (parent == nullptr){ return; }
  size_t child_depth = depth(parent) + 1;
  auto child_it = parent->children.find(child_label);
  if (child_it != parent->children.end()){

    // Assert child is a leaf! 
    node_ptr child = (*child_it).second;
    if (!child->children.empty()){ stop("Tried to remove a non-leaf node! Use remove_subtree instead."); }
    
    // Remove from level map 
    std::string key = std::to_string(child_label) + "-" + std::to_string(child_depth);
    vector< node_ptr >& cousins = level_map[key];
    cousins.erase(std::remove(cousins.begin(), cousins.end(), child), cousins.end());
    
    // Remove from parents children
    parent->children.erase(child_it);
  }
  // Rprintf("removing one %d-simplex\n", child_depth-1);
  record_new_simplexes(child_depth-1, -1);
}
  
// Removes an entire subtree rooted as 'sroot', including 'sroot' itself. This function calls remove_leaf recursively. 
void SimplexTree::remove_subtree(node_ptr sroot){
  if (sroot == nullptr){ return; }
  if (!bool(sroot) || sroot.use_count() <= 1){ return; }
  if (sroot->children.empty()){ remove_leaf(sroot->parent, sroot->label); } 
  else {
    // Remark: make sure to use labels instead of iterator here, otherwise the iterator will be invalidated.
    vector< idx_t > child_labels = get_labels(sroot->children, 0);
    std::for_each(child_labels.begin(), child_labels.end(), [this, &sroot](idx_t label){
      remove_subtree(sroot->children.find(label)->second);
    });
    if (sroot != root){ remove_leaf(sroot->parent, sroot->label); }
  }
}
// First removes all the cofaces of a given simplex, including the simplex itself.
void SimplexTree::remove_simplex(vector< idx_t > simplex){
  node_ptr cn = find(simplex);
  if (cn != nullptr && cn != root){
    vector< node_ptr > cofaces = locate_cofaces(cn); // Note the cofaces contain the simplex itself!
    std::for_each(cofaces.begin(), cofaces.end(), [this](node_ptr coface){
      remove_subtree(coface); // cofaces represent a set of subtrees
    });
  }
}
// void SimplexTree::remove_subtree2(vector< idx_t > simplex){
//   node_ptr cn = find(simplex);
//   if (cn && cn != nullptr){ remove_subtree(cn); }
// }
  
// Inserts a child directly to the parent if the child doesn't already exist.
// Also records the addition of the simplex if the child is indeed created. Otherwise, does nothing. 
node_ptr SimplexTree::insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth){
  if (new_child == nullptr){ return nullptr; }
  std::map<idx_t, node_ptr>& parents_children = c_parent->children;
  std::map<idx_t, node_ptr>::iterator key_lb = parents_children.find(new_child->label); 
  if (key_lb == parents_children.end()) { 
    // node_ptr new_child = node_ptr(new node(child_label, c_parent));
    parents_children.insert(key_lb, std::map<idx_t, node_ptr>::value_type(new_child->label, new_child));
    record_new_simplexes(depth, 1);
    return(new_child);
  }
  return(key_lb->second);// change?
}

// Creates a new set of child nodes for the given parent. createChild check for redundancy with insert. 
// Returns all of the children of the node upon completion. 
// void SimplexTree::add_children(node_ptr c_parent, const std::vector<idx_t>& new_children, idx_t depth){
//   std::for_each(new_children.begin(), new_children.end(), [&](const idx_t& id){
//     node_ptr new_child = node_ptr(new node(id, c_parent));
//     insert_child(c_parent, new_child, depth);
//   });
// }
  
// Rcpp wrapper needs STL for API access
inline void SimplexTree::insert_simplex(std::vector<idx_t> labels){
  std::sort(labels.begin(), labels.end()); // Demand that labels are sorted on insertion! 
  insert((idx_t*) &labels[0], 0, labels.size(), root, 0); // start the recursion from the root
}


void SimplexTree::insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth){
  if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
  // Create a set of (i)-simplexes as children of the current node, if they don't already exist
  size_t child_depth = depth + 1; // depth refers to parent depth, child_depth to child depth
  for (int j = i; j < n_keys; ++j){
    size_t j_exists = c_node->children.count(labels[j]);
    if (!bool(j_exists)){
      // Rcout << "Creating new node " << labels[j] << ", parent: " << c_node->label << std::endl;
      node_ptr new_node = node_ptr(new node(labels[j], c_node));
      insert_child(c_node, new_node, depth);
      if (child_depth > tree_max_depth){ tree_max_depth = child_depth; }
      if (child_depth > 1){ // keep track of nodes which share ids at the same depth
        std::string key = std::to_string(labels[j]) + "-" + std::to_string(child_depth);
        level_map[key].push_back(new_node);
      }
    }
  }
  // Recurse on the subtrees of the current node
  for (int j = i; j < n_keys; ++j){
    // Rcout << "Recursing with child node: " << labels[j] <<  " of parent: " << c_node->label << " and grandparent: " << (c_node->parent == nullptr ? 0 : c_node->parent->label)  << std::endl;
    node_ptr child_node = c_node->children.at(labels[j]);
    insert(labels, j + 1, n_keys, child_node, child_depth);
  }
}

// Rcpp wrapper to the find function
bool SimplexTree::find_simplex(vector< idx_t > simplex){
  std::sort(simplex.begin(), simplex.end());
  if (simplex.size() == 0){ return false; }
  if (simplex.size() == 1){ 
    return(root->children.find(simplex.at(0)) != root->children.end()); 
  } 
  else {
    node_ptr res = find(simplex);
    return(bool(res)); // nullptr conversion to bool guarenteed to be false
  }
}
  
// Overloaded in the case where a single (1-length vector) label is given
node_ptr SimplexTree::find(idx_t label){
  std::map<idx_t, node_ptr>::iterator it = root->children.find(label);
  if (it != root->children.end()){
    return(it->second); 
  } else { return(nullptr); }
}
  
// Given an integer label, searches the tree to see if the simplex exists. If so, the simplex
// (node) is returned, else a nullptr is returned.
node_ptr SimplexTree::find(vector<idx_t> simplex){
  node_ptr c_node = root;
  std::map<idx_t,node_ptr>::iterator node_it;
  for (size_t i = 0; i < simplex.size(); ++i){
    if ((node_it = c_node->children.find(simplex[i])) != c_node->children.end()){
      c_node = node_it->second;
    } else { return nullptr; }
  }
  return(c_node);
}

// Returns the top node associated with a given node
node_ptr SimplexTree::top_node(node_ptr cn){
  if (cn == nullptr || cn == root) { return nullptr; }
  node_ptr res = cn; 
  while (res->parent != root){ res = res->parent; }
  return (res);
}

// Recursively calculate the depth of a node
size_t SimplexTree::depth(node_ptr cn){
  if (cn == nullptr || cn == root){ return 0; }
  size_t d = 1; 
  node_ptr r = cn; // copy 
  while (r->parent != root){ 
    ++d; 
    r = r->parent;
  }
  return d; 
}
  
// Utility to get the maximum height / longest path from any given node.
size_t SimplexTree::max_depth(node_ptr cn){
  size_t max_height = depth(cn); 
  std::for_each(begin_dfs(cn), end_dfs(), [=, &max_height](const node_ptr o){
    size_t c_depth = depth(o); 
    if (c_depth > max_height){ max_height = c_depth; }
  });
  return max_height;
}  
  
// R-wrapper for the whole tree.
void SimplexTree::print_tree(){
  print_subtree(root);
}

// Basic breadth-first printing. Each level is prefixed with '.' <level> number of times, followed by the 
// the ids of the nodes at that breadth-level enclosed within parenthesis, e.g. ..( 2 3 4 ) 
void SimplexTree::print_subtree(node_ptr cn){
  std::map<idx_t, node_ptr>::iterator it;
  for (it = cn->children.begin(); it != cn->children.end(); ++it){ 
    idx_t h = max_depth(it->second)-1; 
    Rcout << it->first << " (h = " << h << "): ";
    for (int i = 1; i <= h; ++i){ 
      for (int j = 1; j <= i; ++j){ Rcout << "."; }
      Rcout << "("; print_level(it->second, i); Rcout << " )";
    }
    Rcout << std::endl;
  }
}


void SimplexTree::print_level(node_ptr cnode, idx_t level){
  if (cnode == nullptr || cnode == NULL) return;
  if (level == 0) { 
    Rprintf(" %d", cnode->label);
    // if (sep) { Rcout << "|"; }
  } else if (level > 0 && (!cnode->children.empty())) {
    std::map<idx_t, node_ptr>::iterator it;
    for (it = cnode->children.begin(); it != cnode->children.end(); ++it){
      print_level(it->second, level-1);
    }
  }
}
  
// Given a map of nodes (e.g. node children) and an offset, retrieves the labels minus the offset
std::vector<idx_t> SimplexTree::get_labels(const std::map<idx_t, node_ptr>& level, const idx_t offset){
  std::vector<idx_t> labels;
  labels.reserve(level.size() - offset);
  std::map<idx_t, node_ptr>::const_iterator it = level.begin(); 
  std::advance(it, offset);
  std::transform (it, level.end(), back_inserter(labels), 
    [&](std::pair<idx_t, node_ptr> const& pair) { return pair.first; });
  return(labels);
}
  
// idx_t SimplexTree::intersection_size(vector<idx_t> v1, vector<idx_t> v2){
//   std::unordered_set<idx_t> s(v1.begin(), v1.end());
//   idx_t res = count_if(v2.begin(), v2.end(), [&](idx_t k) {return s.find(k) != s.end(); });
//   return(res);
// }
// 
// // v1 and v2 should already by sorted!!
// vector<idx_t> SimplexTree::intersection(vector<idx_t> v1, vector<idx_t> v2){
//   vector<idx_t> v3;
//   set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v3));
//   return(v3);
// }

// Experimental k-expansion algorithm.
// Performs an expansion of order k, thus reconstructing a k-skeleton from the 1-skeleton alone.
void SimplexTree::expansion(const idx_t k){
  std::for_each(root->children.begin(), root->children.end(), [&](const std::pair<idx_t, node_ptr>& c_node){
    //idx_t simplex[k];
    //expand(c_node.second->children, k, 0, simplex);
  });
}
  
// Expand operation compares a given 'head' nodes children to its siblings. 
// If they have a non-empty intersection, then the intersection is added as a child to the head node. 
void SimplexTree::expand(std::map<idx_t, node_ptr>& v, const idx_t k, idx_t depth, idx_t* simplex){
  
  // if (v.size() <= 1){ return; } // Current level only has one node; intersection will be empty
  // 
  // // For each child node
  // std::map<idx_t, node_ptr>::const_iterator v_it = v.begin(); 
  // for (v_it = v.begin(); v_it != v.end(); ++v_it){
  //   
  //   // Get the 'head' nodes 
  //   node_ptr rel_head = v_it->second; // *relative* head
  //   node_ptr root_head = find(v_it->first); // *root* head
  //   
  //   // If the root/0-simplex of the head doesn't have children, we're done
  //   if (root_head->children.size() == 0){ return; }
  //   
  //   // Get the (1-offset) siblings of the relative head, and the labels of the children of root head
  //   std::vector<idx_t> siblings = get_labels(v, 1);
  //   std::vector<idx_t> children = get_labels(root_head->children, 0);
  //   std::vector<idx_t> sc_int = intersection(children, siblings);
  //   
  //   IntegerVector sib1 = wrap(siblings);
  //   IntegerVector children1 = wrap(children);
  //   IntegerVector sc_int1 = wrap(sc_int);
  //   Rcout << sib1 << std::endl; 
  //   Rcout << sc_int1 << std::endl; 
  //   if (sc_int.size() > 0){
  //     // Rcout << "Adding children " << sc_int1 << " to node " << v_it->first << " (son of " << v_it->second->parent->label << ")" << std::endl; 
  //     add_children(rel_head, sc_int, depth + 2);
  //     expand(rel_head->children, k, depth + 1, simplex);
  //   }
  // }
}
  
// Exports the 1-skeleton as an adjacency matrix 
IntegerMatrix SimplexTree::as_adjacency_matrix(){
  const size_t n = root->children.size();
  IntegerMatrix res = IntegerMatrix(n, n);
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<idx_t, node_ptr> v){
    const std::map<idx_t, node_ptr> vc = v.second->children;
    std::for_each(vc.begin(), vc.end(), [&](std::pair<idx_t, node_ptr> child){
      const int i = vertex_index(v.first), j = vertex_index(child.first);
      res.at(i, j) = 1; 
      res.at(j, i) = 1; 
    });
  });
  return(res);
}
  
// Exports the 1-skeleton as an adjacency matrix 
List SimplexTree::as_adjacency_list(){
  const size_t n = root->children.size();
  std::unordered_map< std::string, std::vector<idx_t> > res(n); // output
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<idx_t, node_ptr> v){
    const std::map<idx_t, node_ptr> vc = v.second->children;
    std::vector<idx_t> adjacencies = std::vector<idx_t>();
    std::for_each(vc.begin(), vc.end(), [&](std::pair<idx_t, node_ptr> child){
      res[std::to_string(v.first)].push_back(child.first);
      res[std::to_string(child.first)].push_back(v.first);
    });
  });
  return(wrap(res));
}
  
// Exports the 1-skeleton as an edgelist 
IntegerMatrix SimplexTree::as_edge_list(){
  // const size_t n = root->children.size();
  idx_t n_edges = 0; 
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<idx_t, node_ptr> v){
    n_edges += v.second->children.size();
  });
  
  int i = 0; 
  IntegerMatrix res = no_init_matrix(n_edges, 2);
  std::for_each(root->children.begin(), root->children.end(), [&](std::pair<idx_t, node_ptr> v){
    const std::map<idx_t, node_ptr> vc = v.second->children;
    std::for_each(vc.begin(), vc.end(), [&](std::pair<idx_t, node_ptr> child){
      res.row(i++) = IntegerVector::create(v.first, child.first); 
    });
  });
  return(res);
}

// Returns an integer vector of the 0-simplices IDs
IntegerVector SimplexTree::get_vertices(){
  if (n_simplexes.size() == 0){ return(IntegerVector::create()); }
  IntegerVector res = no_init_vector(n_simplexes[0]);
  std::transform(root->children.begin(), root->children.end(), res.begin(), [](std::pair<idx_t, node_ptr> v){
    return(v.first);
  });
  return(res);
}



// Exports the k-skeleton as a list
List SimplexTree::as_list(){
  List res = List();
  vector<idx_t> all = vector<idx_t>(); 
  idx_t d = 1;
  std::for_each(begin_bfs(root), end_bfs(), [this, &d,&res, &all](const node_ptr sigma){
    if ( sigma != nullptr && sigma != root){
      vector<idx_t> si = full_simplex(sigma);
      if (si.size() > d){ 
        const size_t n = all.size() / d;
        IntegerMatrix tmp = IntegerMatrix(n, d);
        for (size_t i = 0; i < n; ++i){
          IntegerVector row = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
          tmp(i, _) = row;
        }
        res.push_back(tmp);
        all.clear(); 
        d = si.size(); 
      }
      all.insert(all.end(), si.begin(), si.end());
    }
  });
  const size_t n = all.size() / d;
  IntegerMatrix tmp = IntegerMatrix(n, d);
  for (size_t i = 0; i < n; ++i){
    IntegerVector row = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
    tmp(i, _) = row;
  }
  res.push_back(tmp);
  return res;
}

// Begin can begin at any subtree in the tree
dfs_iter SimplexTree::begin_dfs(node_ptr cn){
  return dfs_iter(cn);
} 
dfs_iter SimplexTree::end_dfs(){
  return dfs_iter(nullptr);
}

// Begin can begin at any subtree in the tree
bfs_iter SimplexTree::begin_bfs(node_ptr cn){
  return bfs_iter(cn);
} 
bfs_iter SimplexTree::end_bfs(){
  return bfs_iter(nullptr);
}

// Recursive helper to extract the full simplex of a given node
void SimplexTree::full_simplex_r(node_ptr node, vector<idx_t>& res){
  if (node == nullptr || node == root){ return; }
  else {
    res.push_back(node->label);
    if (node->parent == root){ return; }
    else {
      full_simplex_r(node->parent, res);
    }
  }
}

// Given two simplices tau and sigma, checks to see if tau is a face of sigma
bool SimplexTree::is_face(vector<idx_t> tau, vector<idx_t> sigma){
  std::sort(tau.begin(), tau.end());
  std::sort(sigma.begin(), sigma.end());
  return std::includes(sigma.begin(), sigma.end(), tau.begin(), tau.end());
}

// Given a node 'sigma', returns a vector of all the nodes part of the subtree of sigma, 
// including sigma.
vector< node_ptr > SimplexTree::subtree_to_vec(node_ptr sigma){
  vector< node_ptr > subtree_nodes = vector< node_ptr >();
  std::for_each(begin_dfs(sigma), end_dfs(), [&subtree_nodes](const node_ptr tau){
    if (tau->children.empty()){
      subtree_nodes.push_back(tau);
    }
  });
  return(subtree_nodes);
}

// Checks for empty intersection
bool empty_intersection(const vector<idx_t> x, const vector<idx_t> y) {
  vector<idx_t>::const_iterator i = x.begin(), j = y.begin();
  while (i != x.end() && j != y.end()){
    if (*i<*j) ++i; else if(*j<*i) ++j; else return false;
  }
  return true;
}

// Returns a vector of simplices representing the cofaces of a given simplex 'sigma'
// First, all simplices with d > depth(sigma) which end in the same vertex label as sigma are found. 
// Let each of these nodes 'n_j'. There are two condition to test whether n_j is a coface of sigma: 
// 1) n_j is a leaf := n_j is a coface of the current node
// 2) n_j has children := every node in the subtree rooted at n_j is a coface of the current node. 
// In the second case, any node in the subtree rooted at n_j is a coface of sigma. This procedure 
// returns only the roots of these subtrees. 
vector<node_ptr> SimplexTree::locate_cofaces(node_ptr cn){
  size_t h = depth(cn);
  const idx_t j = cn->label;
  vector< idx_t > c_word = full_simplex(cn);
  std::set< node_ptr > cofaces = std::set<node_ptr>();
  cofaces.insert(cn);
  for (size_t i = tree_max_depth; i > h; --i){
    std::string key = std::to_string(j) + "-" + std::to_string(i);
    auto ni = level_map.find(key); 
    if (ni != level_map.end()){ // if exists at some level i > j
      vector<node_ptr> cousins = (*ni).second;
      std::for_each(cousins.begin(), cousins.end(), [this, &c_word, &cofaces](const node_ptr n_j){
        Rcpp::checkUserInterrupt(); // allow user-breaking incase finding cofaces is expensive
        vector<idx_t> word = full_simplex(n_j);
        if (is_face(c_word, word)){ cofaces.insert(n_j); } // insert roots only
      });
    }
  }
  vector< node_ptr > output(cofaces.begin(), cofaces.end()); 
  return output; 
}

// vector<node_ptr> SimplexTree::cofaces(vector<idx_t> simplex){
//   
// }

// Elementary collapse - only capable of collapsing sigma through tau, and only if tau has sigma 
// as its only coface.
void SimplexTree::collapse(node_ptr tau, node_ptr sigma){
  vector< node_ptr > cofaces = locate_cofaces(tau);
  if (cofaces.size() == 1 && cofaces.at(0) == sigma){
    // There are technically two cases, either tau and sigma are both leaves or tau contains 
    // sigma as its unique child. Both can be handled by removing sigma first, then tau. 
    remove_leaf(sigma->parent, sigma->label);
    remove_leaf(tau->parent, tau->label);
    // if (tau->children.empty() && sigma->children.empty()) ...
  } else {
    Function warning("warning");
    warning("Invalid collapse: tau does not have sigma as its only coface.");
  }
}

// Alternative way to specify collapse
// TODO: Add vertex only collapse between vertices not necessarily connected by an edge
// void SimplexTree::collapse(node_ptr u, node_ptr v, node_ptr w){
//   
// }

// R-version
void SimplexTree::collapseR(vector< idx_t > tau, vector< idx_t > sigma){
  node_ptr t = find(tau), s = find(sigma);
  if (t != nullptr && s != nullptr){ collapse(t, s); }
}

// Edge contraction 
void SimplexTree::contract(vector< idx_t > edge){
  if (edge.size() != 2){ stop("Contraction is only possible on edges."); }
  std::set< node_ptr > to_remove;
  std::for_each(begin_dfs(root), end_dfs(), [this, &edge, &to_remove](node_ptr sigma){
    const idx_t va = edge[0], vb = edge[1];
    if (sigma->label == vb){ // only consider simplices which contain v_lb
      vector< idx_t > si = full_simplex(sigma);
      bool includes_a = std::find(si.begin(), si.end(), va) != si.end();
      if (includes_a){ // case 1: sigma includes both v_la and v_lb
        to_remove.insert(sigma); // add whole simplex to remove list, and we're done.
      } else { // case 2: sigma includes v_lb, but not v_la
        // Insert new simplices with v_la --> v_lb, identity otherwise
        std::for_each(begin_dfs(sigma), end_dfs(), [this, va, vb](node_ptr end){
          vector< idx_t > to_insert = full_simplex(end);
          std::replace(to_insert.begin(), to_insert.end(), vb, va);
          insert_simplex(to_insert); // si will be sorted upon insertion
        });
        to_remove.insert(sigma);
      }
    }
  });
  // Remove the simplices containing vb
  std::for_each(to_remove.begin(), to_remove.end(), [this](node_ptr e){
    remove_subtree(e);
  });
}

// Given a node, traces back up to the parent, forming the full simplex.
vector<idx_t> SimplexTree::full_simplex(node_ptr node){
  vector<idx_t> simplex = vector<idx_t>();
  full_simplex_r(node, simplex);
  std::reverse(simplex.begin(), simplex.end());
  return simplex;
}

// Serialize the simplex tree
void SimplexTree::serialize(std::string filename){
  Function saveRDS = Function("saveRDS");
  vector< IntegerVector > minimal;
  std::for_each(begin_dfs(root), end_dfs(), [this, &minimal](node_ptr sigma){
    if (sigma != nullptr && sigma != root && sigma->children.empty()){
      vector< node_ptr > coface_roots = locate_cofaces(sigma);
      if (coface_roots.size() == 1 && coface_roots.at(0) == sigma){
        IntegerVector si = wrap(full_simplex(sigma));
        minimal.push_back(si);
      }
    }
  });
  List res = wrap(minimal);
  saveRDS(_["object"] = res, _["file"] = filename);
}

void SimplexTree::unserialize(std::string filename){
  Function readRDS = Function("readRDS");
  List st = readRDS(_["file"] = filename);
  const size_t n = st.size();
  for (size_t i = 0; i < n; ++i){
    IntegerVector si = st.at(i);
    vector< idx_t > sigma(si.begin(), si.end());
    insert_simplex(sigma);
  }
}

// Link of a simplex
vector< node_ptr > SimplexTree::link(node_ptr sigma){
  vector< idx_t > s = full_simplex(sigma);
  std::set< node_ptr > links; 
  std::for_each(begin_dfs(root), end_dfs(), [this, &links, &s](node_ptr tau){
    vector< idx_t > t = full_simplex(tau); 
    if (t.size() > 0 && empty_intersection(t, s)){
      vector< idx_t > v;
      std::set_union(s.begin(), s.end(), 
                     t.begin(), t.end(), 
                     std::back_inserter(v)); 
      node_ptr v_node = find(v); 
      if (v_node != nullptr){ links.insert(tau); }
    }
  });
  vector< node_ptr > link_res(links.begin(), links.end()); 
  return(link_res);
}

// Generic way to apply function to various types of simplices. 
void SimplexTree::apply(SEXP simp, Function f, std::string type){
  node_ptr sigma;
  if ( Rf_isNull(simp) ){
    sigma = root;
  } else {
    vector< idx_t > si = as< vector< idx_t > >(simp);
    sigma = (si.size() == 0) ? root : find(si);
  }
  if (type == "dfs") {
    for (dfs_iter cn = begin_dfs(sigma); cn != end_dfs(); ++cn){
      vector<idx_t> simplex = full_simplex(*cn);
      f(wrap(simplex));
    }
  }
  else if (type == "bfs") {
    for (bfs_iter cn = begin_bfs(sigma); cn != end_bfs(); ++cn){
      vector<idx_t> simplex = full_simplex(*cn);
      f(wrap(simplex));
    }
  }
  else if (type == "cofaces" || type == "star") {
    if (sigma != nullptr){
      vector< node_ptr > coface_roots = locate_cofaces(sigma);
      vector< node_ptr > cofaces = vector< node_ptr >();
      std::for_each(coface_roots.begin(), coface_roots.end(), [this, &cofaces](const node_ptr coface_root){
        vector< node_ptr > cofaces_tmp = subtree_to_vec(coface_root);
        cofaces.insert(cofaces.end(), cofaces_tmp.begin(), cofaces_tmp.end());
      });
      vector< node_ptr >::iterator co = cofaces.begin();
      for (; co != cofaces.end(); ++co){
        vector<idx_t> simplex = full_simplex(*co);
        f(wrap(simplex));
      }
    }
  } else if (type == "link"){
    if (sigma != nullptr){
      vector< node_ptr > links = link(sigma);
      vector< node_ptr >::iterator li = links.begin();
      for (; li != links.end(); ++li){
        vector<idx_t> simplex = full_simplex(*li);
        f(wrap(simplex));
      }
    }
  } else {
    stop("Iteration 'type' is invalid. Please use one of: dfs, bfs, cofaces, star, link");
  }
}

// Prints the keys for the level map
void SimplexTree::print_cousins(){
  Rcout << "< id >-< depth >: < # of cousins >" << std::endl;
  for(auto kv : level_map) {
    Rcout << kv.first << ": " << kv.second.size() << std::endl;
  } 
}

void SimplexTree::print_simplex(node_ptr cn){
  vector< idx_t > si = full_simplex(cn);
  Rcout << "{ ";
  std::for_each(si.begin(), si.end(), [](const idx_t i){ Rcout << i << " "; });
  Rcout << "}" << std::endl; 
}

void SimplexTree::test(){
  
  // Testing DFS
  // use root or root->children[1]
  // std::for_each(begin_dfs(root), end_dfs(), [=](const node_ptr o){
  //   if (o == root){ return; }
  //   vector<idx_t> simplex = full_simplex(o);
  //   std::for_each(simplex.begin(), simplex.end(), [](const idx_t i){
  //     Rcout << i << " ";
  //   });
  //   Rcout << std::endl;
  // });
  
  // Testing cofaces 
  // vector<idx_t> s = { 4, 5 };
  // node_ptr cn = find(s);
  // if (cn != nullptr){
  //   vector<node_ptr> cofaces = aces(cn);
  //   // Rcout << "# cofaces to (4,5) == " << cofaces.size() << std::endl;  
  //   std::for_each(cofaces.begin(), cofaces.end(), [this](const node_ptr coface){
  //     print_simplex(coface);
  //   });
  // }
  
  // Looking at the ref count
  std::for_each(begin_dfs(root), end_dfs(), [=](const node_ptr& o){
    if (o == root){ return; }
    Rprintf("%d use count: %d, is_leaf: %d\n", o->label, o.use_count(), o->children.empty());
  });
  // Testing removing a subtree
  vector<idx_t> s2 = { 2, 3 };
  remove_subtree(find(s2));

  
  // Testing simplex retrieval
  // dfs_iter it = begin_dfs(root);
  // std::for_each(it, end_dfs(), [&it](const node_ptr o){
  //   // if (o){ Rcout << o->label << std::endl; }
  //   Rcout << "simplex: " << it.c_simplex.size();
  //   std::for_each(it.c_simplex.begin(), it.c_simplex.end(), [](const idx_t i){
  //     Rcout << i << " ";
  //   });
  //   Rcout << std::endl;
  // });
}


// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  Rcpp::class_<SimplexTree>("SimplexTree")
  .constructor()
  .field_readonly("n_simplexes", &SimplexTree::n_simplexes)
  .field_readonly("max_depth", &SimplexTree::tree_max_depth)
  .method( "as_XPtr", &SimplexTree::as_XPtr)
  // .method( "add_vertices", &SimplexTree::add_vertices)
  // .method( "remove_vertices", &SimplexTree::remove_vertices)
  // .method( "vertex_available", &SimplexTree::vertex_available)
  .method( "adjacent_vertices", &SimplexTree::adjacent_vertices)
  .method( "insert_simplex", &SimplexTree::insert_simplex)
  .method( "find_simplex", &SimplexTree::find_simplex)
  // .method( "remove_edge", &SimplexTree::remove_edge)
  .method( "remove_simplex", &SimplexTree::remove_simplex)
  .method( "print_tree", &SimplexTree::print_tree )
  // .method( "print_cofaces", &SimplexTree::print_cofaces )
  // .method( "print_cousins", &SimplexTree::print_cousins )
  // .method( "expansion", &SimplexTree::expansion )
  .method( "as_adjacency_matrix", &SimplexTree::as_adjacency_matrix )
  .method( "as_adjacency_list", &SimplexTree::as_adjacency_list)
  .method( "as_edge_list", &SimplexTree::as_edge_list)
  .method( "as_list", &SimplexTree::as_list)
  // .method( "get_vertices", &SimplexTree::get_vertices)
  // .method( "apply_dfs", &SimplexTree::apply_dfs)
  .method( "collapse", &SimplexTree::collapseR)
  .method( "contract", &SimplexTree::contract)
  .method( "is_face", &SimplexTree::is_face)
  // .method( "remove_subtree", &SimplexTree::remove_subtree2)
  .method( "apply", &SimplexTree::apply)
  .method( "serialize", &SimplexTree::serialize)
  .method( "unserialize", &SimplexTree::unserialize)
  // .method( "test", &SimplexTree::test)
  ;
}

/***R
## Inline tests
stree <- simplextree::simplex_tree()
stree$insert_simplex(c(1, 2, 3))
stree$insert_simplex(c(2, 3, 4, 5))

stree$apply(1L, print, "link")

stree$apply( c(1, 2), function(x){ print(x) }, "bfs")
stree$apply(c(1), function(x){ print(x) }, "cofaces")
stree$remove_simplex(c(2, 3))

stree$remove_simplex(c(1, 2, 3))
stree$print_cofaces(c(1, 2, 3))


# stree$test()

stree$contract(c(1, 3))

stree$insert_simplex(1)
stree$insert_simplex(2)
stree$insert_simplex(3)
stree$insert_simplex(c(1, 2))
stree$insert_simplex(c(2, 3))
stree$insert_simplex(c(1, 3))
stree$insert_simplex(c(1, 2, 3))

*/

