#ifndef SIMPLEXTREE_HPP_
#define SIMPLEXTREE_HPP_

// ============== BEGIN DEFINITIONS ==============
inline SEXP SimplexTree::as_XPtr(){
  Rcpp::XPtr< SimplexTree> p(this, false);
  return(p);
}

// --------- Begin C++ only API --------- 
// These functions are only available through the included header, and thus can only be accessed 
// on the C++ side. R-facing functions are exported through the module. 


// Clear the entire tree by removing all subtrees rooted at the vertices
// Use labels in finding other iterator will be invalidated
inline void SimplexTree::clear(){
  while(!root->children.empty()){ remove_subtree(*root->children.begin()); }
}

inline std::string SimplexTree::get_id_policy(){
  return id_policy == 0 ? std::string("compressed") : std::string("unique");
}

inline void SimplexTree::set_id_policy(std::string policy){
  if (policy == "compressed"){
    id_policy = 0; 
  } else if (policy == "unique"){
    id_policy = 1; 
  } else {
    stop("Invalid ID policy. Must be either 'compressed' or 'unique'."); 
  }
}

// Returns an integer vector of new ids vertex which can be used to insert 0-simplices
// If compress is set to true, the ids are chosen as the first n unoccupied ids found by iterating 
// through the current set of vertices. Otherwise, is compress is false, a maximum id value is 
// maintained, such that new ids generated must exceed that value. 
inline vector< idx_t > SimplexTree::generate_ids(size_t n){
  if (id_policy == 0){
    auto top = root->children;
    vector< idx_t > new_ids = vector< idx_t >();
    idx_t max = top.size() + n;
    for (idx_t cc = 0; cc < max; ++cc){
      if (find_by_id(top, cc) == nullptr){
        new_ids.push_back(cc);
        if (new_ids.size() == n){ break; }
      }
    }
    return(new_ids);
  } else if (id_policy == 1) {
    node_ptr max_v = *std::max_element(begin(root->children), end(root->children), [](const node_ptr v1, const node_ptr v2){
      return v1->label < v2->label; 
    });
    if (max_id < max_v->label){ max_id = max_v->label; }
    vector< idx_t > new_ids(n);
    std::iota(begin(new_ids), end(new_ids), max_id+1); 
    max_id = new_ids.back();
    return(new_ids);
  }
  return vector< idx_t >(0);
}

// Returns the degree of a node with a given id
inline vector< size_t > SimplexTree::degree(vector< idx_t > vids){
  vector< size_t > res = vector< size_t >();
  for (auto id: vids){
    node_ptr cn = find_by_id(root->children, id);
    if (cn == nullptr) { res.push_back(0); }
    else {
      size_t res_deg = 0;
      auto it = level_map.find(std::to_string(id) + "-2"); // Labels with id > v 
      if (it != level_map.end()){ res_deg += (*it).second.size(); }
      res_deg += cn->children.size(); // Labels with id < v 
      res.push_back(res_deg);
    }
  }
  return(res);
}


// // Returns an integer vector of the 0-simplices IDs
// inline vector< idx_t > SimplexTree::get_vertices(){
//   if (n_simplexes.size() == 0){ return vector< idx_t >(); }
//   vector< idx_t > res = vector< idx_t >(n_simplexes.at(0));
//   std::transform(root->children.begin(), root->children.end(), res.begin(), [](std::pair<idx_t, node_ptr> v){
//     return(v.first);
//   });
//   return(res);
// }

// --------- Begin R API --------- 
// These functions are exported through the Rcpp module. 

// Search the level map (cousins) to quickly get the adjacency relations. 
// The set of adjacency relations are the 0-simplexes connected to a given vertex v. 
inline vector< idx_t > SimplexTree::adjacent_vertices(const size_t v){
  
  // Resulting vector to return
  vector< idx_t > res = vector< idx_t >(); 
  
  // First extract the vertices which labels > v by checking edges
  std::string key = std::to_string(v) + "-2";
  std::unordered_map< std::string, vector< node_ptr > >::iterator it = level_map.find(key);
  if (it != level_map.end()){
    vector< node_ptr > cousins = (*it).second;
    for (const node_ptr& cn: cousins){ res.push_back(cn->parent->label); }
  }
  
  // Then get the vertices with labels < v
  node_ptr cn = find_by_id(root->children, v); 
  if (cn != nullptr){
    for (const node_ptr& ch: cn->children){ 
      res.push_back(ch->label); 
    }
  }
  
  // Return 
  vector< idx_t >::iterator tmp = std::unique(res.begin(), res.end());
  res.resize( std::distance(res.begin(), tmp) );
  return(res);
}

// Modifies the number of simplices at dimension k by +/- n
inline void SimplexTree::record_new_simplexes(const idx_t k, const idx_t n){
  if (n_simplexes.size() < k+1){ n_simplexes.resize(k+1); }
  n_simplexes.at(k) += n;
  while(n_simplexes.back() == 0 && n_simplexes.size() > 0){ n_simplexes.resize(n_simplexes.size() - 1); }
}
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
//         vector< node_ptr >& adj_vertices = level_map[key];
//         
//         // If the query vertex is the parent of any such vertices, remove it 
//         vector< node_ptr >::iterator c_node; 
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
inline void SimplexTree::remove_leaf(node_ptr parent, idx_t child_label){
  if (parent == nullptr){ return; }
  size_t child_depth = depth(parent) + 1;
  node_ptr child = find_by_id(parent->children, child_label);
  if (child){ // does not equal nullptr
    if (!child->children.empty()){ stop("Tried to remove a non-leaf node! Use remove_subtree instead."); }
    
    // Remove from level map 
    std::string key = std::to_string(child_label) + "-" + std::to_string(child_depth);
    vector< node_ptr >& cousins = level_map[key];
    cousins.erase(std::remove(cousins.begin(), cousins.end(), child), cousins.end());
    
    // Remove from parents children
    parent->children.erase(child);
    record_new_simplexes(child_depth-1, -1);
  }
}

// Removes an entire subtree rooted as 'sroot', including 'sroot' itself. This function calls remove_leaf recursively. 
inline void SimplexTree::remove_subtree(node_ptr sroot){
  if (sroot == nullptr){ return; }
  if (!bool(sroot) || sroot.use_count() <= 1){ return; }
  if (sroot->children.empty()){ remove_leaf(sroot->parent, sroot->label); } 
  else {
    // Remark: make sure to use labels instead of iterator here, otherwise the iterator will be invalidated.
    vector< idx_t > child_labels = get_labels(sroot->children, 0);
    std::for_each(child_labels.begin(), child_labels.end(), [this, &sroot](idx_t label){
      remove_subtree(find_by_id(sroot->children, label)); 
    });
    if (sroot != root){ remove_leaf(sroot->parent, sroot->label); }
  }
}

// Inserts multiple simplices specified by a container of integer vectors 
inline void SimplexTree::remove_simplices(vector< vector< idx_t > > simplices){
  for (simplex_t& sigma: simplices){ remove_simplex(sigma); }
}

// First removes all the cofaces of a given simplex, including the simplex itself.
inline void SimplexTree::remove_simplex(vector< idx_t > simplex){
  std::sort(simplex.begin(), simplex.end()); // Demand that labels are sorted on insertion! 
  node_ptr cn = find_node(simplex);
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
inline node_ptr SimplexTree::insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth){
  if (new_child == nullptr){ return nullptr; }
  auto& pc = c_parent->children; // parents children
  node_ptr cn = find_by_id(pc, new_child->label);
  if (cn == nullptr) { // child doesn't exist! create it
    pc.insert(new_child);
    record_new_simplexes(depth, 1);
    return(new_child);
  }
  return(cn); 
}

// Creates a new set of child nodes for the given parent. createChild check for redundancy with insert. 
// Returns all of the children of the node upon completion. 
// void SimplexTree::add_children(node_ptr c_parent, const vector<idx_t>& new_children, idx_t depth){
//   std::for_each(new_children.begin(), new_children.end(), [&](const idx_t& id){
//     node_ptr new_child = node_ptr(new node(id, c_parent));
//     insert_child(c_parent, new_child, depth);
//   });
// }
// inline void SimplexTree::insert(IntegerVector simplex){
//   vector< idx_t > sigma(begin(simplex), end(simplex));
//   insert_simplex(sigma);
// }

template <typename ... Args>
constexpr bool return_void(void(Args ...)) { return true; }
template <typename R, typename ... Args>
constexpr bool return_void(R(Args ...)) { return false; }

template <typename Lambda>
void vector_handler(SEXP sigma, Lambda&& f){
  const unsigned int s_type = TYPEOF(sigma);
  const auto check_valid = [](SEXP v) -> bool { 
    IntegerVector iv = v;
    return std::any_of(begin(iv), end(iv), [](int el) -> bool { 
      return(el < 0 || el > std::numeric_limits< idx_t >::max()); 
    });
  };
  if (s_type == INTSXP || s_type == REALSXP){
    if (check_valid(sigma)){ stop("Only unsigned integer simplices are supported."); }
    vector< idx_t > simplex = as< vector< idx_t > >(sigma);
    f(simplex);
  } else if (s_type == LISTSXP || s_type == VECSXP){
    List simplices = List(sigma);
    const size_t n = simplices.size(); 
    for (size_t i = 0; i < n; ++i){
      if (check_valid(simplices.at(i))){ stop("Only unsigned integer simplices are supported."); }
      vector< idx_t > simplex = as< vector< idx_t > >(simplices.at(i));
      f(simplex);
    }
  } else { stop("Unknown type passed, must be list type or vector type."); }
}

// R-facing insert wrapper
inline void SimplexTree::insert(SEXP sigma){
  vector_handler(sigma, [this](vector< idx_t > simplex){
    insert_simplex(simplex);
  });
}
// R-facing remove wrapper
inline void SimplexTree::remove(SEXP sigma){
  vector_handler(sigma, [this](vector< idx_t > simplex){
    remove_simplex(simplex);
  });
}
// R-facing find wrapper
inline LogicalVector SimplexTree::find(SEXP sigma){
  LogicalVector res = LogicalVector();
  vector_handler(sigma, [this, &res](vector< idx_t > simplex){
    res.push_back(find_simplex(simplex));
  });
  return(res);
}

// Inserts multiple simplices specified by a container of integer vectors 
inline void SimplexTree::insert_simplices(vector< vector< idx_t > > simplices){
  for (simplex_t& sigma: simplices){ insert_simplex(sigma); }
}

// Inserts a simplex 'sigma' specified by an integer vector 
inline void SimplexTree::insert_simplex(vector< idx_t > sigma){
  std::sort(sigma.begin(), sigma.end()); // Demand that labels are sorted on insertion! 
  insert((idx_t*) &sigma[0], 0, sigma.size(), root, 0); // start the recursion from the root
}

inline void SimplexTree::insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth){
  if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
  // Create a set of (i)-simplexes as children of the current node, if they don't already exist
  size_t child_depth = depth + 1; // depth refers to parent depth, child_depth to child depth
  for (int j = i; j < n_keys; ++j){
    node_ptr cn = find_by_id(c_node->children, labels[j]);
    if (!bool(cn)){ // cn is nullptr
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
    node_ptr child_node = find_by_id(c_node->children, labels[j]);
    insert(labels, j + 1, n_keys, child_node, child_depth);
  }
}

// Generates a node id equality predicate 
inline std::function<bool(const node_ptr)> SimplexTree::eq_node_id(const idx_t label){ 
  return [label](const node_ptr cn)->bool{ return(cn->label == label); };
}

// Finds the 0-based index of the vertex in the top nodes, or -1 otherwise.
inline size_t SimplexTree::vertex_index(const idx_t id){
  auto it = std::find_if(begin(root->children), end(root->children), eq_node_id(id));
  return (it == end(root->children) ? -1 : std::distance(begin(root->children), it));
}

// Overloaded in the case where a single (1-length vector) label is given
inline node_ptr SimplexTree::find_by_id(const node_set_t& level, idx_t label){
  auto it = std::find_if(begin(level), end(level), eq_node_id(label));
  return it != end(level) ? (*it) : nullptr; 
}

// Wrapper to find a vertex from the top nodes
inline node_ptr SimplexTree::find_vertex(const idx_t v_id){
  return find_by_id(root->children, v_id);
}

// Given an integer label, searches the tree to see if the simplex exists. If so, the simplex
// (node) is returned, else a nullptr is returned.
inline node_ptr SimplexTree::find_node(vector< idx_t > simplex){
  node_ptr c_node = root;
  const size_t d = simplex.size();
  for (size_t i = 0; i < d; ++i){
    c_node = find_by_id(c_node->children, simplex[i]);
    if (c_node == nullptr){ return nullptr; }
  }
  return(c_node);
}

// Applies find_simplex to all the simplices in the container
inline vector< bool > SimplexTree::find_simplices(vector< vector< idx_t > > simplices){
  vector< bool > si_exists(simplices.size());
  const auto fs = [this](const simplex_t sigma) -> bool { return find_simplex(sigma); };
  std::transform(begin(simplices), end(simplices), begin(si_exists), fs);
  return si_exists; 
}

// Finds the simplex 'sigma', specified by an integer vector 
inline bool SimplexTree::find_simplex(vector< idx_t > sigma){
  std::sort(sigma.begin(), sigma.end());
  if (sigma.size() == 0){ return false; }
  if (sigma.size() == 1){ 
    return(bool(find_vertex(sigma.at(0)))); 
  } 
  else {
    node_ptr res = find_node(sigma);
    return(bool(res)); // nullptr conversion to bool guarenteed to be false
  }
}


// Returns the top node associated with a given node
inline node_ptr SimplexTree::top_node(node_ptr cn){
  if (cn == nullptr || cn == root) { return nullptr; }
  node_ptr res = cn; 
  while (res->parent != root){ res = res->parent; }
  return (res);
}

// Recursively calculate the depth of a node
inline size_t SimplexTree::depth(node_ptr cn){
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
inline size_t SimplexTree::max_depth(node_ptr cn){
  size_t max_height = depth(cn); 
  std::for_each(begin_dfs(cn), end_dfs(), [=, &max_height](const node_ptr o){
    size_t c_depth = depth(o); 
    if (c_depth > max_height){ max_height = c_depth; }
  });
  return max_height;
}  

// R-wrapper for the whole tree.
inline void SimplexTree::print_tree(){
  print_subtree(root);
}

// Basic breadth-first printing. Each level is prefixed with '.' <level> number of times, followed by the 
// the ids of the nodes at that breadth-level enclosed within parenthesis, e.g. ..( 2 3 4 ) 
inline void SimplexTree::print_subtree(node_ptr cn){
  for (const node_ptr& ch: cn->children){
    idx_t h = max_depth(ch)-1; 
    Rcout << ch->label << " (h = " << h << "): ";
    for (int i = 1; i <= h; ++i){ 
      for (int j = 1; j <= i; ++j){ Rcout << "."; }
      Rcout << "("; print_level(ch, i); Rcout << " )";
    }
    Rcout << std::endl;
  }
}

// Prints a given level of the tree
inline void SimplexTree::print_level(node_ptr cn, idx_t level){
  if (cn == nullptr || cn == NULL) return;
  if (level == 0) { Rprintf(" %d", cn->label); } 
  else if (level > 0 && (!cn->children.empty())) {
    for (const node_ptr& ch: cn->children){
      print_level(ch, level-1);
    }
  }
}

// Given a set of nodes (e.g. node children) and an offset, retrieves the labels past the offset
inline vector< idx_t > SimplexTree::get_labels(const node_set_t& level, idx_t offset){
  if (level.size() == 0){ return(vector< idx_t >()); }
  if (offset >= level.size()){ return(vector< idx_t >()); }
  vector< idx_t > labels = vector< idx_t >();
  labels.reserve(level.size() - offset);
  auto it = begin(level);
  std::advance(it, offset);
  std::transform(it, end(level), back_inserter(labels), [&](const node_ptr cn) { return cn->label; });
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
// Performs an expansion of order k, reconstructing the k-skeleton flag complex via an in-depth expansion of the 1-skeleton.
inline void SimplexTree::expansion(const idx_t k){
  if ((tree_max_depth-1) >= 2){ stop("Can only perform k-expansion on the 1-skeleton."); }
  for (node_ptr cn: root->children){ 
    if (!cn->children.empty()){ expand(cn->children, k-1); }
  }
}

// Expand operation checks A \cap N^+(vj) \neq \emptyset
// If they have a non-empty intersection, then the intersection is added as a child to the head node. 
inline void SimplexTree::expand(set< node_ptr, ptr_comp< node > >& c_set, const idx_t k){
  if (k == 0){ return; }
  
  // Traverse the children
  auto siblings = begin(c_set);
  std::advance(siblings, 1);
  for (node_ptr cn: c_set){
    node_ptr top_v = find_vertex(cn->label);
    if (top_v != nullptr && (!top_v->children.empty())){
      vector< idx_t > sibling_labels = get_labels(c_set, std::distance(begin(c_set), siblings)); // A 
      vector< idx_t > top_children_labels = get_labels(top_v->children, 0); // N^+(top_v)
      
      // Get the intersection
      vector< idx_t > intersection;
      std::set_intersection(begin(sibling_labels), end(sibling_labels), 
                            begin(top_children_labels), end(top_children_labels), 
                            std::back_inserter(intersection));
    
      // Insert and recursively expand 
      if (intersection.size() > 0){
        for (idx_t new_label: intersection){
          vector< idx_t > sigma = full_simplex(cn);
          sigma.push_back(new_label);
          insert_simplex(sigma);
          expand(cn->children, k-1); // recurse
        }
      }
    }
    if (siblings != end(c_set)){ ++siblings; }
  }
}

// Exports the 1-skeleton as an adjacency matrix 
inline IntegerMatrix SimplexTree::as_adjacency_matrix(){
  const size_t n = root->children.size();
  IntegerMatrix res = IntegerMatrix(n, n);
  
  // Fill in the adjacency matrix
  for (const node_ptr& vi: root->children){
    const size_t i = vertex_index(vi->label);
    for (const node_ptr& vj: vi->children){
      const size_t j = vertex_index(vj->label);
      res.at(i, j) = res.at(j, i) = 1; 
    }
  }
  return(res);
}

// Exports the 1-skeleton as an adjacency list 
inline List SimplexTree::as_adjacency_list(){
  const size_t n = root->children.size();
  std::unordered_map< std::string, vector< idx_t > > res(n);
  for (const node_ptr& vi: root->children){
    for (const node_ptr& vj: vi->children){
      res[std::to_string(vi->label)].push_back(vj->label);
      res[std::to_string(vj->label)].push_back(vi->label);
    }
  }
  return(wrap(res));
}

// Exports the 1-skeleton as an edgelist 
inline IntegerMatrix SimplexTree::as_edge_list(){
  return get_k_simplices(1);
}

// Exports the k-skeleton as a list
inline List SimplexTree::as_list(){
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
inline dfs_iter SimplexTree::begin_dfs(node_ptr cn){
  return dfs_iter(cn);
} 
inline dfs_iter SimplexTree::end_dfs(){
  return dfs_iter(nullptr);
}

// Begin can begin at any subtree in the tree
inline bfs_iter SimplexTree::begin_bfs(node_ptr cn){
  return bfs_iter(cn);
} 
inline bfs_iter SimplexTree::end_bfs(){
  return bfs_iter(nullptr);
}

// Recursive helper to extract the full simplex of a given node
inline void SimplexTree::full_simplex_r(node_ptr node, vector<idx_t>& res){
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
inline bool SimplexTree::is_face(vector<idx_t> tau, vector<idx_t> sigma){
  std::sort(tau.begin(), tau.end());
  std::sort(sigma.begin(), sigma.end());
  return std::includes(sigma.begin(), sigma.end(), tau.begin(), tau.end());
}

// Checks for empty intersection
inline bool empty_intersection(const vector<idx_t> x, const vector<idx_t> y) {
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
inline vector< node_ptr > SimplexTree::locate_cofaces(node_ptr cn){
  const size_t h = depth(cn);
  vector< idx_t > c_word = full_simplex(cn);
  set< node_ptr > cofaces = { cn }; // a simplex cofaces include the simplex itself
  for (size_t i = tree_max_depth; i > h; --i){
    std::string key = std::to_string(cn->label) + "-" + std::to_string(i);
    auto ni = level_map.find(key); 
    if (ni != level_map.end()){ // if exists at some level i > j
      vector< node_ptr > cousins = (*ni).second;
      for (const node_ptr& n_j: cousins){
        Rcpp::checkUserInterrupt(); // allow user-breaking incase finding cofaces is expensive
        vector< idx_t > word = full_simplex(n_j);
        if (is_face(c_word, word)){ cofaces.insert(n_j); } // insert roots only
      }
    }
  }
  vector< node_ptr > output(cofaces.begin(), cofaces.end()); 
  return output; 
}

// Given a node 'sigma', returns a vector of all the nodes part of the subtree of sigma, 
// including sigma.
inline vector< node_ptr > SimplexTree::expand_subtree(node_ptr sigma){
  vector< node_ptr > subtree_nodes = vector< node_ptr >();
  std::for_each(begin_dfs(sigma), end_dfs(), [&subtree_nodes](const node_ptr tau){
    subtree_nodes.push_back(tau);
  });
  return(subtree_nodes);
}

// Expand a given vector of subtrees, collected all of the simplices under these trees. 
inline vector< node_ptr > SimplexTree::expand_subtrees(vector< node_ptr > roots){
  vector< node_ptr > faces = vector< node_ptr >();
  std::for_each(roots.begin(), roots.end(), [this, &faces](const node_ptr subtree_root){
    vector< node_ptr > tmp = expand_subtree(subtree_root);
    faces.insert(faces.end(), tmp.begin(), tmp.end());
  });
  return(faces);
}

inline bool SimplexTree::vertex_collapseR(idx_t v1, idx_t v2, idx_t v3){
  node_ptr vp1 = find_vertex(v1), vp2 = find_vertex(v2), vt = find_vertex(v3);
  if (vp1 == nullptr || vp2 == nullptr || vt == nullptr){ stop("Invalid vertices specified. Not found."); }
  return vertex_collapse(vp1, vp2, vt); // collapse the free pair (vp1, vp2) --> vt
}

template <typename T>
vector<T> get_unique(vector<T> v){
  auto it = std::unique(begin(v), end(v)); 
  v.resize(std::distance(begin(v), it));
  return v;
}

inline void SimplexTree::get_cousins(){
  Rprintf("< id >-< depth >: < number of cousins >\n");
  for (auto kv: level_map){
    Rprintf("%s: %d\n", kv.first.c_str(), kv.second.size()); 
  }
}

// Vertex collapse - A vertex collapse, in this sense, is the result of applying a 
// peicewise map f to all vertices sigma \in K, where given a pair (u,v) -> w, 
// f is defined as: 
// f(x) = { (1) w if x in { u, v }, (2) x o.w. }
inline bool SimplexTree::vertex_collapse(node_ptr vp1, node_ptr vp2, node_ptr vt){
  using simplex_v = vector< idx_t >;

  // Get cofaces of each vertex in the free pair
  vector< node_ptr > v1_cofaces = expand_subtrees( locate_cofaces(vp1) );
  vector< node_ptr > v2_cofaces = expand_subtrees( locate_cofaces(vp2) );
  
  // Lambda that returns a lambda to map vp -> vt and insert resulting unique simplices 
  // into the simplex tree. 
  auto insert_mapped_si = [this](const node_ptr vp, const node_ptr vt){
    return [this, &vp, &vt](const node_ptr coface){
      simplex_v si = full_simplex(coface); // retrieve the full simplex 
      std::replace(begin(si), end(si), vp->label, vt->label); // do the mapping
      insert_simplex(get_unique(si)); // insert the (unique) mapped simplices
    };
  };
  
  // Perform the mapping with both vertices
  std::for_each(begin(v1_cofaces), end(v1_cofaces), insert_mapped_si(vp1, vt));
  std::for_each(begin(v2_cofaces), end(v2_cofaces), insert_mapped_si(vp2, vt));

  // Remove the original pair 
  if (vp1 != vt) { remove_simplex(simplex_v(1, vp1->label)); }
  if (vp2 != vt) { remove_simplex(simplex_v(1, vp2->label)); }
  return true; 
}

// Elementary collapse - only capable of collapsing sigma through tau, and only if tau has sigma 
// as its only coface.
inline bool SimplexTree::collapse(node_ptr tau, node_ptr sigma){
  vector< node_ptr > cofaces = expand_subtrees(locate_cofaces(tau));
  bool tau_is_coface = std::find(cofaces.begin(), cofaces.end(), tau) != cofaces.end();
  bool sigma_is_coface = std::find(cofaces.begin(), cofaces.end(), sigma) != cofaces.end();
  if (cofaces.size() == 2 && (tau_is_coface && sigma_is_coface)){
    // There are technically two cases, either tau and sigma are both leaves or tau contains 
    // sigma as its unique child. Both can be handled by removing sigma first, then tau. 
    remove_leaf(sigma->parent, sigma->label);
    remove_leaf(tau->parent, tau->label);
    return(true);
  } else {
    Function warning("warning");
    warning("Invalid collapse: tau does not have sigma as its only coface.");
  }
  return(false);
}

// R-version
inline bool SimplexTree::collapseR(vector< idx_t > tau, vector< idx_t > sigma){
  std::sort(tau.begin(), tau.end());
  std::sort(sigma.begin(), sigma.end());
  node_ptr t = find_node(tau), s = find_node(sigma);
  if (t != nullptr && s != nullptr){ return collapse(t, s); }
  return false; 
}

// Returns the connected components given by the simplicial complex
inline vector< idx_t > SimplexTree::connected_components(){
  
  // Provide means of mapping vertex ids to index values
  vector< idx_t > v = get_vertices(); // vertices are ordered, so lower_bound is valid
  const auto idx_of = [&v](const idx_t val) { return(std::distance(begin(v), std::lower_bound(begin(v), end(v), val))); };
  
  // Traverse the edges, unioning vertices 
  UnionFind uf = UnionFind(root->children.size());
  traverse_max_skeleton(root, [&idx_of, &uf](const node_ptr cn, const size_t d){
    const idx_t a = idx_of(cn->label), b = idx_of(cn->parent->label);
    uf.Union(a, b);
  }, 1);
  
  // Create the connected components
  size_t i = 0; 
  vector< idx_t > cc = vector< idx_t >(v.size());
  for (idx_t cv: v){ cc[i++] = uf.Find(idx_of(cv)); }
  return(cc);
}

// Alternative way to specify collapse
// TODO: Add vertex only collapse between vertices not necessarily connected by an edge
// void SimplexTree::collapse(node_ptr u, node_ptr v, node_ptr w){
//   
// }

// Edge contraction 
inline void SimplexTree::contract(vector< idx_t > edge){
  if (edge.size() != 2){ stop("Contraction is only possible on edges."); }
  if (!find_simplex(edge)){ stop("Edge not found."); }
  set< node_ptr > to_remove;
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
inline vector<idx_t> SimplexTree::full_simplex(node_ptr node){
  vector<idx_t> simplex = vector<idx_t>();
  full_simplex_r(node, simplex);
  std::reverse(simplex.begin(), simplex.end());
  return simplex;
}

// Serialize the simplex tree
inline vector< vector< idx_t > > SimplexTree::serialize(){
  using simplex_t = vector< idx_t >;
  vector< simplex_t > minimal;
  std::for_each(begin_dfs(root), end_dfs(), [this, &minimal](node_ptr sigma){
    if (sigma != nullptr && sigma != root && sigma->children.empty()){
      vector< node_ptr > coface_roots = locate_cofaces(sigma);
      if (coface_roots.size() == 1 && coface_roots.at(0) == sigma){
        minimal.push_back(full_simplex(sigma));
      }
    }
  });
  return(minimal);
}

// Deserialization (C++ side) 
inline void SimplexTree::deserialize(vector< vector< idx_t > > simplices){
  using simplex_t = vector< idx_t >;
  std::for_each(begin(simplices), end(simplices), [this](const simplex_t& sigma){
    insert_simplex(sigma);
  });
}

// Saves the minimal set of simplices needed to restore the complex to an RDS save
inline void SimplexTree::save(std::string filename){
  using simplex_t = vector< idx_t >;
  Function saveRDS = Function("saveRDS");
  vector< simplex_t > minimal = serialize();
  List res = wrap(minimal);
  saveRDS(_["object"] = res, _["file"] = filename);
}

// Deserialization (R interface) 
// inline void SimplexTree::deserializeR(List simplices){
//   const size_t n = simplices.size();
//   for (size_t i = 0; i < n; ++i){
//     IntegerVector si = as<IntegerVector>(simplices.at(i));
//     vector< idx_t > sigma(si.begin(), si.end());
//     insert_simplex(sigma);
//   }
// }

inline void SimplexTree::load(std::string filename){
  Function readRDS = Function("readRDS");
  List st = readRDS(_["file"] = filename);
  const size_t n = st.size();
  for (size_t i = 0; i < n; ++i){
    IntegerVector si = st.at(i);
    vector< idx_t > sigma(si.begin(), si.end());
    insert_simplex(sigma);
  }
}

inline void SimplexTree::reindex(SEXP target_ids){
  const unsigned int s_type = TYPEOF(target_ids);
  if (s_type == INTSXP || s_type == REALSXP){
    vector< idx_t > t_ids = as< vector< idx_t > >(target_ids);
    reindex(t_ids);
  } else if (s_type == LISTSXP || s_type == VECSXP){
    List t_ids_lst = as< List >(target_ids); 
    CharacterVector nm = t_ids_lst.names();
    if (Rf_isNull(nm) || Rf_length(nm) == 0){ stop("target ids must be named if given as a list."); }
    
    // Do the mapping, then send to regular reindexing function
    const vector< idx_t > base_vids = get_vertices();
    vector< idx_t > new_vids = vector< idx_t >(begin(base_vids), end(base_vids));
    for (size_t i = 0; i < nm.size(); ++i){
      idx_t src_id = std::stoi(as< std::string >(nm.at(i)));
      idx_t tgt_id = as< idx_t >(t_ids_lst.at(i));
      const size_t idx = std::distance(begin(base_vids), std::lower_bound(begin(base_vids), end(base_vids), src_id));
      new_vids.at(idx) = tgt_id;
    }
    reindex(new_vids);
  }
}

inline void SimplexTree::reindex(vector< idx_t > target_ids){
  if (n_simplexes.at(0) != target_ids.size()){ stop("target id vector must match the size of the number of 0-simplices."); }
  vector< vector< idx_t > > minimal = serialize();
  vector< idx_t > vids = get_vertices();
  
  // Check the target ids are unique 
  vector< idx_t > v_check(vids.size());
  auto it = std::unique_copy(begin(target_ids), end(target_ids), begin(v_check));
  if (std::distance(begin(v_check), it) != vids.size()){
    stop("target ids must all unique.");
  }
  
  clear(); // clear the tree now that it's been serialized
  for (simplex_t sigma: minimal){
    const size_t n = sigma.size(); 
    for (size_t i = 0; i < n; ++i){
      const size_t idx = std::distance(begin(vids), std::lower_bound(begin(vids), end(vids), sigma.at(i)));
      sigma.at(i) = target_ids.at(idx);
    }
    insert_simplex(sigma);
  }
}

// Link of a simplex
inline vector< node_ptr > SimplexTree::link(node_ptr sigma){
  vector< idx_t > s = full_simplex(sigma);
  set< node_ptr > links; 
  std::for_each(begin_dfs(root), end_dfs(), [this, &links, &s](node_ptr tau){
    vector< idx_t > t = full_simplex(tau); 
    if (t.size() > 0 && empty_intersection(t, s)){
      vector< idx_t > potential_link;
      std::set_union(s.begin(), s.end(), t.begin(), t.end(), std::back_inserter(potential_link)); 
      node_ptr link_node = find_node(potential_link); 
      if (link_node != nullptr){ links.insert(tau); }
    }
  });
  vector< node_ptr > link_res(links.begin(), links.end()); 
  return(link_res);
}

// Applies condition DFS, evaluating the first predicate to determine whether to call 
// 'f' on a given element, and applying the second predicate to determine if the children
// of a given element should be added. 
template <typename Lambda, typename P1, typename P2> 
inline void SimplexTree::traverse_dfs_if(node_ptr s, Lambda f, P1 p1, P2 p2){
  using d_node = std::pair< node_ptr, idx_t >; 
  
  // Prepare to iteratively do DFS 
  d_node current = std::make_pair(s, depth(s)); 
  std::stack< d_node > node_stack; 
  node_stack.push(current);
  
  // Also track depth
  while (!node_stack.empty()){
    d_node cn = node_stack.top();
    if (p1(cn.first, cn.second)){ 
      f(cn.first, cn.second); // apply function
    }
    node_stack.pop();
    if (p2(cn.first, cn.second) && !cn.first->children.empty()){
      std::for_each(cn.first->children.rbegin(), cn.first->children.rend(), [&node_stack, &cn](const node_ptr ch){
        node_stack.push(std::make_pair(ch, cn.second+1));
      });
    }
  }
}

template <typename Lambda> 
inline void SimplexTree::traverse_dfs(node_ptr s, Lambda f){
  using d_node = std::pair< node_ptr, idx_t >; // node ptr + depth marker
  
  // Prepare to iteratively do DFS 
  d_node current = std::make_pair(s, depth(s)); 
  std::stack< d_node > node_stack; 
  node_stack.push(current);
  
  // Also track depth
  while (!node_stack.empty()){
    current = node_stack.top();
    f(current.first, current.second); // apply function
    node_stack.pop();
    if (!current.first->children.empty()){
      std::for_each(current.first->children.rbegin(), current.first->children.rend(), [&node_stack, &current](const node_ptr ch){
        node_stack.push(make_pair(ch, current.second+1));
      });
    }
  }
}
template <typename Lambda> 
inline void SimplexTree::traverse_bfs(node_ptr s, Lambda f){
  std::for_each(begin_bfs(s), end_bfs(), f);
}
template <typename Lambda> 
inline void SimplexTree::traverse_cofaces(node_ptr s, Lambda f){
  if (s != nullptr){
    vector< node_ptr > coface_roots = locate_cofaces(s);
    vector< node_ptr > cofaces = expand_subtrees( coface_roots );
    for (node_ptr& co: cofaces){ f(co); }
  }
}
template <typename Lambda> 
inline void SimplexTree::traverse_link(node_ptr s, Lambda f){
  if (s != nullptr){
    vector< node_ptr > links = link(s);
    for (node_ptr& li: links){ f(li); }
  }
}

// Applies the lambda function 'f' to every simplex in the k-skeleton given by 'k'
template <typename Lambda> 
inline void SimplexTree::traverse_skeleton(node_ptr s, Lambda f, size_t k){
  const auto valid_eval = [k](const node_ptr cn, const size_t d){ return d <= k + 1; };
  const auto valid_children = [k](const node_ptr cn, const size_t d){ return d < k + 1; };
  traverse_dfs_if(s, f, valid_eval, valid_children);
}

// Applies the lambda function 'f' to the maximal faces of the k-skeleton given by 'k'
template <typename Lambda> 
inline void SimplexTree::traverse_max_skeleton(node_ptr s, Lambda f, size_t k){
  const auto valid_eval = [k](const node_ptr cn, const size_t d){ return d == k + 1; };
  const auto valid_children = [k](const node_ptr cn, const size_t d){ return d < k + 1; };
  traverse_dfs_if(s, f, valid_eval, valid_children);
}

template <typename Lambda>
inline void SimplexTree::trav_switch(node_ptr sigma, Lambda f, std::string type, List args){
  if (type == "dfs") { traverse_dfs(sigma, f); }
  else if (type == "bfs") { traverse_bfs(sigma, f); }
  else if (type == "cofaces" || type == "star") { traverse_cofaces(sigma, f); } 
  else if (type == "link"){ traverse_link(sigma, f); } 
  else if (type == "skeleton" || type == "maximal-skeleton"){
    if (args.size() == 0){ stop("Expecting dimension 'k' to be passed."); }
    vector< std::string > args_str = as< vector< std::string > >(args.names());
    bool contains_k = std::any_of(args_str.begin(), args_str.end(), [](const std::string arg){
      return arg == "k";
    });
    if (!contains_k){ stop("Expecting dimension 'k' to be passed."); }
    size_t k = args["k"];
    if (type == "skeleton"){ traverse_skeleton(sigma, f, k); }
    if (type == "maximal-skeleton"){ traverse_max_skeleton(sigma, f, k); }
  } else { stop("Iteration 'type' is invalid. Please use one of: dfs, bfs, cofaces, star, link, skeleton, or maximal-skeleton"); }
}

// Generic way to apply function to various types of simplices. 
// This acts the generic R-facing version.
inline List SimplexTree::traverse_int(SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args, bool save_res){
  node_ptr sigma = bool(Rf_isNull(simp)) ? root : find_node( as< vector<idx_t> >(simp) );
  if (sigma == nullptr){ return List::create(); }
  
  // Get the arguments
  List actual_args = args.isNotNull() ? List(args) : List();
  // vector< std::string > args_str;as< vector< std::string > >(actual_args.names());
  
  // Store results in a list
  List res = List(); 
  
  // Default eval function: retrieve the full simplex, and call R function with that simplex
  auto default_f = [this, &f, &res](const node_ptr cn, const size_t d = 0){
    vector< idx_t > simplex = full_simplex(cn);
    res.push_back(f(wrap(simplex)));
  };
  
  // Do the traversal, possibly saving the results
  trav_switch(sigma, default_f, type, actual_args);
  return res; 
}

// Traversal overload 1
inline void SimplexTree::traverse(Function f, std::string type){
  // Rcout << "here" << std::endl; 
  traverse_int(R_NilValue, f, type, R_NilValue, false);
}
// Traversal overload 2
inline void SimplexTree::traverse(SEXP simp, Function f, std::string type){
  traverse_int(simp, f, type, R_NilValue, false);
}
// Traversal overload 3
inline void SimplexTree::traverse(SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args){
  traverse_int(simp, f, type, args, false);
}

// List-traversal overload 1
inline List SimplexTree::ltraverse(Function f, std::string type){
  return traverse_int(R_NilValue, f, type, R_NilValue, true);
}
// List-traversal overload 2
inline List SimplexTree::ltraverse(SEXP simp, Function f, std::string type){
  return traverse_int(simp, f, type, R_NilValue, true);
}
// List-traversal overload 3
inline List SimplexTree::ltraverse(SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args){
  return traverse_int(simp, f, type, args, true);
}

// Prints a simplex
inline void SimplexTree::print_simplex(node_ptr cn){
  vector< idx_t > si = full_simplex(cn);
  Rcout << "{ ";
  std::for_each(si.begin(), si.end(), [](const idx_t i){ Rcout << i << " "; });
  Rcout << "}" << std::endl; 
}

// Returns whether the graph is acyclic.
inline bool SimplexTree::is_tree(){
  UnionFind ds = UnionFind(n_simplexes.at(0));

  // Traverse the 1-skeleton, unioning all edges. If any of them are part of the same CC, there is a cycle. 
  bool has_cycle = false; 
  auto valid_eval = [&has_cycle](const node_ptr cn, const size_t d){ return (!has_cycle && d == 2); };
  auto valid_children = [&has_cycle](const node_ptr cn, const size_t d){ return (!has_cycle && d < 2); };
  const vector< idx_t > v = get_vertices();
  const auto index_of = [&v](const idx_t vid) -> size_t{ return std::distance(begin(v), std::find(begin(v), end(v), vid)); };
  
  // Apply DFS w/ UnionFind. If a cycle is detected, no more recursive evaluations are performed. 
  traverse_dfs_if(root, [this, &ds, &has_cycle, &index_of](const node_ptr sigma, const size_t d){
    vector< idx_t > si = full_simplex(sigma);
    idx_t i1 = index_of(si.at(0)), i2 = index_of(si.at(1)); 
    if (ds.Find(i1) == ds.Find(i2)){ has_cycle = true; }
    ds.Union(i1, i2);
  }, valid_eval, valid_children);
  return !has_cycle;
}

// inline bool SimplexTree::is_cycle(vector< idx_t > v){
//   UnionFind ds = UnionFind(vertices.size());
//   
//   // Traverse the 1-skeleton, unioning all edges. If any of them are part of the same CC, there is a cycle. 
//   bool is_cycle = false; 
//   auto valid_eval = [&has_cycle](const node_ptr cn, const size_t d){ return (!is_cycle && d == 2); };
//   auto valid_children = [&has_cycle](const node_ptr cn, const size_t d){ return (!is_cycle && d < 2); };
//   const auto index_of = [&v](const idx_t vid) -> size_t{ return std::distance(begin(v), std::find(begin(v), end(v), vid)); };
//   
//   // Apply DFS w/ UnionFind. If a cycle is detected, no more recursive evaluations are performed. 
//   const size_t n_vertex = vertices.size(); 
//   traverse_dfs_if(root, [this, &ds, &is_cycle, &index_of, &n_vertex](const node_ptr sigma, const size_t d){
//     vector< idx_t > si = full_simplex(sigma);
//     idx_t i1 = index_of(si.at(0)), i2 = index_of(si.at(1)); 
//     if (i1 >= n_vertex || i2 >= n_vertex){ return(); }
//     if (ds.Find(i1) == ds.Find(i2)){ is_cycle = true; }
//     ds.Union(i1, i2);
//   }, valid_eval, valid_children);
//   return !is_cycle;
// }

inline IntegerMatrix SimplexTree::get_k_simplices(const size_t k){
  if (n_simplexes.size() <= k){ return IntegerMatrix(0, k+1); }
  IntegerMatrix res = IntegerMatrix(n_simplexes.at(k), k+1);
  size_t i = 0; 
  traverse_max_skeleton(root, [this, &res, &i](const node_ptr s, const size_t d){
    IntegerVector x = wrap(full_simplex(s));
    res(i++, _) = x;
  }, k);
  return(res);
}

// Retrieve the vertices by their label
inline vector< idx_t > SimplexTree::get_vertices(){
  if (n_simplexes.size() == 0){ return vector< idx_t >(0); }
  vector< idx_t > v;
  v.reserve(n_simplexes.at(0));
  for (const node_ptr& cn: root->children){ v.push_back(cn->label); }
  return v;
}
inline IntegerMatrix SimplexTree::get_edges(){ return get_k_simplices(1); }
inline IntegerMatrix SimplexTree::get_triangles(){ return get_k_simplices(2); }
inline IntegerMatrix SimplexTree::get_quads(){ return get_k_simplices(3); }

#endif