#ifndef SIMPLEXTREE_HPP_
#define SIMPLEXTREE_HPP_

// Szudziks pairing function. Takes as input two unsigned integral types (a, b), and uniquely 
// maps (a, b) to a number c, where c is possibly a different integral type 
template <typename T1, typename T2> 
inline T2 szudzik_pair(T1 a, T1 b){
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  return static_cast<T2>(a >= b ? a * a + a + b : a + b * b);
}

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
  level_map.clear(); 
  n_simplexes.clear();
  tree_max_depth = 0; 
  max_id = 0; 
}

inline std::string SimplexTree::get_id_policy() const{
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
inline size_t SimplexTree::degree(idx_t vid) const{
  node_ptr cn = find_by_id(root->children, vid);
  if (cn == nullptr) { return(0); }
  else {
    size_t res_deg = cn->children.size(); // Labels with id < v 
    // auto it = level_map.find(std::to_string(vid) + "-2"); // Labels with id > v 
    auto it = level_map.find(szudzik_pair< idx_t, size_t >(vid, 2));
    if (it != level_map.end()){ res_deg += (*it).second.size(); }
    return(res_deg);
  }
}

// Returns the degree of a node with a given id
inline vector< size_t > SimplexTree::degree(vector< idx_t > vids) const{
  vector< size_t > res = vector< size_t >();
  for (auto id: vids){
    node_ptr cn = find_by_id(root->children, id);
    if (cn == nullptr) { res.push_back(0); }
    else {
      size_t res_deg = 0;
      // auto it = level_map.find(std::to_string(id) + "-2"); // Labels with id > v 
      auto it = level_map.find(szudzik_pair< idx_t, size_t >(id, 2));
      if (it != level_map.end()){ res_deg += (*it).second.size(); }
      res_deg += cn->children.size(); // Labels with id < v 
      res.push_back(res_deg);
    }
  }
  return(res);
}

// --------- Begin R API --------- 
// These functions are exported through the Rcpp module. 

// Search the level map (cousins) to quickly get the adjacency relations. 
// The set of adjacency relations are the 0-simplexes connected to a given vertex v. 
inline vector< idx_t > SimplexTree::adjacent_vertices(const size_t v) const {
  
  // Resulting vector to return
  vector< idx_t > res = vector< idx_t >(); 
  
  // First extract the vertices which labels > v by checking edges
  //std::string key = std::to_string(v) + "-2";
  auto it = level_map.find(szudzik_pair< idx_t, size_t >(v, 2));
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

// Remove a child node directly from the parent, if it exists
// This will check that the child is a leaf, and if so, then it will: 
// 1) Remove the child from the parents children map
// 2) Remove the child from the level map 
inline void SimplexTree::remove_leaf(node_ptr parent, idx_t child_label){
  if (parent == nullptr){ return; }
  idx_t child_depth = depth(parent) + 1;
  node_ptr child = find_by_id(parent->children, child_label);
  if (child){ // does not equal nullptr
    if (!child->children.empty()){ stop("Tried to remove a non-leaf node! Use remove_subtree instead."); }
    
    // Remove from level map 
    //std::string key = std::to_string(child_label) + "-" + std::to_string(child_depth);
    const size_t key = szudzik_pair< idx_t, size_t >(child_label, child_depth);
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

template <typename Lambda>
void vector_handler(SEXP sigma, Lambda&& f){
  const unsigned int s_type = TYPEOF(sigma);
  const auto check_valid = [](SEXP v) -> bool { 
    IntegerVector iv = v;
    return std::any_of(begin(iv), end(iv), [](int el) -> bool { 
      return(el < 0 || el > std::numeric_limits< idx_t >::max()); 
    });
  };
  if (!Rf_isNull(Rf_getAttrib(sigma, R_DimSymbol))){
    NumericMatrix m = as< NumericMatrix >(sigma);
    const size_t n = m.nrow();
    for (size_t i = 0; i < n; ++i){
      NumericVector cr = m(i,_);
      f(as< vector< idx_t > >(cr));
    }
  } else if (s_type == INTSXP || s_type == REALSXP){
    if (check_valid(sigma)){ stop("Only unsigned integer simplices are supported."); }
    vector< idx_t > simplex = as< vector< idx_t > >(sigma);
    f(simplex);
  } else if (s_type == LISTSXP || s_type == VECSXP){
    List simplices = List(sigma);
    const size_t n = simplices.size(); 
    for (size_t i = 0; i < n; ++i){
      if (check_valid(simplices.at(i))){ stop("Only unsigned integer simplices are supported."); }
      f(as< vector< idx_t > >(simplices.at(i)));
    }
  } else { stop("Unknown type passed, must be list type or vector type."); }
}

// R-facing insert wrapper
inline void SimplexTree::insert(SEXP sigma){
  // RProgress::RProgress pb("Inserting [:bar] ETA: :eta", List(sigma).size());
  // pb.tick(0);
  vector_handler(sigma, [this](vector< idx_t > simplex){
    insert_simplex(simplex);
    // pb.tick();
  });
}
// R-facing remove wrapper
inline void SimplexTree::remove(SEXP sigma){
  vector_handler(sigma, [this](vector< idx_t > simplex){
    remove_simplex(simplex);
  });
}
// R-facing find wrapper
inline LogicalVector SimplexTree::find(SEXP sigma) const {
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
  
  // Require sigma be unique.
  auto it = std::unique(sigma.begin(), sigma.end());
  sigma.resize(std::distance(sigma.begin(), it)); 
  
  // Begin recursive routine
  insert((idx_t*) &sigma[0], 0, sigma.size(), root, 0); // start the recursion from the root
}

inline void SimplexTree::insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth){
  if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
  // Create a set of (i)-simplexes as children of the current node, if they don't already exist
  idx_t child_depth = depth + 1; // depth refers to parent depth, child_depth to child depth
  for (int j = i; j < n_keys; ++j){
    node_ptr cn = find_by_id(c_node->children, labels[j]);
    if (!bool(cn)){ // cn is nullptr
      node_ptr new_node = node_ptr(new node(labels[j], c_node));
      insert_child(c_node, new_node, depth);
      if (child_depth > tree_max_depth){ tree_max_depth = child_depth; }
      if (child_depth > 1){ // keep track of nodes which share ids at the same depth
        //std::string key = std::to_string(labels[j]) + "-" + std::to_string(child_depth);
        const size_t key = szudzik_pair< idx_t, size_t >(labels[j], child_depth);
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
inline std::function<bool(const node_ptr)> SimplexTree::eq_node_id(const idx_t label) const{ 
  return [label](const node_ptr cn)->bool{ return(cn->label == label); };
}

// Finds the 0-based index of the vertex in the top nodes, or -1 otherwise.
inline size_t SimplexTree::vertex_index(const idx_t id) const{
  auto it = std::find_if(begin(root->children), end(root->children), eq_node_id(id));
  return (it == end(root->children) ? -1 : std::distance(begin(root->children), it));
}

// Overloaded in the case where a single (1-length vector) label is given
inline node_ptr SimplexTree::find_by_id(const node_set_t& level, idx_t label) const{
  auto it = std::find_if(begin(level), end(level), eq_node_id(label));
  return it != end(level) ? (*it) : nullptr; 
}

// Wrapper to find a vertex from the top nodes
inline node_ptr SimplexTree::find_vertex(const idx_t v_id) const{
  return find_by_id(root->children, v_id);
}

// Given an integer label, searches the tree to see if the simplex exists. If so, the simplex
// (node) is returned, else a nullptr is returned.
inline node_ptr SimplexTree::find_node(vector< idx_t > simplex) const{
  node_ptr c_node = root;
  const size_t d = simplex.size();
  for (size_t i = 0; i < d; ++i){
    c_node = find_by_id(c_node->children, simplex[i]);
    if (c_node == nullptr){ return nullptr; }
  }
  return(c_node);
}

// Applies find_simplex to all the simplices in the container
inline vector< bool > SimplexTree::find_simplices(vector< vector< idx_t > > simplices) const{
  vector< bool > si_exists(simplices.size());
  const auto fs = [this](const simplex_t sigma) -> bool { return find_simplex(sigma); };
  std::transform(begin(simplices), end(simplices), begin(si_exists), fs);
  return si_exists; 
}

// Finds the simplex 'sigma', specified by an integer vector 
inline bool SimplexTree::find_simplex(vector< idx_t > sigma) const{
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
inline node_ptr SimplexTree::top_node(node_ptr cn) const{
  if (cn == nullptr || cn == root) { return nullptr; }
  node_ptr res = cn; 
  while (res->parent != root){ res = res->parent; }
  return (res);
}

// Recursively calculate the depth of a node
inline size_t SimplexTree::depth(node_ptr cn) const{
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
inline size_t SimplexTree::max_depth(node_ptr cn) const{
  size_t max_height = depth(cn); 
  std::for_each(begin_dfs(cn), end_dfs(), [=, &max_height](const node_ptr o){
    size_t c_depth = depth(o); 
    if (c_depth > max_height){ max_height = c_depth; }
  });
  return max_height;
}  

// R-wrapper for the whole tree.
inline void SimplexTree::print_tree() const{
  print_subtree(root);
}

// Basic breadth-first printing. Each level is prefixed with '.' <level> number of times, followed by the 
// the ids of the nodes at that breadth-level enclosed within parenthesis, e.g. ..( 2 3 4 ) 
inline void SimplexTree::print_subtree(node_ptr cn) const{
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
inline void SimplexTree::print_level(node_ptr cn, idx_t level) const{
  if (cn == nullptr || cn == NULL) return;
  if (level == 0) { Rprintf(" %d", cn->label); } 
  else if (level > 0 && (!cn->children.empty())) {
    for (const node_ptr& ch: cn->children){
      print_level(ch, level-1);
    }
  }
}

// Given a set of nodes (e.g. node children) and an offset, retrieves the labels past the offset
inline vector< idx_t > SimplexTree::get_labels(const node_set_t& level, idx_t offset) const{
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
  // if ((tree_max_depth-1) >= 2){ stop("Can only perform k-expansion on the 1-skeleton."); }
  for (node_ptr cn: root->children){ 
    if (!cn->children.empty()){ expand(cn->children, k-1); }
  }
}

// Expand operation checks A \cap N^+(vj) \neq \emptyset
// If they have a non-empty intersection, then the intersection is added as a child to the head node. 
inline void SimplexTree::expand(node_set_t& c_set, const idx_t k){
  if (k == 0){ return; }
  
  // Traverse the children
  auto siblings = begin(c_set);
  std::advance(siblings, 1);
  vector< node_ptr > intersection;
  for (node_ptr cn: c_set){
    node_ptr top_v = find_vertex(cn->label);
    if (top_v != nullptr && (!top_v->children.empty())){
      
      // Get the intersection
      intersection.clear();
      std::set_intersection(
        siblings, end(c_set), 
        begin(top_v->children), end(top_v->children), 
        std::back_inserter(intersection), 
        [](const node_ptr sib_n, const node_ptr child_n) -> bool {
          return sib_n->label < child_n->label; 
        }
      );
    
      // Insert and recursively expand 
      if (intersection.size() > 0){
        vector< idx_t > sigma = full_simplex(cn);
        sigma.resize(sigma.size() + 1);
        for (auto int_node: intersection){
          sigma[sigma.size() - 1] = int_node->label;
          // insert_simplex(sigma); // slowest version w/ sorting + checking 
          // insert((idx_t*) &sigma[0], 0, sigma.size(), root, 0); // faster version w/o sorting
          insert((idx_t*) &sigma[0], sigma.size()-1, sigma.size(), cn, sigma.size() - 1); // fastest version
        }
        expand(cn->children, k-1); // recurse
      }
    }
    if (siblings != end(c_set)){ ++siblings; }
  }
}

// Exports the 1-skeleton as an adjacency matrix 
inline IntegerMatrix SimplexTree::as_adjacency_matrix() const{
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
inline List SimplexTree::as_adjacency_list() const{
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
inline IntegerMatrix SimplexTree::as_edge_list() const{
  return get_k_simplices(1);
}

// Exports the k-skeleton as a list
inline List SimplexTree::as_list() const{
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
inline dfs_iter SimplexTree::begin_dfs(node_ptr cn) const{
  return dfs_iter(cn);
} 
inline dfs_iter SimplexTree::end_dfs() const{
  return dfs_iter(nullptr);
}

// Begin can begin at any subtree in the tree
inline bfs_iter SimplexTree::begin_bfs(node_ptr cn) const{
  return bfs_iter(cn);
} 
inline bfs_iter SimplexTree::end_bfs() const{
  return bfs_iter(nullptr);
}

// Recursive helper to extract the full simplex of a given node
inline void SimplexTree::full_simplex_r(node_ptr node, vector< idx_t >& res) const{
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
// Assumes tau and sigma are both sorted.
inline bool SimplexTree::is_face(vector< idx_t > tau, vector< idx_t > sigma) const{
  // std::sort(tau.begin(), tau.end());
  // std::sort(sigma.begin(), sigma.end());
  return std::includes(sigma.begin(), sigma.end(), tau.begin(), tau.end());
}

// Checks for empty intersection
inline bool empty_intersection(const vector<idx_t> x, const vector<idx_t> y){
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
inline vector< node_ptr > SimplexTree::locate_cofaces(node_ptr cn) const {
  const size_t h = depth(cn);
  vector< idx_t > c_word = full_simplex(cn);
  set< node_ptr > cofaces = { cn }; // a simplex cofaces include the simplex itself
  for (idx_t i = tree_max_depth; i > h; --i){
    const size_t key = szudzik_pair< idx_t, size_t >(cn->label, i);
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
inline vector< node_ptr > SimplexTree::expand_subtree(node_ptr sigma) const{
  vector< node_ptr > subtree_nodes = vector< node_ptr >();
  std::for_each(begin_dfs(sigma), end_dfs(), [&subtree_nodes](const node_ptr tau){
    subtree_nodes.push_back(tau);
  });
  return(subtree_nodes);
}

// Expand a given vector of subtrees, collected all of the simplices under these trees. 
inline vector< node_ptr > SimplexTree::expand_subtrees(vector< node_ptr > roots) const {
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
vector<T> get_unique(vector<T> v) {
  auto it = std::unique(begin(v), end(v)); 
  v.resize(std::distance(begin(v), it));
  return v;
}

inline void SimplexTree::get_cousins() const{
  Rprintf("< id >-< depth >: < number of cousins >\n");
  for (auto kv: level_map){
    Rprintf("%d: %d\n", kv.first, kv.second.size()); 
  }
}

// Vertex collapse - A vertex collapse, in this sense, is the result of applying a 
// peicewise map f to all vertices sigma \in K, where given a pair (u,v) -> w, 
// f is defined as: 
// f(x) = { (1) w if x in { u, v }, (2) x o.w. }
inline bool SimplexTree::vertex_collapse(node_ptr vp1, node_ptr vp2, node_ptr vt){
  using simplex_v = vector< idx_t >;

  // Get cofaces of each vertex in the free pair
  auto v1_cofaces = expand_subtrees( locate_cofaces(vp1) );
  auto v2_cofaces = expand_subtrees( locate_cofaces(vp2) );
  
  // Lambda factory to create the map vp -> vt where vp is unique
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
  bool tau_is_coface = ::find(begin(cofaces), end(cofaces), tau) != end(cofaces);
  bool sigma_is_coface = ::find(begin(cofaces), end(cofaces), sigma) != end(cofaces);
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
inline vector< idx_t > SimplexTree::connected_components() const{
  
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
inline vector<idx_t> SimplexTree::full_simplex(node_ptr node) const{
  vector<idx_t> simplex = vector< idx_t >();
  full_simplex_r(node, simplex);
  std::reverse(simplex.begin(), simplex.end());
  return simplex;
}

// Serialize the simplex tree
inline vector< vector< idx_t > > SimplexTree::serialize() const{
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
inline void SimplexTree::save(std::string filename) const{
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
  vector< idx_t > v_check(begin(target_ids), end(target_ids));
  std::sort(begin(v_check), end(v_check));
  auto it = std::unique(begin(v_check), end(v_check));
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
inline vector< node_ptr > SimplexTree::link(node_ptr sigma) const{
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
inline void SimplexTree::traverse_dfs_if(node_ptr s, Lambda f, P1 p1, P2 p2) const{
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
inline void SimplexTree::traverse_dfs(node_ptr s, Lambda f) const{
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

// Applies condition BFS, evaluating the first predicate to determine whether to call 
// 'f' on a given element, and applying the second predicate to determine if the children
// of a given element should be added. 
template <typename Lambda, typename P1, typename P2> 
inline void SimplexTree::traverse_bfs_if(node_ptr s, Lambda f, P1 p1, P2 p2) const{
  using d_node = std::pair< node_ptr, idx_t >; // node ptr + depth marker
  
  // Prepare to iteratively do BFS 
  d_node current = std::make_pair(s, depth(s));
  std::queue< d_node > node_queue; 
  node_queue.push(current);
  
  // Also track depth
  while(!node_queue.empty()){
    d_node cn = node_queue.front();
    if (p2(cn.first, cn.second) &&  !cn.first->children.empty()){
      for (auto nv = cn.first->children.begin(); nv != cn.first->children.end(); ++nv){ 
        node_queue.emplace(make_pair(*nv, cn.second+1)); 
      }
    }
    
    if (p1(cn.first, cn.second)){
      f(cn.first, cn.second);
    }
    node_queue.pop(); 
  }
}

template < typename Lambda > 
inline void SimplexTree::traverse_bfs(node_ptr s, Lambda f) const{
  using d_node = std::pair< node_ptr, idx_t >; // node ptr + depth marker
  
  // Prepare to iteratively do BFS 
  d_node current = std::make_pair(s, depth(s));
  std::queue< d_node > node_queue; 
  node_queue.push(current);
  
  // Also track depth
  while(!node_queue.empty()){
    current = node_queue.front();
    if (current.first != nullptr && !current.first->children.empty()){
      for (auto nv = current.first->children.begin(); nv != current.first->children.end(); ++nv){ 
        node_queue.emplace(make_pair(*nv, current.second+1)); 
      }
    }
    f(current.first, current.second);
    node_queue.pop(); 
  }
}

template < typename Lambda >
inline void SimplexTree::traverse_faces_s(node_ptr s, Lambda f) const{
  const simplex_t sigma = full_simplex(s);
  const size_t k = sigma.size(); 
  const auto valid_eval = [k](const node_ptr cn, const size_t d) -> bool { return d <= k + 1; };
  
  // Before recursing into a child, ensure label is at least as big as sigma's
  const auto valid_children = [sigma, k](const node_ptr cn, const size_t d){ 
    return ((d < k + 1) && cn->label >= sigma.at(d-1));
  };
  
  // Perform top-down DFS at each vertex.
  simplex_t tau = sigma;
  for (size_t v_i = 0; v_i < sigma.size(); ++v_i){
    auto v = find_vertex(sigma.at(v_i));
    traverse_dfs_if(v, [this, &sigma, &tau, &f](node_ptr t, idx_t d){
      tau[d-1] = t->label;
      simplex_t c_tau = simplex_t(begin(tau), begin(tau)+d);
      if (is_face(c_tau, sigma)){
        f(c_tau); 
      }
    }, valid_eval, valid_children);
  }
}

template <typename Lambda>
inline void SimplexTree::traverse_facets(node_ptr s, Lambda f) const{

  // Constants
  const simplex_t sigma = full_simplex(s);
  const size_t sigma_depth = sigma.size();
  
  // Conditions for recursion
  const auto valid_eval = [sigma_depth](const node_ptr cn, const size_t d) -> bool { return d <= (sigma_depth - 1); };
  const auto valid_children = [&sigma, sigma_depth](const node_ptr cn, const size_t d){ 
    return (d <= (sigma_depth - 2)) && cn->label >= sigma.at(d-1);
  };
  
  // Facet search 
  if (sigma_depth <= 1){ return; }
  else if (sigma_depth == 2){ 
    // The facet of an edge is just its two vertices
    f(find_vertex(s->label), 1);
    f(s->parent, 1);
    return; 
  } else {
    size_t c_depth = sigma_depth - 1; 
    node_ptr cn = s->parent; 
    simplex_t tau = sigma;
    
    while (cn != root){
      // Check self 
      tau[c_depth-1] = cn->label;
      if (c_depth == sigma_depth-1 && std::includes(begin(sigma), end(sigma), begin(tau), begin(tau)+c_depth)){
        f(cn, c_depth); 
      }
      
      // Look at the siblings 
      auto tn = std::find(begin(cn->parent->children), end(cn->parent->children), cn);
      std::advance(tn,1);
      
      // Check siblings + their children up to facet depth 
      for (; tn != end(cn->parent->children); ++tn){ // TODO: end loop if tn != face
        simplex_t c_tau = full_simplex(*tn);
        c_tau.resize(sigma_depth);
        traverse_dfs_if(*tn, [&c_tau, &sigma, &f](node_ptr t, idx_t d){
          c_tau[d-1] = t->label;
          if ((d == c_tau.size()-1) && std::includes(begin(sigma), end(sigma), begin(c_tau), begin(c_tau)+d)){
            f(t, d); 
          }
        }, valid_eval, valid_children);
      }
      
      // Move up
      cn = cn->parent; 
      c_depth--;
    }
  }
}

template <typename Lambda> 
inline void SimplexTree::traverse_cofaces(node_ptr s, Lambda f) const{
  if (s != nullptr){
    vector< node_ptr > coface_roots = locate_cofaces(s);
    vector< node_ptr > cofaces = expand_subtrees( coface_roots );
    for (node_ptr& co: cofaces){ f(co); }
  }
}

template <typename Lambda> 
inline void SimplexTree::traverse_link(node_ptr s, Lambda f) const{
  if (s != nullptr){
    vector< node_ptr > links = link(s);
    for (node_ptr& li: links){ f(li); }
  }
}

// Applies the lambda function 'f' to every simplex in the k-skeleton given by 'k'
template <typename Lambda> 
inline void SimplexTree::traverse_skeleton(node_ptr s, Lambda f, size_t k) const{
  const auto valid_eval = [k](const node_ptr cn, const size_t d){ return d <= k + 1; };
  const auto valid_children = [k](const node_ptr cn, const size_t d){ return d < k + 1; };
  traverse_bfs_if(s, f, valid_eval, valid_children);
}

// Applies the lambda function 'f' to the maximal faces of the k-skeleton given by 'k'
template <typename Lambda> 
inline void SimplexTree::traverse_max_skeleton(node_ptr s, Lambda f, size_t k) const{
  const auto valid_eval = [k](const node_ptr cn, const size_t d){ return d == k + 1; };
  const auto valid_children = [k](const node_ptr cn, const size_t d){ return d < k + 1; };
  traverse_bfs_if(s, f, valid_eval, valid_children);
}

template <typename Lambda>
inline void SimplexTree::trav_switch(node_ptr sigma, Lambda f, std::string type, List args) const{
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
  } else if(type == "facets"){
    traverse_facets(sigma, f);
  } else { stop("Iteration 'type' is invalid. Please use one of: dfs, bfs, cofaces, star, link, skeleton, or maximal-skeleton"); }
}

// Generic way to apply function to various types of simplices. 
// This acts the generic R-facing version.
inline List SimplexTree::traverse_int(SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args, bool save_res) const{
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

// Prints a simplex
inline void SimplexTree::print_simplex(node_ptr cn) const{
  vector< idx_t > si = full_simplex(cn);
  Rcout << "{ ";
  std::for_each(si.begin(), si.end(), [](const idx_t i){ Rcout << i << " "; });
  Rcout << "}" << std::endl; 
}

// Intermediate struct to enable faster filtration building
struct weighted_simplex {
  node_ptr np;
  double diam; 
  double face_diam;
  size_t dim; 
  // weighted_simplex(node_ptr cn, double d1, double d2, size_t di) 
  //   : np(cn), diam(d1), face_diam(d2), dim(di){}
};

// TODO: rename to flag complex; a Flag complex is a complex built atop a weighted graph. This function 
// essentially constructs a filtration of a weighted flag complex. 
// TODO: Allow two approaches: one independent storage option, which allows arbitrary insertions before and 
// and computing the filtration, or one static option, which *moves* the unique node_ptrs into a container 
// and allows constant-time thresholding, but is invalidated by arbitrary insertions/removals afterwards. 
// Given a dimension k and set of weighted edges (u,v) representing the weights of the ordered edges in the trie, 
// constructs a std::function object which accepts as an argument some weight value 'epsilon' and 
// returns the simplex tree object.  
inline void SimplexTree::rips(vector< double > weights, const size_t k){
  if (weights.size() != this->n_simplexes.at(1)){ stop("Must have one weight per edge."); return; }
  
  // 1. Perform k-expansion.
  expansion(k);

  // 2. Assign edge weights to a map; relies on precondition that 
  size_t i = 0;
  std::map< node_ptr, double > simplex_weights;
  traverse_max_skeleton(root, [&simplex_weights, &weights, &i](const node_ptr cn, const size_t d){
    if (cn != nullptr){
      simplex_weights.emplace(cn, weights.at(i++)); 
    }
  }, 1);

  // 2. Define way getting 'weight' of any simplex.
  std::function< double(const node_ptr, const idx_t) > sigma_weight; 
  sigma_weight = [this, &simplex_weights, &sigma_weight](const node_ptr sigma, const idx_t d) -> double {
    switch (d){ 
      case 0: return(-std::numeric_limits< double >::infinity());
      case 1: return(0.0);
      case 2: return( double(simplex_weights[sigma]) );
      default:
        auto it = simplex_weights.find(sigma);
        if (it != end(simplex_weights)){
          return(it->second); 
        } else {
          double max_weight = 0.0;
          traverse_facets(sigma, [&max_weight, &sigma_weight](const node_ptr& cn, idx_t d){
            const double ew = sigma_weight(cn, d);
            if (ew > max_weight){ max_weight = ew; }
          });
          simplex_weights.emplace_hint(it, sigma, max_weight);
          return(max_weight);
        }
    }
  };

  // 3. Sort simplices by weight in ascending order
  vector< weighted_simplex > w_simplices;
  traverse_skeleton(root, [&w_simplices, &sigma_weight](const node_ptr cn, const size_t d){
    double w = sigma_weight(cn, d); 
    double pw = d == 0 ? -std::numeric_limits< double >::infinity() : sigma_weight(cn->parent, d-1);
    w_simplices.push_back(weighted_simplex{cn, w, pw, d});
  }, k);
  
  // Break weight ties w/ dimension
  std::sort(begin(w_simplices), end(w_simplices), [](const weighted_simplex& s1, const weighted_simplex& s2) -> bool {
    return(s1.diam == s2.diam ? s1.dim < s2.dim : s1.diam < s2.diam);
  });
    
  // 4. Create indexed filtration 
  fc.clear();
  fc.reserve(w_simplices.size());
  for (auto sigma_it = begin(w_simplices); sigma_it != end(w_simplices); ++sigma_it){
    auto& sigma = *sigma_it; 
    indexed_simplex tau;
    tau.label = sigma.np->label;
    tau.index = sigma.diam; 
    // Find the index of sigma's parent, or 0 if sigma if the parent is the empty face. 
    if (sigma.np->parent == root || sigma.np == root){
      tau.parent_idx = 0;
    } else {
      // Search for the lower bound on where sigma's parents weight is
      auto lb = std::lower_bound(begin(w_simplices), sigma_it, sigma.face_diam, 
        [](const weighted_simplex& si, const double& val) -> bool {
        return(si.diam < val);
      });
      const auto sigma_parent = sigma.np->parent;
      auto p_it = std::find_if(lb, sigma_it, [&sigma_parent](const weighted_simplex& si){
        return(si.np == sigma_parent);
      });
      // if (p_it == sigma_it){ stop("sigma detected itself as the parent!"); }
      tau.parent_idx = std::distance(begin(w_simplices), p_it);
    }
    fc.push_back(tau);
  }
  
  // Set state of the filtration to the max
  included = vector< bool >(fc.size(), true);
}
 
// Returns the indices of where the labels that make up the simplex 
// at index 'idx' are in the filtration in ascending order.
inline vector< size_t > SimplexTree::simplex_idx(const size_t idx) const {
  vector< size_t > si_idx;
  size_t c_idx = idx; 
  while(c_idx > 0){
    si_idx.push_back(c_idx);
    c_idx = fc.at(c_idx).parent_idx;
  }
  std::reverse(begin(si_idx), end(si_idx));
  return(si_idx);
}
 
// Given a vector of indices, expand the indices to form the simplex
inline vector< idx_t > SimplexTree::expand_simplex(const vector< size_t > idx) const {
  simplex_t sigma(idx.size());
  //const vector< indexed_simplex >& fc_cache = this->fc;
  std::transform(begin(idx), end(idx), begin(sigma), [this](const size_t i){
    return(fc[i].label);
  });
  return(sigma);
}

// Get index corresponding to eps (inclusive)
inline void SimplexTree::threshold_function(double eps){
  using IS = indexed_simplex; 
  const auto eps_it = std::lower_bound(begin(fc), end(fc), eps, [](const IS s, double val) -> bool {
    return(s.index <= val); 
  });
  const size_t eps_idx = std::distance( begin(fc), (eps_it-1) );
  threshold_index(eps_idx);
}

// Given a function value 'eps', 
inline void SimplexTree::threshold_index(size_t eps_idx){
  if (eps_idx >= fc.size()) { eps_idx = fc.size() - 1; };
    
  // Get current index complex is at in filtration (inclusive)
  const size_t current_idx = std::distance(
    begin(included),
    std::find(included.begin(), included.end(), false)
  )-1;
  
  // If they match, we're done
  if (current_idx == eps_idx){ return; }
 
  // If current index is less than the eps index, for each from the eps index 
  // to i=0, expand each simplex to all cofaces for which it is maximal
  if (current_idx < eps_idx){
    // Add all simplices up to the root
    for (size_t i = eps_idx; i > 0; --i){ // TODO: i > current_idx? 
      if (!included.at(i)){
        // Collect the indices comprising the simplex. Mark all faces as included.
        vector< size_t > inc = simplex_idx(i);
        for (auto inc_idx: inc){
          included[inc_idx] = true; 
        }
        // Insert maximal face into trie
        insert_simplex(expand_simplex(inc));
      }
    }
  // Otherwise if the current index is greater, just remove 
  } else if (current_idx > eps_idx){
    // Rprintf("Current index: %d, eps_idx: %d\n", current_idx, eps_idx);
    for (size_t i = current_idx; i > eps_idx; --i){
      if (included.at(i)){
        // Remove higher order faces, but not lower-level ones
        vector< size_t > c_face_idx = simplex_idx(i);
        
        // Find lowest face that should be excluded
        auto low_idx = std::find_if(begin(c_face_idx), end(c_face_idx), [eps_idx](const size_t face_idx){
          return(face_idx > eps_idx);
        });
        
        // Expand the lowest face and remove. By definition of a filtration, all cofaces of 
        // the lowest face need also be removed. 
        auto maximal_idx = vector< size_t >(c_face_idx.begin(), low_idx+1);
        remove_simplex(expand_simplex(maximal_idx));
        
        // Mark all cofaces as excluded, since they are also removed in the process
        vector< size_t > exclude(low_idx, end(c_face_idx));
        for (auto exc: exclude){
          included[exc] = false;
        }
      }
    }
  }
}

// Returns the current index in the filtration
inline size_t SimplexTree::rips_index() const {
  if (included.size() == 0){ return 0; }
  const size_t current_idx = std::distance(
    begin(included), std::find(begin(included), end(included), false)
  )-1;
  return(current_idx);
}
inline double SimplexTree::rips_epsilon() const {
  if (included.size() == 0){ return -std::numeric_limits< double >::infinity(); }
  return(fc[rips_index()].index);
}

// Returns the simplices in the filtration in a list
inline vector< vector< idx_t > > SimplexTree::rips_simplices() const {
  const size_t n = fc.size();
  vector< vector< idx_t > > simplices(n);
  for (size_t i = 0; i < n; ++i){
    simplices[i] = expand_simplex(simplex_idx(i)); 
  }
  return simplices;
}

// Retrieves the filtration weights
inline vector< double > SimplexTree::rips_weights() const{
  const size_t n = fc.size();
  vector< double > weights = vector< double >(n);
  for (size_t i = 0; i < n; ++i){
    weights[i] = fc[i].index; 
  }
  return weights;
}

// Returns whether the graph is acyclic.
inline bool SimplexTree::is_tree() const{
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

inline IntegerMatrix SimplexTree::get_k_simplices(const size_t k) const{
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
inline vector< idx_t > SimplexTree::get_vertices() const{
  if (n_simplexes.size() == 0){ return vector< idx_t >(0); }
  vector< idx_t > v;
  v.reserve(n_simplexes.at(0));
  for (const node_ptr& cn: root->children){ v.push_back(cn->label); }
  return v;
}
inline IntegerMatrix SimplexTree::get_edges() const{ return get_k_simplices(1); }
inline IntegerMatrix SimplexTree::get_triangles() const{ return get_k_simplices(2); }
inline IntegerMatrix SimplexTree::get_quads() const{ return get_k_simplices(3); }

#endif
