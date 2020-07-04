#ifndef SIMPLEXTREE_HPP_
#define SIMPLEXTREE_HPP_

#include "simplextree.h"

// --------- Begin C++ only API --------- 
// These functions are only available through the included header, and thus can only be accessed 
// on the C++ side. R-facing functions are exported through the module. 
using node_ptr = SimplexTree::node_ptr;
using simplex_t = SimplexTree::simplex_t; 


// Copy constructor
inline SimplexTree::SimplexTree(const SimplexTree& sc) : root(new node(-1, nullptr)), tree_max_depth(0), max_id(0), id_policy(0) {
  auto max_tr = st::maximal< true >(&sc, sc.root.get());
  traverse(max_tr, [this](node_ptr cn, idx_t depth, simplex_t sigma){
    insert_it(begin(sigma), end(sigma), root.get(), 0);
    return true; 
  });
  id_policy = sc.id_policy;
};

// Assignment operator
inline SimplexTree& SimplexTree::operator=(const SimplexTree& sc) {
  auto max_tr = st::maximal< true >(&sc, sc.root.get());
  traverse(max_tr, [this](node_ptr cn, idx_t depth, simplex_t sigma){
    insert_it(begin(sigma), end(sigma), root.get(), 0);
    return true; 
  });
  id_policy = sc.id_policy;
  return *this;
};

// Clear the entire tree by removing all subtrees rooted at the vertices
// Use labels in finding other iterator will be invalidated
inline void SimplexTree::clear(){
  //while(!root->children.empty()){ remove_subtree((*begin(root->children)).get()); }
  root.reset(new node(-1, nullptr));
	level_map.clear(); 
  n_simplexes.clear();
  tree_max_depth = 0; 
  max_id = 0; 
}

inline std::string SimplexTree::get_id_policy() const{
  return id_policy == 0 ? std::string("compressed") : std::string("unique");
}

inline void SimplexTree::set_id_policy(std::string policy){
  if (policy == "compressed"){ id_policy = 0; } 
  else if (policy == "unique"){ id_policy = 1; } 
}

// Returns an integer vector of new ids vertex which can be used to insert 0-simplices
// If compress is set to true, the ids are chosen as the first n unoccupied ids found by iterating 
// through the current set of vertices. Otherwise, is compress is false, a maximum id value is 
// maintained, such that new ids generated must exceed that value. 
inline vector< idx_t > SimplexTree::generate_ids(size_t n){
  if (id_policy == 0){
    vector< idx_t > new_ids = vector< idx_t >();
    idx_t max = root->children.size() + n;
    for (idx_t cc = 0; cc < max && new_ids.size() < n; ++cc){
      if (find_by_id(root->children, cc) == nullptr){
        new_ids.push_back(cc);
      }
    }
    return(new_ids);
  } else if (id_policy == 1) {
		// TODO: get rid of max_lementn 
    auto vid = node_label(*std::max_element(begin(root->children), end(root->children)));
    if (max_id < vid){ max_id = vid; }
    vector< idx_t > new_ids(n);
    std::iota(begin(new_ids), end(new_ids), max_id+1); 
    max_id = new_ids.back();
    return(new_ids);
  }
  return vector< idx_t >(0);
}

// Returns the degree of a node with a given id
inline size_t SimplexTree::degree(idx_t vid) const{
  auto cn = find_by_id(root->children, vid);
  if (cn == nullptr) { return(0); }
  else {
    size_t res_deg = cn->children.size(); // Labels with id < v 
    traverse_cousins(vid, 2, [&res_deg](node_ptr cousin){ res_deg += 1; });
//     auto it = level_map.find(encode_node(vid, 2));
//     if (it != level_map.end()){ 
// 			const auto& cousins = (*it).second;
// 			for (const auto& ch: cousins){
// 				res_deg += node_children(ch).size(); 
// 			}
// 		}
    return(res_deg);
  }
}

// Returns the degree of a node with a given id
// inline vector< size_t > SimplexTree::degree(vector< idx_t > vids) const{
//   vector< size_t > res = vector< size_t >();
//   for (auto id: vids){
//     node_ptr cn = find_by_id(root->children, id);
//     if (cn == nullptr) { res.push_back(0); }
//     else {
//       size_t res_deg = 0;
//       // auto it = level_map.find(std::to_string(id) + "-2"); // Labels with id > v 
//       auto it = level_map.find(encode_node(id, 2));
//       if (it != level_map.end()){ 
// 				const auto& cousins = (*it).second;
// 				for (const auto& ch: cousins){
// 					res_deg += node_children(ch).size(); 
// 				}
// 			}
//       res_deg += node_children(cn).size(); // Labels with id < v 
//       res.push_back(res_deg);
//     }
//   }
//   return(res);
// }

// Search the level map (cousins) to quickly get the adjacency relations. 
// The set of adjacency relations are the 0-simplexes connected to a given vertex v. 
inline vector< idx_t > SimplexTree::adjacent_vertices(const size_t v) const {
  
  // Resulting vector to return
  vector< idx_t > res = vector< idx_t >(); 
  
  // First extract the vertices which labels > v by checking edges
  //std::string key = std::to_string(v) + "-2";
  traverse_cousins(v, 2, [&res](node_ptr cousin){
    res.push_back(node_label(cousin->parent)); 
  });
//   auto it = level_map.find(encode_node(v, 2));
//   if (it != level_map.end()){
// 		const auto& cousins = (*it).second;
//     for (const auto& cn: cousins){ 
// 			res.push_back(node_label(cn->parent)); 
// 		}
//   }
  
  // Then get the vertices with labels < v
  node_ptr cn = find_by_id(root->children, v); 
  if (cn != nullptr){
    for (const auto& ch: node_children(cn)){ 
      res.push_back(node_label(ch)); 
    }
  }
  
  // Return 
  vector< idx_t >::iterator tmp = std::unique(res.begin(), res.end());
  res.resize( std::distance(res.begin(), tmp) );
  return(res);
}

// Modifies the number of simplices at dimension k by +/- n. Shrinks the n_simplexes array as needed. 
inline void SimplexTree::record_new_simplexes(const idx_t k, const idx_t n){
  if (n_simplexes.size() < k+1){ n_simplexes.resize(k+1); }
  n_simplexes.at(k) += n;
  while(n_simplexes.back() == 0 && n_simplexes.size() > 0){ n_simplexes.resize(n_simplexes.size() - 1); }
  tree_max_depth = n_simplexes.size();
}

// Remove a child node directly from the parent, if it exists
// This will check that the child is a leaf, and if so, then it will: 
// 1) Remove the child from the parents children map
// 2) Remove the child from the level map 
inline void SimplexTree::remove_leaf(node_ptr parent, idx_t child_label){
  if (parent == nullptr){ return; }
  const idx_t child_depth = depth(parent) + 1;
	const auto eq_node_id_lambda = [child_label](auto& cn){ return(child_label == node_label(cn)); };
  auto child_it = std::find_if(begin(parent->children), end(parent->children), eq_node_id_lambda);
	if (child_it != end(parent->children)){ 
    // Remove from level map 
    auto child = (*child_it).get(); // copy regular node_ptr
    remove_cousin(child, child_depth);
    // const size_t key = encode_node(child_label, child_depth);
    // auto& cousins = level_map[key];

    // cousins.erase(std::remove(begin(cousins), end(cousins), child), end(cousins));
    
    // // If that was the last cousin in the map, erase the key
    // if (cousins.empty()){ level_map.erase(key); }
    
    // Remove from parents children
		parent->children.erase(child_it);
    record_new_simplexes(child_depth-1, -1);
  }
}

// Removes an entire subtree rooted as 'sroot', including 'sroot' itself; calls 'remove_leaf' recursively. 
inline void SimplexTree::remove_subtree(node_ptr sroot){
  if (sroot == nullptr){ return; }
  if (sroot->children.empty()){ remove_leaf(sroot->parent, sroot->label); }  // remove self 
  else {
    // Remark: make sure to use labels instead of iterator here, otherwise the iterator will be invalidated.
    vector< node_ptr > nc(sroot->children.size());
    std::transform(begin(node_children(sroot)), end(node_children(sroot)), begin(nc), [](auto& u_np){
      return u_np.get(); 
    });
    for (auto cn: nc){
			remove_subtree(find_by_id(node_children(sroot), node_label(cn))); 
		}
    // Remove self
    if (sroot && sroot != root.get()){ remove_leaf(sroot->parent, sroot->label); }
  }
}

// Inserts multiple simplices specified by a container of integer vectors 
inline void SimplexTree::remove_simplices(vector< vector< idx_t > > simplices){
  for (simplex_t& sigma: simplices){ remove_simplex(sigma); }
}

// First removes all the cofaces of a given simplex, including the simplex itself.
inline void SimplexTree::remove_simplex(vector< idx_t > simplex){
  if (simplex.size() == 0){ return; }
  std::sort(simplex.begin(), simplex.end()); // Demand that labels are sorted prior to removal  
  node_ptr cn = find_node(simplex);
  if (cn != nullptr && cn != root.get()){
    auto cr = st::coface_roots(this, cn);
    vector< node_ptr > co_v;
    std::transform(cr.begin(), cr.end(), std::back_inserter(co_v), [](auto& cn){ return(get< 0 >(cn)); });
    for (auto cn: co_v){
      remove_subtree(cn);
    }
  }
};

// First removes all the cofaces of a given simplex, including the simplex itself.
template < typename Iter >
inline void SimplexTree::remove_it(Iter s, Iter e){
  if (s == e){ return; }
  node_ptr cn = find_it(s, e, root.get());
  if (cn != nullptr && cn != root.get()){
    auto cr = st::coface_roots(this, cn);
    SmallVector< node_ptr >::allocator_type::arena_type arena;
    SmallVector< node_ptr > co_v{ arena };
    //vector< node_ptr > co_v;
    std::transform(cr.begin(), cr.end(), std::back_inserter(co_v), [](auto& cn){ return(get< 0 >(cn)); });
    for (auto co_n: co_v){
      remove_subtree(co_n);
    }
  }
};

// Inserts multiple simplices specified by a container of integer vectors 
inline void SimplexTree::insert_simplices(vector< vector< idx_t > > simplices){
  for (simplex_t& sigma: simplices){ insert_simplex(sigma); }
}

// Inserts a simplex 'sigma' specified by an integer vector 
inline void SimplexTree::insert_simplex(vector< idx_t > sigma){
  if (sigma.size() == 0){ return; }
  std::sort(sigma.begin(), sigma.end()); // Demand that labels are sorted on insertion! 
  
  // Require sigma be unique.
  auto it = std::unique(sigma.begin(), sigma.end());
  sigma.resize(std::distance(sigma.begin(), it)); 
  
  // Begin recursive routine
  // insert((idx_t*) &sigma[0], 0, sigma.size(), root.get(), 0); // start the recursion from the root
  insert_it(begin(sigma), end(sigma), root.get(), 0);
}

// Create a set of (i)-simplexes as children of the current node, if they don't already exist
// depth == (depth of c_node)
template< bool lex_order, typename Iter >
inline void SimplexTree::insert_it(Iter s, Iter e, node_ptr c_node, const idx_t depth){
  if (s == e || c_node == nullptr){ return; }
  
  const idx_t child_depth = depth+1;
  std::for_each(s, e, [this, &c_node, child_depth](auto label){
    if constexpr (lex_order){
      auto new_it = c_node->children.emplace_hint(c_node->children.end(), std::make_unique< node >(label, c_node));
      if (child_depth > 1){ // keep track of nodes which share ids at the same depth
        // level_map[encode_node(label, child_depth)].push_back((*new_it).get());
        add_cousin((*new_it).get(), child_depth);
      }
      record_new_simplexes(child_depth-1, 1);
    } else {
      auto it = std::find_if(begin(node_children(c_node)), end(node_children(c_node)), [label](auto& cn){
        return(cn->label == label);
      });
      if (it == end(c_node->children)){
        auto new_it = c_node->children.emplace_hint(it, std::make_unique< node >(label, c_node));
        if (child_depth > 1){ // keep track of nodes which share ids at the same depth
          // level_map[encode_node(label, child_depth)].push_back((*new_it).get());
          add_cousin((*new_it).get(), child_depth);
        }
        record_new_simplexes(child_depth-1, 1);
      }
    }
  });
  
  // Recurse on the subtrees of the current node
  idx_t j = 1;
  std::for_each(s, e, [&](auto label){
    insert_it(std::next(s, j), e, find_by_id(c_node->children, label), depth+1);
    ++j;
  });
}


// Create a set of (i)-simplexes as children of the current node, if they don't already exist
inline void SimplexTree::insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth){
  if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
  idx_t child_depth = depth + 1; // depth refers to parent depth, child_depth to child depth
  for (int j = i; j < n_keys; ++j){
		using ut = decltype(*begin(node_children(c_node)));
		const auto compare_node_id = [&labels, &j](ut& cn) -> bool {
			return(cn->label == labels[j]);
		};
		auto it = std::find_if(begin(node_children(c_node)), end(node_children(c_node)), compare_node_id);
 		// auto it = std::find_if(begin(node_children(c_node)), end(node_children(c_node)), eq_node_id(labels[j]));
    if (it == end(c_node->children)){ // doesn't exist yet  
      auto new_it = c_node->children.emplace_hint(it, std::make_unique< node >(labels[j], c_node));
			record_new_simplexes(depth, 1);
      if (child_depth > tree_max_depth){ tree_max_depth = child_depth; }
      if (child_depth > 1){ // keep track of nodes which share ids at the same depth
        add_cousin((*new_it).get(), child_depth);
        // level_map[encode_node(labels[j], child_depth)].push_back((*new_it).get());
      }
    }
  }
  // Recurse on the subtrees of the current node
  for (int j = i; j < n_keys; ++j){
    node_ptr child_node = find_by_id(c_node->children, labels[j]);
    insert(labels, j + 1, n_keys, child_node, child_depth);
  }
}

// Finds the 0-based index of the vertex in the top nodes, or -1 otherwise.
inline size_t SimplexTree::vertex_index(const idx_t id) const{
	using ut = decltype(*begin(node_children(root)));
	const auto compare_node_id = [&id](ut& cn) -> bool {
		return(cn->label == id);
	};
  auto it = std::find_if(begin(root->children), end(root->children), compare_node_id);
  return (it == end(root->children) ? -1 : std::distance(begin(root->children), it));
}

// Overloaded in the case where a single (1-length vector) label is given
inline SimplexTree::node_ptr SimplexTree::find_by_id(const node_set_t& level, idx_t label) const{
  // level.find(label);
  // auto it = std::find_if(begin(level), end(level), eq_node_id(label));
  auto it = std::lower_bound(begin(level), end(level), label, [](auto& np, auto& id){
    return np->label < id;
  });
  return (it != end(level) && (*it)->label == label) ? (*it).get() : nullptr; 
}

// Wrapper to find a vertex from the top nodes
inline SimplexTree::node_ptr SimplexTree::find_vertex(const idx_t v_id) const{
  return find_by_id(root->children, v_id);
}

// Find iterator version
template< typename Iter >
inline SimplexTree::node_ptr SimplexTree::find_it(Iter s, Iter e, node_ptr cn) const {
  for (; s != e && cn != nullptr; ++s){
    cn = find_by_id(cn->children, *s);
    if (cn == nullptr){ return nullptr; }
  }
  return(cn);
}

// Given an integer label, searches the tree to see if the simplex exists. If so, the simplex
// (node) is returned, else a nullptr is returned.
inline SimplexTree::node_ptr SimplexTree::find_node(vector< idx_t > simplex) const{
  node_ptr c_node = root.get();
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

// template< Iter >
// inline bool SimplexTree::find_it(s, e) const{
//   if (s == e){ return false;  }
//   if (sigma.size() == 1){ 
//     return(bool(find_vertex(sigma.at(0)))); 
//   } 
//   else {
//     node_ptr res = find_node(sigma);
//     return(bool(res)); // nullptr conversion to bool guarenteed to be false
//   }
// }

// Returns the top node associated with a given node
inline SimplexTree::node_ptr SimplexTree::top_node(node_ptr cn) const{
  if (cn == nullptr || cn == root.get()) { return nullptr; }
  node_ptr res = cn; 
  while (res->parent != root.get()){ res = res->parent; }
  return (res);
}

// Recursively calculate the depth of a node
inline size_t SimplexTree::depth(node_ptr cn) const{
  if (cn == nullptr || cn == root.get()){ return 0; }
  size_t d;
  for (d = 1; cn && cn->parent != root.get(); ++d){ 
    cn = cn->parent;
  }
  return d; 
}

// Utility to get the maximum height / longest path from any given node.
inline size_t SimplexTree::max_depth(node_ptr cn) const {
	auto dfs = st::preorder(this, cn);
	auto me = std::max_element(dfs.begin(), dfs.end(), [](auto& n1, auto& n2){
		return(get< 1 >(n1) < get< 1 >(n2));
	});
  return get< 1 >(*me);
}  

// Print the whole tree.
template < typename OutputStream >
inline void SimplexTree::print_tree(OutputStream& os) const {
  print_subtree(os, root.get());
}

template < typename OutputStream >
inline void SimplexTree::print_cousins(OutputStream& os) const {
	auto labels = get_vertices();
	for (idx_t c_depth = 2; c_depth <= tree_max_depth; ++c_depth){
		for (auto &label: labels){
		  if (cousins_exist(label, c_depth)){
		    os << "(last=" << label << ", depth=" << c_depth << "): ";
		    traverse_cousins(label, c_depth, [this, &os](node_ptr cousin){
  				print_simplex(os, cousin, false);
  				os << " ";
  		  });
  		  os << std::endl;
		  }
		}
	}
};

// Basic breadth-first printing. Each level is prefixed with '.' <level> number of times, followed by the 
// the ids of the nodes at that breadth-level enclosed within parenthesis, e.g. ..( 2 3 4 ) 
template < typename OutputStream >
inline void SimplexTree::print_subtree(OutputStream& os, node_ptr cn) const {
  for (const auto& ch: cn->children){
    idx_t h = max_depth(ch.get())-1; 
    os << ch->label << " (h = " << h << "): ";
    for (int i = 1; i <= h; ++i){ 
      for (int j = 1; j <= i; ++j){ os << "."; }
      os << "("; 
      print_level(os, ch.get(), i); 
      os << " )";
    }
    os << std::endl;
  }
}

// Prints a given level of the tree
template < typename OutputStream >
inline void SimplexTree::print_level(OutputStream& os, node_ptr cn, idx_t level) const{
  if (cn == nullptr || cn->parent == nullptr) return;
  if (level == 0) { os << " " << cn->label; } 
  else if (level > 0 && (!cn->children.empty())) {
    for (const auto& ch: cn->children){
      print_level(os, ch.get(), level-1);
    }
  }
}

// Prints an individual simplex
template < typename OutputStream >
inline void SimplexTree::print_simplex(OutputStream& os, node_ptr cn, bool newline) const {
  simplex_t si = full_simplex(cn);
  os << "{ ";
  std::for_each(si.begin(), si.end(), [&os](const idx_t i){ os << i << " "; });
  os << "}";
	if (newline){ os << std::endl; }
}

// Given a set of nodes (e.g. node children) and an offset, retrieves the labels past the offset
inline vector< idx_t > SimplexTree::get_labels(const node_set_t& level, idx_t offset) const{
  if (level.size() == 0){ return(vector< idx_t >()); }
  if (offset >= level.size()){ return(vector< idx_t >()); }
  vector< idx_t > labels = vector< idx_t >();
  labels.reserve(level.size() - offset);
  auto it = begin(level);
  std::advance(it, offset);
  std::transform(it, end(level), back_inserter(labels), [&](auto& cn) { return node_label(cn); });
  return(labels);
}

// Performs an expansion of order k, reconstructing the k-skeleton flag complex via an in-depth expansion of the 1-skeleton.
inline void SimplexTree::expansion(const idx_t k){
  expansion_f(k, [this](auto parent, idx_t depth, idx_t label){
    std::array< idx_t, 1 > int_label = { label };
    insert_it(begin(int_label), end(int_label), parent, depth);
  });
}

template < typename Lambda >
inline void SimplexTree::expansion_f(const idx_t k, Lambda&& f){
  for (auto& cn: node_children(root)){ 
    if (!node_children(cn).empty()){ 
			expand_f(cn->children, k-1, 2, f); 
		}
  }
}

// Expand operation checks A \cap N^+(vj) \neq \emptyset
// If they have a non-empty intersection, then the intersection is added as a child to the head node. 
template < typename Lambda >
inline void SimplexTree::expand_f(node_set_t& c_set, const idx_t k, size_t depth, Lambda&& f){
  if (k == 0 || c_set.empty()){ return; }
  // Traverse the children
  auto siblings = std::next(begin(c_set), 1);
  SmallVector< node_ptr >::allocator_type::arena_type arena1;
  SmallVector< node_ptr > intersection { arena1 };
  for (auto& cn: c_set){
    node_ptr top_v = find_vertex(cn->label);
    if (top_v != nullptr && (!top_v->children.empty())){
      
			// Temporary 
			SmallVector< node_ptr >::allocator_type::arena_type arena2;
			SmallVector< node_ptr > sib_ptrs { arena2 } ;
			std::transform(siblings, end(c_set), std::back_inserter(sib_ptrs), [](const auto& n){
				return n.get();
			});

      // Get the intersection
      intersection.clear();
      std::set_intersection(
        begin(sib_ptrs), end(sib_ptrs),
        begin(top_v->children), end(top_v->children), 
        std::back_inserter(intersection), 
        [](auto& sib_n, auto& child_n) -> bool {
          return node_label(sib_n) < node_label(child_n); 
        }
      );
    
      // Insert and recursively expand 
      if (intersection.size() > 0){
        for (auto& int_node: intersection){
          auto face = find_by_id(cn->children, int_node->label);
          if (face == nullptr){
            f(cn.get(), depth, int_node->label);
          }
        }
        expand_f(cn->children, k-1, depth+1, f); // recurse
      }
    }
    if (siblings != end(c_set)){ ++siblings; }
  }
}


inline void SimplexTree::reindex(vector< idx_t > target_ids){
  if (n_simplexes.at(0) != target_ids.size()){ throw std::invalid_argument("target id vector must match the size of the number of 0-simplices."); }
  if (!std::is_sorted(begin(target_ids), end(target_ids))){ throw std::invalid_argument("target ids must be ordered."); }
  if (std::unique(begin(target_ids), end(target_ids)) != end(target_ids)){ throw std::invalid_argument("target ids must all unique."); }

  // Create the map between vertex ids -> target ids
  auto id_map = std::map< idx_t, idx_t >();
  auto vertex_ids = get_vertices();
  for (size_t i = 0; i < vertex_ids.size(); ++i){
    id_map.emplace_hint(end(id_map), vertex_ids[i], target_ids[i]);
  }

  // Apply the map
  auto tr = st::preorder< false >(this);
  st::traverse(tr, [&id_map](node_ptr cn, idx_t depth){
    cn->label = id_map[cn->label];
    return true;
  });

  // Remap the cousins
  for (size_t d = 2; d < tree_max_depth; ++d){
    auto& cousins = level_map.at(depth_index(d));
    for (auto v_id: vertex_ids){
      auto it = cousins.find(v_id);
      if (it != cousins.end()){
        auto kv = std::make_pair(std::move(it->first), std::move(it->second)); // copy 
        cousins.erase(it);
        kv.first = id_map[v_id];
        cousins.insert(kv);
      }
      // TODO: Include this when support for std::extract improves
      // auto node = cousins.extract(idx_t(v_id));
      // if (!node.empty()){
      //   node.key() = id_map[v_id];
      //   cousins.insert(std::move(node));
      // }
    }
  }
}


// Expand operation checks A \cap N^+(vj) \neq \emptyset
// If they have a non-empty intersection, then the intersection is added as a child to the head node. 
// inline void SimplexTree::expand(node_set_t& c_set, const idx_t k){
//   if (k == 0){ return; }
//   // Traverse the children
//   auto siblings = ++begin(c_set);
//   SmallVector< node_ptr >::allocator_type::arena_type arena1;
//   SmallVector< node_ptr > intersection { arena1 };
//   for (auto& cn: c_set){
//     node_ptr top_v = find_vertex(cn->label);
//     if (top_v != nullptr && (!top_v->children.empty())){
//       
// 			// Temporary 
// 			SmallVector< node_ptr >::allocator_type::arena_type arena2;
// 			SmallVector< node_ptr > sib_ptrs { arena2 } ;
// 			std::transform(siblings, end(c_set), std::back_inserter(sib_ptrs), [](const auto& n){
// 				return n.get();
// 			});
// 
//       // Get the intersection
//       intersection.clear();
//       std::set_intersection(
//         begin(sib_ptrs), end(sib_ptrs),
//         begin(top_v->children), end(top_v->children), 
//         std::back_inserter(intersection), 
//         [](auto& sib_n, auto& child_n) -> bool {
//           return node_label(sib_n) < node_label(child_n); 
//         }
//       );
//     
//       // Insert and recursively expand 
//       if (intersection.size() > 0){
//         std::array< idx_t, 1 > int_label;
//         for (auto& int_node: intersection){
//           int_label[0] = int_node->label;
//           insert_it(begin(int_label), end(int_label), cn.get(), depth(cn.get()));
//         }
//         expand(cn->children, k-1); // recurse
//       }
//     }
//     if (siblings != end(c_set)){ ++siblings; }
//   }
// }

// Given two simplices tau and sigma, checks to see if tau is a face of sigma
// Assumes tau and sigma are both sorted.
inline bool SimplexTree::is_face(vector< idx_t > tau, vector< idx_t > sigma) {
  // std::sort(tau.begin(), tau.end());
  // std::sort(sigma.begin(), sigma.end());
  return std::includes(sigma.begin(), sigma.end(), tau.begin(), tau.end());
}

// Returns a vector of simplices representing the cofaces of a given simplex 'sigma'
// First, all simplices with d > depth(sigma) which end in the same vertex label as sigma are found. 
// Let each of these nodes 'n_j'. There are two condition to test whether n_j is a coface of sigma: 
// 1) n_j is a leaf := n_j is a coface of the current node
// 2) n_j has children := every node in the subtree rooted at n_j is a coface of the current node. 
// In the second case, any node in the subtree rooted at n_j is a coface of sigma. 
// Note that this procedure returns only the roots of these subtrees. 
// inline vector< SimplexTree::node_ptr > SimplexTree::locate_cofaces(node_ptr cn) const {
// 	vector< idx_t > c_word = full_simplex(cn);
// 	const size_t h = c_word.size();
//   set< node_ptr > cofaces = { cn }; // a simplex cofaces include the simplex itself
//   for (idx_t i = tree_max_depth; i > h; --i){
// 		for (auto& n_j: node_cousins(cn, i)){
// 			if (is_face(c_word, full_simplex(n_j))){ 
// 				cofaces.insert(n_j); // insert roots only
// 			} 
// 		}
// 	}
//   vector< node_ptr > output(cofaces.begin(), cofaces.end()); 
//   return output; 
// }

// Given a node 'sigma', returns a vector of all the nodes part of the subtree of sigma, 
// including sigma.
// inline vector< SimplexTree::node_ptr > SimplexTree::expand_subtree(node_ptr sigma) const {
//   vector< node_ptr > subtree_nodes = vector< node_ptr >();
// 	for (auto& node: dfs< false >(this, sigma)){
// 		subtree_nodes.push_back(get< 0 >(node));
// 	}
//   return(subtree_nodes);
// }

// Expand a given vector of subtrees, collected all of the simplices under these trees. 
// inline vector< node_ptr > SimplexTree::expand_subtrees(vector< node_ptr > roots) const {
//   vector< node_ptr > faces = vector< node_ptr >();
// 	for (auto& subtree_root: roots){
// 		vector< node_ptr > tmp = expand_subtree(subtree_root);
//     faces.insert(faces.end(), tmp.begin(), tmp.end());
// 	}
//   return(faces);
// }

inline bool SimplexTree::vertex_collapse(idx_t v1, idx_t v2, idx_t v3){
  node_ptr vp1 = find_vertex(v1), vp2 = find_vertex(v2), vt = find_vertex(v3);
  return vertex_collapse(vp1, vp2, vt); // collapse the free pair (vp1, vp2) --> vt
}

template <typename T>
vector<T> get_unique(vector<T> v) {
  auto it = std::unique(begin(v), end(v)); 
  v.resize(std::distance(begin(v), it));
  return v;
}

// Vertex collapse - A vertex collapse, in this sense, is the result of applying a 
// peicewise map f to all vertices sigma \in K, where given a pair (u,v) -> w, 
// f is defined as: 
// f(x) = { (1) w if x in { u, v }, (2) x o.w. }
inline bool SimplexTree::vertex_collapse(node_ptr vp1, node_ptr vp2, node_ptr vt){
	// Lambda to do the mapping
	const auto map_collapse = [vt](simplex_t& si, node_ptr vp){
		std::replace(begin(si), end(si), vp->label, vt->label); 
		return(get_unique(si));
	};
	std::vector< simplex_t > to_insert; 
	for (auto& cn: st::cofaces< true >(this, vp1)){ to_insert.push_back(map_collapse(get< 2 >(cn), vp1)); }
	for (auto& cn: st::cofaces< true >(this, vp2)){ to_insert.push_back(map_collapse(get< 2 >(cn), vp2)); }
	for (auto& si: to_insert){
		insert_simplex(si);
	}
	if (vp1 != vt) { remove_simplex(simplex_t(1, vp1->label)); }
  if (vp2 != vt) { remove_simplex(simplex_t(1, vp2->label)); }

  // // Get cofaces of each vertex in the free pair
	// auto v1_cofaces = generate(cofaces(this, vp1), [](node_ptr cn, idx_t depth){ return cn; });
	// auto v2_cofaces = generate(cofaces(this, vp2), [](node_ptr cn, idx_t depth){ return cn; });
  
  // // Lambda factory to create the map vp -> vt where vp is unique
  // auto insert_mapped_si = [this](const node_ptr vp, const node_ptr vt){
  //   return [this, &vp, &vt](const node_ptr coface){
  //     simplex_v si = full_simplex(coface); // retrieve the full simplex 
  //     std::replace(begin(si), end(si), vp->label, vt->label); // do the mapping
  //     insert_simplex(get_unique(si)); // insert the (unique) mapped simplices
  //   };
  // };
  
  // // Perform the mapping with both vertices
  // std::for_each(begin(v1_cofaces), end(v1_cofaces), insert_mapped_si(vp1, vt));
  // std::for_each(begin(v2_cofaces), end(v2_cofaces), insert_mapped_si(vp2, vt));

  // // Remove the original pair 

  return true; 
}

// Elementary collapse - only capable of collapsing sigma through tau, and only if tau has sigma 
// as its only coface. There are technically two cases, either tau and sigma are both leaves or 
// tau contains sigma as its unique child. Both can be handled by removing sigma first, then tau.
inline bool SimplexTree::collapse(node_ptr tau, node_ptr sigma){
  // vector< node_ptr > cofaces = expand_subtrees(locate_cofaces(tau));
  // bool tau_is_coface = std::find(begin(cofaces), end(cofaces), tau) != end(cofaces);
  // bool sigma_is_coface = std::find(begin(cofaces), end(cofaces), sigma) != end(cofaces);
  // if (cofaces.size() == 2 && (tau_is_coface && sigma_is_coface)){ 
  //   remove_leaf(sigma->parent, sigma->label);
  //   remove_leaf(tau->parent, tau->label);
  //   return(true);
  // } 
  // return(false);
	auto tau_cofaces = st::cofaces(this, tau);
	bool sigma_only_coface = false; 
	size_t n_cf = 0; 
	traverse(tau_cofaces, [&tau, &n_cf, &sigma_only_coface](node_ptr coface, idx_t depth){
		sigma_only_coface |= (coface == tau);
		return((++n_cf) >= 2); 
	});
	if (sigma_only_coface){
		remove_leaf(sigma->parent, sigma->label);
  	remove_leaf(tau->parent, tau->label);
		return(true);
	}
	return(false);

}

inline bool SimplexTree::collapse(simplex_t tau, simplex_t sigma){
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
  traverse(st::k_simplices(this, root.get(), 1), [&idx_of, &uf](node_ptr cn, idx_t d){
    uf.Union(idx_of(cn->label), idx_of(cn->parent->label));
		return true; 
  });
  
  // Create the connected components
  std::transform(begin(v), end(v), begin(v), idx_of);
  return(uf.FindAll(v));
}

// Edge contraction 
inline void SimplexTree::contract(vector< idx_t > edge){
  set< node_ptr > to_remove;
  vector< simplex_t > to_insert; 
  traverse(st::preorder< true >(this, root.get()), [this, edge, &to_remove, &to_insert](node_ptr np, idx_t depth, simplex_t sigma){
    const idx_t va = edge[0], vb = edge[1];
    if (np->label == vb){ // only consider simplices which contain v_lb
      bool includes_a = std::find(sigma.begin(), sigma.end(), va) != sigma.end();
      if (includes_a){ // case 1: sigma includes both v_la and v_lb
        to_remove.insert(np); // add whole simplex to remove list, and we're done.
      } else { // case 2: sigma includes v_lb, but not v_la
        // Insert new simplices with v_la --> v_lb, identity otherwise
        const auto local_preorder = st::preorder< true >(this, np);
        traverse(local_preorder, [&to_insert, va, vb](node_ptr end, idx_t depth, simplex_t tau){
          std::replace(tau.begin(), tau.end(), vb, va);
          to_insert.push_back(tau); // si will be sorted upon insertion
          return true; 
        });
        to_remove.insert(np);
      }
    }
    return true; 
  });
  // Remove the simplices containing vb
	for (auto& edge: to_remove){ remove_subtree(edge); }
	for (auto& edge: to_insert){ insert_simplex(edge); }
}

// Recursive version output 
template < size_t I, typename OutputIt >
inline void SimplexTree::full_simplex_r(node_ptr cn, OutputIt out) const noexcept {
	if (cn == nullptr || cn->parent == nullptr){ return; };
	if constexpr(I > 0 && I <= array_threshold){
		auto dispatcher = make_index_dispatcher< I >();
		dispatcher([&out, &cn](auto idx) { *out = node_label_r< (I - idx - 1) >(cn); ++out; });
	} else {
	  splex_alloc_t a; 
	  splex_t buffer{a};
		auto buffer_it = std::back_inserter(buffer);
		while (cn->parent != nullptr){
			*(buffer_it)++ = cn->label;
			cn = cn->parent;
		}
		std::move(buffer.rbegin(), buffer.rend(), out);
	}
	return;
}

template< typename OutputIt >
inline void SimplexTree::full_simplex_out(node_ptr cn, const idx_t depth, OutputIt out) const noexcept {
	switch(depth){ 
		case 1: 
			return(full_simplex_r< 1 >(cn, out));
			break; 
		case 2: 
			return(full_simplex_r< 2 >(cn, out));
			break; 
		case 3: 
			return(full_simplex_r< 3 >(cn, out));
			break; 
		case 4: 
			return(full_simplex_r< 4 >(cn, out));
			break; 
		case 5: 
			return(full_simplex_r< 5 >(cn, out));
			break; 
		case 6: 
			return(full_simplex_r< 6 >(cn, out));
			break; 
		case 7: 
			return(full_simplex_r< 7 >(cn, out));
			break; 
		case 8: 
			return(full_simplex_r< 8 >(cn, out));
			break; 
		case 9: 
			return(full_simplex_r< 9 >(cn, out));
			break; 
		default:
			return(full_simplex_r< 0 >(cn, out));
			break;
	}
}

inline simplex_t SimplexTree::full_simplex(node_ptr cn, const idx_t depth) const noexcept {
  simplex_t result;
  result.reserve(depth);
  full_simplex_out(cn, depth, std::back_inserter(result));
  return(result);
}

// // Serialize the simplex tree
// inline vector< vector< idx_t > > SimplexTree::serialize() const{
//   using simplex_t = vector< idx_t >;
//   vector< simplex_t > minimal;
// 	traverse(st::maximal< true >(this, root.get()), [&minimal](node_ptr cn, idx_t depth, simplex_t sigma){
// 		minimal.push_back(sigma);
// 		return true;
// 	});
//   return(minimal);
// }

// Deserialization 
// inline void SimplexTree::deserialize(vector< vector< idx_t > > simplices){
//   for (auto& sigma: simplices){ insert_simplex(sigma); }
// }

template <typename Lambda>
inline void SimplexTree::traverse_facets(node_ptr s, Lambda f) const{

  // Constants
  const simplex_t sigma = full_simplex(s);
  const size_t sigma_depth = sigma.size();
  
  // Conditions for recursion
  const auto valid_eval = [sigma_depth](auto& cn) -> bool { return get< 1 >(cn) <= (sigma_depth - 1); };
  const auto valid_children = [&sigma, sigma_depth](auto& cn){ 
    return (get< 1 >(cn) <= (sigma_depth - 2)) && get< 0 >(cn)->label >= sigma.at(get< 1 >(cn)-1);
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
    
    while (cn != root.get()){
      // Check self 
      tau[c_depth-1] = cn->label;
      if (c_depth == sigma_depth-1 && std::includes(begin(sigma), end(sigma), begin(tau), begin(tau)+c_depth)){
        f(cn, c_depth); 
      }
      
      // Look at the siblings 
			using ut = decltype(*begin(node_children(root)));
			const auto compare_node_id = [&cn](ut& on) -> bool {
				return(cn->label == on->label);
			};
      auto tn = std::find_if(begin(cn->parent->children), end(cn->parent->children), compare_node_id);
      std::advance(tn,1);
      
      // Check siblings + their children up to facet depth 
      for (; tn != end(cn->parent->children); ++tn){ // TODO: end loop if tn != face
        simplex_t c_tau = full_simplex((*tn).get());
        c_tau.resize(sigma_depth);
				auto tr = st::preorder(this, (*tn).get(), valid_eval, valid_children);
        traverse(tr, [&c_tau, &sigma, &f](node_ptr t, idx_t d){
          c_tau[d-1] = t->label;
          if ((d == c_tau.size()-1) && std::includes(begin(sigma), end(sigma), begin(c_tau), begin(c_tau)+d)){
            f(t, d); 
          }
					return true; 
        });
      }
      
      // Move up
      cn = cn->parent; 
      c_depth--;
    }
  }
}

inline vector< idx_t > SimplexTree::get_vertices() const{
  if (n_simplexes.size() == 0){ return vector< idx_t >(0); }
  vector< idx_t > v;
  v.reserve(n_simplexes[0]);
  for (auto& cn: node_children(root)){ v.push_back(node_label(cn)); }
  return v;
}

// Returns whether the graph is acyclic.
inline bool SimplexTree::is_tree() const{
	if (n_simplexes.size() == 0){ return false; }
  UnionFind ds = UnionFind(n_simplexes.at(0));

  // Traverse the 1-skeleton, unioning all edges. If any of them are part of the same CC, there is a cycle. 
  const vector< idx_t > v = get_vertices();
  const auto index_of = [&v](const idx_t vid) -> size_t{ return std::distance(begin(v), std::find(begin(v), end(v), vid)); };
  
  // Apply DFS w/ UnionFind. If a cycle is detected, no more recursive evaluations are performed. 
	bool has_cycle = false; 
	auto st_dfs = st::k_skeleton< true >(this, root.get(), 1);
	for (auto& cn: st_dfs){
		const auto si = get< 2 >(cn);
		idx_t i1 = index_of(si.at(0)), i2 = index_of(si.at(1)); 
		if (ds.Find(i1) == ds.Find(i2)){ 
			has_cycle = true; 
			break; 
		}
    ds.Union(i1, i2);
	}
  return !has_cycle;
}



#endif
