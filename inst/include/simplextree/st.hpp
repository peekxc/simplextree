#ifndef SIMPLEXTREE_HPP_
#define SIMPLEXTREE_HPP_

#include "simplextree.h"


// --------- Begin C++ only API --------- 
// These functions are only available through the included header, and thus can only be accessed 
// on the C++ side. R-facing functions are exported through the module. 
using node_ptr = SimplexTree::node_ptr;
using simplex_t = SimplexTree::simplex_t; 
using namespace st;

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
  root.reset(new node(-1, nullptr));
	level_map.clear(); 
  n_simplexes.fill(0);
  // n_simplexes.clear();
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
  if (cousins_exist(v, 2)){
    traverse_cousins(v, 2, [&res](node_ptr cousin){
      res.push_back(node_label(cousin->parent)); 
    });
  }
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
inline void SimplexTree::record_new_simplexes(const idx_t k, const int n){
  if (k >= 32){ std::invalid_argument("Invalid dimension to record."); }
  n_simplexes[k] += n; 
  auto first_zero = std::find(n_simplexes.begin(), n_simplexes.end(), 0);
  tree_max_depth = std::distance(n_simplexes.begin(), first_zero);
  // if (n_simplexes.size() < k+1){ n_simplexes.resize(k+1); }
  // n_simplexes.at(k) += n;
  // while(n_simplexes.back() == 0 && n_simplexes.size() > 0){ n_simplexes.resize(n_simplexes.size() - 1); }
  // tree_max_depth = n_simplexes.size();
}

// Remove a child node directly from the parent, if it exists
// This will check that the child is a leaf, and if so, then it will: 
// 1) Remove the child from the parents children map
// 2) Remove the child from the level map 
inline void SimplexTree::remove_leaf(node_ptr parent, idx_t child_label){
  if (parent == nullptr){ return; }
  const idx_t child_depth = depth(parent) + 1;
  auto child_it = std::find_if(begin(parent->children), end(parent->children), [child_label](const node_uptr& cn)->bool{ return(cn->label == child_label); });
	if (child_it != end(parent->children)){ 
    // Remove from level map 
    auto child = (*child_it).get(); // copy regular node_ptr
    remove_cousin(child, child_depth);
    
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
    std::transform(begin(node_children(sroot)), end(node_children(sroot)), begin(nc), [](const node_uptr& u_np){
      return u_np.get(); 
    });
    for (auto cn: nc){
			remove_subtree(find_by_id(node_children(sroot), node_label(cn))); 
		}
    // Remove self
    if (sroot && sroot != root.get()){ remove_leaf(sroot->parent, sroot->label); }
  }
}

// First removes all the cofaces of a given simplex, including the simplex itself.
inline void SimplexTree::remove(node_ptr cn){
  if (cn != nullptr && cn != root.get()){
    auto cr = st::coface_roots< false >(this, cn);
    SmallVector< node_ptr >::allocator_type::arena_type arena;
    SmallVector< node_ptr > co_v{ arena };
    std::transform(cr.begin(), cr.end(), std::back_inserter(co_v), [](std::tuple< node_ptr, idx_t >& cn){ return(get< 0 >(cn)); });
    for (auto co_n: co_v){
      remove_subtree(co_n);
    }
  }
};

template< typename Iter >
inline auto SimplexTree::append_node(Iter pos, node_ptr cn, idx_t label, size_t depth) -> node_set_t::iterator {
  auto new_it = cn->children.emplace_hint(pos, make_unique< node >(label, cn));
  add_cousin((*new_it).get(), depth);
  record_new_simplexes(depth-1, 1);
  return(new_it);
}

// Inserts man simplices of a fixed dimension
// Assumes Iterator to int types are correct type to avoid casts + sorted
// template< typename Iter, size_t d >
// inline void SimplexTree::insert_fast(Iter s, Iter e){
// 
//   if constexpr(d == 0){
//     while(s != e){
//       node_children(root).emplace(make_unique< node >(*s, root))
//       ++s;
//     }
//   } else if constexpr ( d == 1 ){
//     const auto label1 = *s;
//     const auto label2 = *(s+1);
//     // const auto val_end = std::upper_bound(s, e, label); // wrong [1, 2, 1, 3, 1, 4]...
//     auto v_it = node_children(root).find(label);
//     if (v_it == end(node_children(root))){ v_it = append_node(v_it, root, label1, 1); }
//     if (v_it != end(node_children(root))){
//       v_it->
//     }
//   }
//   if (it == end(c_node->children)){
//     auto new_it = c_node->children.emplace_hint(it, make_unique< node >(label, c_node));
//     if (child_depth > 1){ // keep track of nodes which share ids at the same depth
//       add_cousin((*new_it).get(), child_depth);
//     }
//     record_new_simplexes(child_depth-1, 1);
//   }
// }

// Create a set of (i)-simplexes as children of the current node, if they don't already exist
// depth == (depth of c_node)
template< bool lex_order, typename Iter >
inline void SimplexTree::insert_it(Iter s, Iter e, node_ptr c_node, const idx_t depth){
  if (s == e || c_node == nullptr){ return; }
  // using it_t = typename Iter::value_t; 
  
  const idx_t child_depth = depth+1;
  std::for_each(s, e, [this, &c_node, child_depth](idx_t label){
    // if constexpr (lex_order){
    //   auto new_it = c_node->children.emplace_hint(c_node->children.end(), make_unique< node >(label, c_node));
    //   if (child_depth > 1){ // keep track of nodes which share ids at the same depth
    //     // level_map[encode_node(label, child_depth)].push_back((*new_it).get());
    //     add_cousin((*new_it).get(), child_depth);
    //   }
    //   record_new_simplexes(child_depth-1, 1);
    // } else {
      auto it = std::find_if(begin(node_children(c_node)), end(node_children(c_node)), [label](const node_uptr& cn){
        return(cn->label == label);
      });
      if (it == end(c_node->children)){
        auto new_it = c_node->children.emplace_hint(it, make_unique< node >(label, c_node));
        if (child_depth > 1){ // keep track of nodes which share ids at the same depth
          add_cousin((*new_it).get(), child_depth);
        }
        record_new_simplexes(child_depth-1, 1);
      }
    // }
  });
  
  // Recurse on the subtrees of the current node
  idx_t j = 1;
  std::for_each(s, e, [&](idx_t label){
    insert_it(std::next(s, j), e, find_by_id(c_node->children, label), depth+1);
    ++j;
  });
}

// Wrapper to find a vertex from the top nodes
template< typename Iterable >
inline void SimplexTree::insert(Iterable v) {
  static_assert(st::detail::has_begin< Iterable >::value, "Must be iterable object.");
  auto b = v.begin(), e = v.end();
  std::sort(b, e);          // Demand sorted labels
  e = std::unique(b, e);    // Demand unique labels
  insert_it(b, e, root.get(), 0);
}

// Create a set of (i)-simplexes as children of the current node, if they don't already exist
// inline void SimplexTree::insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth){
//   if (i >= n_keys || labels == nullptr || c_node == nullptr){ return; } // base case + safety checks
//   idx_t child_depth = depth + 1; // depth refers to parent depth, child_depth to child depth
//   for (int j = i; j < n_keys; ++j){
// 		using ut = decltype(*begin(node_children(c_node)));
// 		const auto compare_node_id = [&labels, &j](ut& cn) -> bool {
// 			return(cn->label == labels[j]);
// 		};
// 		auto it = std::find_if(begin(node_children(c_node)), end(node_children(c_node)), compare_node_id);
//  		// auto it = std::find_if(begin(node_children(c_node)), end(node_children(c_node)), eq_node_id(labels[j]));
//     if (it == end(c_node->children)){ // doesn't exist yet  
//       auto new_it = c_node->children.emplace_hint(it, std::make_unique< node >(labels[j], c_node));
// 			record_new_simplexes(depth, 1);
//       if (child_depth > tree_max_depth){ tree_max_depth = child_depth; }
//       if (child_depth > 1){ // keep track of nodes which share ids at the same depth
//         add_cousin((*new_it).get(), child_depth);
//         // level_map[encode_node(labels[j], child_depth)].push_back((*new_it).get());
//       }
//     }
//   }
//   // Recurse on the subtrees of the current node
//   for (int j = i; j < n_keys; ++j){
//     node_ptr child_node = find_by_id(c_node->children, labels[j]);
//     insert(labels, j + 1, n_keys, child_node, child_depth);
//   }
// }

// Overloaded in the case where a single (1-length vector) label is given
inline SimplexTree::node_ptr SimplexTree::find_by_id(const node_set_t& level, idx_t label) const{
  auto it = std::lower_bound(begin(level), end(level), label, [](const node_uptr& np, const idx_t id){
    return np->label < id;
  });
  return (it != end(level) && (*it)->label == label) ? (*it).get() : nullptr; 
}

// Wrapper to find a vertex from the top nodes
template< typename Iterable >
inline SimplexTree::node_ptr SimplexTree::find(Iterable v) const {
  // static_assert(std::is_integral<>)
  auto b = v.begin(), e = v.end();
  std::sort(b, e);          // Demand sorted labels
  e = std::unique(b, e);    // Demand unique labels
  return find_it(b, e, root.get());
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

// Recursively calculate the depth of a node
inline size_t SimplexTree::depth(node_ptr cn) const {
  if (cn == nullptr || cn == root.get()){ return 0; }
  size_t d;
  for (d = 1; cn && cn->parent != root.get(); ++d){ 
    cn = cn->parent;
  }
  return d; 
}

// Utility to get the maximum height / longest path from any given node.
inline size_t SimplexTree::max_depth(node_ptr cn) const {
	auto dfs = st::preorder< false >(this, cn);
  idx_t max_d = 0; 
  traverse(dfs, [&max_d](node_ptr np, idx_t depth){
    if (depth > max_d){ max_d = depth; }
    return true; 
  });
  // traverse_node_pairs(dfs, [&max_d](std::pair< node_ptr, idx_t > np){
  //   if (np.second > max_d){ max_d = np.second; }
  // });
	return(max_d);
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
    for (size_t i = 1; i <= h; ++i){ 
      for (size_t j = 1; j <= i; ++j){ os << "."; }
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

// Performs an expansion of order k, reconstructing the k-skeleton flag complex via an in-depth expansion of the 1-skeleton.
inline void SimplexTree::expansion(const idx_t k){
  expansion_f(k, [this](node_ptr parent, idx_t depth, idx_t label){
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
    node_ptr top_v = find_by_id(root->children, cn->label);
    if (top_v != nullptr && (!top_v->children.empty())){
      
			// Temporary 
			SmallVector< node_ptr >::allocator_type::arena_type arena2;
			SmallVector< node_ptr > sib_ptrs { arena2 } ;
			std::transform(siblings, end(c_set), std::back_inserter(sib_ptrs), [](const node_uptr& n){
				return (node_ptr) n.get();
			});

      // Get the intersection
      intersection.clear();
      std::set_intersection(
        begin(sib_ptrs), end(sib_ptrs),
        begin(top_v->children), end(top_v->children), 
        std::back_inserter(intersection), 
        less_np_label()
      );
    
      // Insert and recursively expand 
      if (intersection.size() > 0){
        for (auto& int_node: intersection){
          auto face = find_by_id(cn->children, int_node->label);
          if (face == nullptr){
            f((node_ptr) cn.get(), depth, int_node->label);
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

// Given two simplices tau and sigma, checks to see if tau is a face of sigma
// Assumes tau and sigma are both sorted.
inline bool SimplexTree::is_face(vector< idx_t > tau, vector< idx_t > sigma) const {
  auto tau_np = find(tau);
  auto sigma_np = find(sigma);
  if (tau_np != nullptr && sigma_np != nullptr){
    return std::includes(sigma.begin(), sigma.end(), tau.begin(), tau.end());
  }
  return false; 
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
  node_ptr vp1 = find_by_id(root->children, v1); 
  node_ptr vp2 = find_by_id(root->children, v2); 
  node_ptr vt = find_by_id(root->children, v3);
  return vertex_collapse(vp1, vp2, vt); // collapse the free pair (vp1, vp2) --> vt
}

// Vertex collapse - A vertex collapse, in this sense, is the result of applying a 
// peicewise map f to all vertices sigma \in K, where given a pair (u,v) -> w, 
// f is defined as: 
// f(x) = { (1) w if x in { u, v }, (2) x o.w. }
inline bool SimplexTree::vertex_collapse(node_ptr vp1, node_ptr vp2, node_ptr vt){
	// Lambda to do the mapping
	vector< simplex_t > to_insert; 
	const auto map_collapse = [&to_insert, vt](simplex_t si, node_ptr vp){
		std::replace(begin(si), end(si), vp->label, vt->label); 
	  to_insert.push_back(si);
	};
	
	// Enumerates the cofaces of each node, performing the maps on each one
	for (auto& cn: cofaces< true >(this, vp1)){ map_collapse(get< 2 >(cn), vp1); }
	for (auto& cn: cofaces< true >(this, vp2)){ map_collapse(get< 2 >(cn), vp2); }
	
	// Insert all the mapped simplices, remove vertices if they exist
	for (auto& sigma: to_insert){ insert(sigma); }
	if (vp1 != vt) { remove(find_by_id(root->children, vp1->label)); }
  if (vp2 != vt) { remove(find_by_id(root->children, vp2->label)); }
  return true; 
}

// Elementary collapse - only capable of collapsing sigma through tau, and only if tau has sigma 
// as its only coface. There are technically two cases, either tau and sigma are both leaves or 
// tau contains sigma as its unique child. Both can be handled by removing sigma first, then tau.
inline bool SimplexTree::collapse(node_ptr tau, node_ptr sigma){
  if (tau == nullptr || sigma == nullptr){ return false; }
	auto tau_cofaces = st::cofaces< false >(this, tau);
	bool sigma_only_coface = true; 
	traverse(tau_cofaces, [&tau, &sigma, &sigma_only_coface](node_ptr coface, idx_t depth){
		sigma_only_coface &= (coface == tau) || (coface == sigma);
		return(sigma_only_coface); 
	});
	if (sigma_only_coface){
		remove_leaf(sigma->parent, sigma->label);
  	remove_leaf(tau->parent, tau->label);
		return(true);
	}
	return(false);
}

// Returns the connected components given by the simplicial complex
inline vector< idx_t > SimplexTree::connected_components() const{
  
  // Provide means of mapping vertex ids to index values
  vector< idx_t > v = get_vertices(); // vertices are ordered, so lower_bound is valid
  const auto idx_of = [&v](const idx_t val) { return(std::distance(begin(v), std::lower_bound(begin(v), end(v), val))); };
  
  // Traverse the edges, unioning vertices 
  UnionFind uf = UnionFind(root->children.size());
  traverse(st::k_simplices< false >(this, root.get(), 1), [&idx_of, &uf](node_ptr cn, idx_t d){
    uf.Union(idx_of(cn->label), idx_of(cn->parent->label));
		return true; 
  });
  
  // Create the connected components
  std::transform(begin(v), end(v), begin(v), idx_of);
  return(uf.FindAll(v));
}

// Edge contraction 
inline void SimplexTree::contract(vector< idx_t > edge){
  vector< simplex_t > to_remove;
  vector< simplex_t > to_insert; 
  traverse(st::preorder< true >(this, root.get()), [this, edge, &to_remove, &to_insert](node_ptr np, idx_t depth, simplex_t sigma){
    const idx_t va = edge[0], vb = edge[1];
    if (np->label == vb){ // only consider simplices which contain v_lb
      bool includes_a = std::find(sigma.begin(), sigma.end(), va) != sigma.end();
      if (includes_a){ // case 1: sigma includes both v_la and v_lb
        to_remove.push_back(sigma); // add whole simplex to remove list, and we're done.
      } else { // case 2: sigma includes v_lb, but not v_la
        // Insert new simplices with v_la --> v_lb, identity otherwise
        const auto local_preorder = st::preorder< true >(this, np);
        traverse(local_preorder, [&to_insert, va, vb](node_ptr end, idx_t depth, simplex_t tau){
          std::replace(tau.begin(), tau.end(), vb, va);
          to_insert.push_back(tau); // si will be sorted upon insertion
          return true; 
        });
        to_remove.push_back(sigma);
      }
    }
    return true; 
  });
  
  // for (auto& edge: to_remove){ print_simplex(std::cout, edge, true); }
  
  // Remove the simplices containing vb
	for (auto& edge: to_remove){ remove(find(edge)); }
	for (auto& edge: to_insert){ insert(edge); }
}

template < typename Lambda > // Assume lambda is boolean return 
void SimplexTree::traverse_up(node_ptr cn, const size_t depth, Lambda&& f) const noexcept {
  if (cn == nullptr || cn->parent == nullptr){ return; };
  switch(depth){
	  case 6:
			f(cn);
	    cn = cn->parent; 
		case 5:
		  f(cn);
	    cn = cn->parent; 
		case 4:
			f(cn);
	    cn = cn->parent; 
		case 3:
			f(cn);
	    cn = cn->parent; 
		case 2:
			f(cn);
	    cn = cn->parent; 
	  case 1:
			f(cn);
	    break; 
    default: 
      idx_t d = 0; 
      while (cn != root.get() && cn->parent != nullptr && d <= tree_max_depth){
    		f(cn);
    		cn = cn->parent;
    		++d;
    	}
      break;
  }
}
  
template< typename OutputIt >
inline void SimplexTree::full_simplex_out(node_ptr cn, const idx_t depth, OutputIt out) const noexcept {
   if (cn == nullptr || cn == root.get()){ return; }
   if (depth == 0){
    std::deque< idx_t > labels; 
    traverse_up(cn, depth, [&labels](node_ptr np){ labels.push_front(np->label); });
    std::move(begin(labels), end(labels), out);
  } else {
    splex_alloc_t a; 
    splex_t labels{a};
    labels.resize(depth);
    size_t i = 1; 
    traverse_up(cn, depth, [&depth, &i, &labels](node_ptr np){ labels.at(depth - (i++)) = np->label; });
    std::move(begin(labels), end(labels), out);
  }
}

inline simplex_t SimplexTree::full_simplex(node_ptr cn, const idx_t depth) const noexcept {
  simplex_t labels; 
  labels.reserve(depth);
  full_simplex_out(cn, depth, std::back_inserter(labels));
  return(labels);
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

// template <typename Lambda>
// inline void SimplexTree::traverse_facets(node_ptr s, Lambda f) const{
// 
//   // Constants
//   const simplex_t sigma = full_simplex(s);
//   const size_t sigma_depth = sigma.size();
//   
//   // Conditions for recursion
//   const auto valid_eval = [sigma_depth](auto& cn) -> bool { return get< 1 >(cn) <= (sigma_depth - 1); };
//   const auto valid_children = [&sigma, sigma_depth](auto& cn){ 
//     return (get< 1 >(cn) <= (sigma_depth - 2)) && get< 0 >(cn)->label >= sigma.at(get< 1 >(cn)-1);
//   };
//   
//   // Facet search 
//   if (sigma_depth <= 1){ return; }
//   else if (sigma_depth == 2){ 
//     // The facet of an edge is just its two vertices
//     f(find_by_id(root->children, s->label));
//     f(s->parent, 1);
//     return; 
//   } else {
//     size_t c_depth = sigma_depth - 1; 
//     node_ptr cn = s->parent; 
//     simplex_t tau = sigma;
//     
//     while (cn != root.get()){
//       // Check self 
//       tau[c_depth-1] = cn->label;
//       if (c_depth == sigma_depth-1 && std::includes(begin(sigma), end(sigma), begin(tau), begin(tau)+c_depth)){
//         f(cn, c_depth); 
//       }
//       
//       // Look at the siblings 
// 			using ut = decltype(*begin(node_children(root)));
//       auto tn = std::find_if(begin(cn->parent->children), end(cn->parent->children), eq_node_id(cn->label));
//       std::advance(tn,1);
//       
//       // Check siblings + their children up to facet depth 
//       for (; tn != end(cn->parent->children); ++tn){ // TODO: end loop if tn != face
//         simplex_t c_tau = full_simplex((*tn).get());
//         c_tau.resize(sigma_depth);
// 				auto tr = st::preorder(this, (*tn).get(), valid_eval, valid_children);
//         traverse(tr, [&c_tau, &sigma, &f](node_ptr t, idx_t d){
//           c_tau[d-1] = t->label;
//           if ((d == c_tau.size()-1) && std::includes(begin(sigma), end(sigma), begin(c_tau), begin(c_tau)+d)){
//             f(t, d); 
//           }
// 					return true; 
//         });
//       }
//       
//       // Move up
//       cn = cn->parent; 
//       c_depth--;
//     }
//   }
// }

inline vector< idx_t > SimplexTree::get_vertices() const{
  if (tree_max_depth == 0){ return vector< idx_t >(0); }
  vector< idx_t > v;
  v.reserve(n_simplexes[0]);
  for (auto& cn: node_children(root)){ v.push_back(node_label(cn)); }
  return v;
}

// Returns whether the graph is acyclic.
inline bool SimplexTree::is_tree() const{
	if (tree_max_depth == 0){ return false; }
  UnionFind ds = UnionFind(n_simplexes.at(0));

  // Traverse the 1-skeleton, unioning all edges. If any of them are part of the same CC, there is a cycle. 
  const vector< idx_t > v = get_vertices();
  const auto index_of = [&v](const idx_t vid) -> size_t{ return std::distance(begin(v), std::find(begin(v), end(v), vid)); };
  
  // Apply DFS w/ UnionFind. If a cycle is detected, no more recursive evaluations are performed. 
	bool has_cycle = false; 
	auto st_dfs = st::k_simplices< true >(this, root.get(), 1);
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
