#ifndef SIMPLEXTREE_HPP_
#define SIMPLEXTREE_HPP_

#include "simplextree.h"


// --------- Begin C++ only API --------- 
// These functions are only available through the included header, and thus can only be accessed 
// on the C++ side. R-facing functions are exported through the module. 
using node_ptr = SimplexTree::node_ptr;
using simplex_t = SimplexTree::simplex_t; 

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
    for (idx_t cc = 0; cc < max; ++cc){
      if (find_by_id(root->children, cc) == nullptr){
        new_ids.push_back(cc);
        if (new_ids.size() == n){ break; }
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
    // auto it = level_map.find(std::to_string(vid) + "-2"); // Labels with id > v 
    auto it = level_map.find(encode_node(vid, 2));
    if (it != level_map.end()){ 
			const auto& cousins = (*it).second;
			for (const auto& ch: cousins){
				res_deg += node_children(ch).size(); 
			}
		}
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
      auto it = level_map.find(encode_node(id, 2));
      if (it != level_map.end()){ 
				const auto& cousins = (*it).second;
				for (const auto& ch: cousins){
					res_deg += node_children(ch).size(); 
				}
			}
      res_deg += node_children(cn).size(); // Labels with id < v 
      res.push_back(res_deg);
    }
  }
  return(res);
}

// Search the level map (cousins) to quickly get the adjacency relations. 
// The set of adjacency relations are the 0-simplexes connected to a given vertex v. 
inline vector< idx_t > SimplexTree::adjacent_vertices(const size_t v) const {
  
  // Resulting vector to return
  vector< idx_t > res = vector< idx_t >(); 
  
  // First extract the vertices which labels > v by checking edges
  //std::string key = std::to_string(v) + "-2";
  auto it = level_map.find(encode_node(v, 2));
  if (it != level_map.end()){
		const auto& cousins = (*it).second;
    for (const auto& cn: cousins){ 
			res.push_back(node_label(cn->parent)); 
		}
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
  // node_ptr child = find_by_id(parent->children, child_label);
	const auto eq_node_id_lambda = [child_label](auto& cn){
		return(child_label == node_label(cn));
	};
  auto child_it = std::find_if(begin(parent->children), end(parent->children), eq_node_id_lambda);
	if (*child_it){ // does not equal nullptr
    // Remove from level map 
    //std::string key = std::to_string(child_label) + "-" + std::to_string(child_depth);
    const size_t key = encode_node(child_label, child_depth);
    auto& cousins = level_map[key];
		auto child = (*child_it).get();
    cousins.erase(std::remove(begin(cousins), end(cousins), child), end(cousins));
    
    // Remove from parents children
		parent->children.erase(child_it);
    // parent->children.erase(std::remove(begin(parent->children), end(parent->children), child_it));
    record_new_simplexes(child_depth-1, -1);
  }
}

// Removes an entire subtree rooted as 'sroot', including 'sroot' itself. This function calls SimplexTree::remove_leaf(node_ptr parent, idx_t child_label) recursively. 
inline void SimplexTree::remove_subtree(node_ptr sroot){
  if (sroot == nullptr){ return; }
  if (sroot->children.empty()){ remove_leaf(sroot->parent, sroot->label); }  // TODO: figure out how to overload constructor
  else {
    // Remark: make sure to use labels instead of iterator here, otherwise the iterator will be invalidated.
    for (auto& cn: node_children(sroot)){
			remove_subtree(find_by_id(node_children(sroot), node_label(cn))); 
		}
		// std::for_each(begin(sroot->children), end(sroot->children), [this, &sroot](node_ptr cn){
		// 	remove_subtree(find_by_id(sroot->children, cn->label)); 
    // });
    if (sroot != root.get()){ remove_leaf(sroot->parent, sroot->label); }
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
  auto co = st::cofaces(this, cn);
  vector< node_ptr > co_v;
  std::transform(co.begin(), co.end(), std::back_inserter(co_v), [](auto& cn){ return(get< 0 >(cn)); });
  for (auto& cn: co_v){
    remove_subtree(cn);
  }
//   if (cn != nullptr && cn != root.get()){
// 		for (auto& coface: locate_cofaces(cn)){ 
// 			remove_subtree(coface); 
// 		}
// 	  // std::for_each(cofaces.begin(), cofaces.end(), [this](node_ptr coface){
//     //   remove_subtree(coface); // cofaces represent a set of subtrees
//     // });
//   }
};

// Inserts a child directly to the parent if the child doesn't already exist.
// Also records the addition of the simplex if the child is indeed created. Otherwise, does nothing. 
// inline SimplexTree::node_ptr SimplexTree::insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth){
//   if (new_child == nullptr){ return nullptr; }
//   auto& pc = c_parent->children; // parents children
//   node_ptr cn = find_by_id(pc, new_child->label);
//   if (cn == nullptr) { // child doesn't exist! create it
//     pc.insert(new_child);
//     record_new_simplexes(depth, 1);
//     return(new_child);
//   }
//   return(cn); 
// }

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
  insert((idx_t*) &sigma[0], 0, sigma.size(), root.get(), 0); // start the recursion from the root
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
        level_map[encode_node(labels[j], child_depth)].push_back((*new_it).get());
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
	const auto eq_node_id_lambda = [label](auto& cn){
		return(label == node_label(cn));
	};
  auto it = std::find_if(begin(level), end(level), eq_node_id_lambda);
  return it != end(level) ? (*it).get() : nullptr; 
}

// Wrapper to find a vertex from the top nodes
inline SimplexTree::node_ptr SimplexTree::find_vertex(const idx_t v_id) const{
  return find_by_id(root->children, v_id);
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
  size_t d = 1; 
  node_ptr r = cn; // copy 
  while (r && r->parent != root.get()){ 
    ++d; 
    r = r->parent;
		// std::cout << "1: " << r << ", " << r->parent << std::flush << std::endl;
		// if (r->parent == nullptr || r->parent == NULL){
		// 	std::cout << "here" << std::endl;
		// 	break; 
		// }
  }
	// std::cout << "got depth" << std::flush << std::endl;
  return d; 
}

// Utility to get the maximum height / longest path from any given node.
inline size_t SimplexTree::max_depth(node_ptr cn) const {
  // size_t max_height = depth(cn);
	auto dfs = st::preorder(this, cn);
	auto me = std::max_element(dfs.begin(), dfs.end(), [](auto& n1, auto& n2){
		return(get< 1 >(n1) < get< 1 >(n2));
	});
  return get< 1 >(*me);
}  

// Print the whole tree.
inline void SimplexTree::print_tree() const {
  print_subtree(root.get());
}

inline void SimplexTree::print_cousins() const {
	idx_t c_depth = 1;
	auto labels = get_vertices();
	while(c_depth <= tree_max_depth){
		for (auto &label: labels){
			auto ni = level_map.find(encode_node(label, c_depth)); 
			if (ni != level_map.end()){
				std::cout << "(last=" << label << ", depth=" << c_depth << "): ";
				for (auto& cousin_np: (*ni).second){
					// std::cout << cousin_np->label << std::flush << std::endl; 
					print_simplex(cousin_np, false);
					std::cout << " ";
				}
				std::cout << std::endl;
			}
		}
		++c_depth;
	}
};

// Basic breadth-first printing. Each level is prefixed with '.' <level> number of times, followed by the 
// the ids of the nodes at that breadth-level enclosed within parenthesis, e.g. ..( 2 3 4 ) 
inline void SimplexTree::print_subtree(node_ptr cn) const {
  for (const auto& ch: cn->children){
    idx_t h = max_depth(ch.get())-1; 
    std::cout << ch->label << " (h = " << h << "): ";
    for (int i = 1; i <= h; ++i){ 
      for (int j = 1; j <= i; ++j){ std::cout << "."; }
      std::cout << "("; print_level(ch.get(), i); std::cout << " )";
    }
    std::cout << std::endl;
  }
}

// Prints a given level of the tree
inline void SimplexTree::print_level(node_ptr cn, idx_t level) const{
  if (cn == nullptr || cn->parent == nullptr) return;
  if (level == 0) { printf(" %lu", cn->label); } 
  else if (level > 0 && (!cn->children.empty())) {
    for (const auto& ch: cn->children){
      print_level(ch.get(), level-1);
    }
  }
}

// Prints an individual simplex
inline void SimplexTree::print_simplex(node_ptr cn, bool newline) const {
  simplex_t si = full_simplex(cn);
  std::cout << "{ ";
  std::for_each(si.begin(), si.end(), [](const idx_t i){ std::cout << i << " "; });
  std::cout << "}";
	if (newline){ std::cout << std::endl; }
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

// Experimental k-expansion algorithm.
// Performs an expansion of order k, reconstructing the k-skeleton flag complex via an in-depth expansion of the 1-skeleton.
inline void SimplexTree::expansion(const idx_t k){
  for (auto& cn: node_children(root)){ 
    if (!node_children(cn).empty()){ 
			expand(cn->children, k-1); 
		}
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
  for (auto& cn: c_set){
    node_ptr top_v = find_vertex(cn->label);
    if (top_v != nullptr && (!top_v->children.empty())){
      
			// Temporary 
			vector< node_ptr > sib_ptrs;
			std::transform(siblings, end(c_set), std::back_inserter(sib_ptrs), [](const auto& n){
				return n.get();
			});

      // Get the intersection
      intersection.clear();
      std::set_intersection(
        begin(sib_ptrs), end(sib_ptrs),
				// siblings, end(c_set), 
        begin(top_v->children), end(top_v->children), 
        std::back_inserter(intersection), 
        [](auto& sib_n, auto& child_n) -> bool {
          return node_label(sib_n) < node_label(child_n); 
        }
      );
    
      // Insert and recursively expand 
      if (intersection.size() > 0){
        vector< idx_t > sigma = full_simplex(cn.get());
        sigma.resize(sigma.size() + 1);
        for (auto& int_node: intersection){
          sigma[sigma.size() - 1] = int_node->label;
          // insert_simplex(sigma); // slowest version w/ sorting + checking 
          // insert((idx_t*) &sigma[0], 0, sigma.size(), root, 0); // faster version w/o sorting
        	insert((idx_t*) &sigma[0], sigma.size()-1, sigma.size(), cn.get(), sigma.size() - 1); // fastest version
        }
        expand(cn->children, k-1); // recurse
      }
    }
    if (siblings != end(c_set)){ ++siblings; }
  }
}

// Recursive helper to extract the full simplex of a given node
// inline void SimplexTree::full_simplex_r(node_ptr node, vector< idx_t >& res) const{
//   if (node == nullptr || node == root.get()){ return; }
//   else {
//     res.push_back(node->label);
//     if (node->parent == root.get()){ return; }
//     else {
//       full_simplex_r(node->parent, res);
//     }
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

inline bool SimplexTree::vertex_collapseR(idx_t v1, idx_t v2, idx_t v3){
  node_ptr vp1 = find_vertex(v1), vp2 = find_vertex(v2), vt = find_vertex(v3);
  return vertex_collapse(vp1, vp2, vt); // collapse the free pair (vp1, vp2) --> vt
}

template <typename T>
vector<T> get_unique(vector<T> v) {
  auto it = std::unique(begin(v), end(v)); 
  v.resize(std::distance(begin(v), it));
  return v;
}

// inline void SimplexTree::get_cousins() const{
//   Rprintf("< id >-< depth >: < number of cousins >\n");
//   for (auto kv: level_map){
//     Rprintf("%d: %d\n", kv.first, kv.second.size()); 
//   }
// }

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

// R-version: TODO move
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
  traverse(st::max_skeleton(this, root.get(), 1), [&idx_of, &uf](node_ptr cn, idx_t d){
    uf.Union(idx_of(cn->label), idx_of(cn->parent->label));
		return true; 
  });
  
  // Create the connected components
  std::transform(begin(v), end(v), begin(v), idx_of);
  return(uf.FindAll(v));
}

// Alternative way to specify collapse
// TODO: Add vertex only collapse between vertices not necessarily connected by an edge
// void SimplexTree::collapse(node_ptr u, node_ptr v, node_ptr w){
//   
// }

// Edge contraction 
// TODO: fix
inline void SimplexTree::contract(vector< idx_t > edge){
  set< node_ptr > to_remove;
	// traverse_dfs(root.get(), [this, &edge, &to_remove](node_ptr sigma){
  //   const idx_t va = edge[0], vb = edge[1];
  //   if (sigma->label == vb){ // only consider simplices which contain v_lb
  //     vector< idx_t > si = full_simplex(sigma);
  //     bool includes_a = std::find(si.begin(), si.end(), va) != si.end();
  //     if (includes_a){ // case 1: sigma includes both v_la and v_lb
  //       to_remove.insert(sigma); // add whole simplex to remove list, and we're done.
  //     } else { // case 2: sigma includes v_lb, but not v_la
  //       // Insert new simplices with v_la --> v_lb, identity otherwise
  //       std::for_each(begin_dfs(sigma), end_dfs(), [this, va, vb](node_ptr end){
  //         vector< idx_t > to_insert = full_simplex(end);
  //         std::replace(to_insert.begin(), to_insert.end(), vb, va);
  //         insert_simplex(to_insert); // si will be sorted upon insertion
  //       });
  //       to_remove.insert(sigma);
  //     }
  //   }
  // });
  // std::for_each(begin_dfs(root.get()), end_dfs(), );
  // Remove the simplices containing vb
	for (auto& edge: to_remove){
		remove_subtree(edge);
	}
}

// template < size_t I >
// inline std::array< idx_t, I > SimplexTree::full_simplex(node_ptr node) const {
// //	if constexpr TODO: 
// 	// vector<idx_t> simplex = vector< idx_t >();
//   // full_simplex_r(node, simplex);
//   // std::reverse(simplex.begin(), simplex.end());
//   // return simplex;
// }

// Given a node, traces back up to the parent, forming the full simplex.
// inline vector<idx_t> SimplexTree::full_simplex(node_ptr node) const{
//   vector<idx_t> simplex = vector< idx_t >();
//   full_simplex_r(node, simplex);
//   std::reverse(simplex.begin(), simplex.end());
//   return simplex;
// }

// Recursive version
template < size_t I >
inline simplex_t SimplexTree::full_simplex_r(node_ptr cn) const noexcept {
	if (cn == nullptr || cn->parent == nullptr){ return(simplex_t()); };
	splex_alloc_t a; 
	splex_t buffer{a};
	if constexpr(I > 0 && I <= array_threshold){
		buffer.resize(I);
		auto dispatcher = make_index_dispatcher< I >();
		dispatcher([&buffer, &cn](auto idx) { buffer[idx] = node_label_r< (I - idx - 1) >(cn); });
	} else {
		auto out = std::back_inserter(buffer);
		while (cn->parent != nullptr){
			*(out)++ = cn->label;
			cn = cn->parent;
		}
		std::reverse(begin(buffer), end(buffer));
	}
	simplex_t result; 
	result.insert(end(result), std::make_move_iterator(begin(buffer)), std::make_move_iterator(end(buffer)));
	return result;
}

// Run-time version that supports compile-time optimizations
inline simplex_t SimplexTree::full_simplex(node_ptr cn, const idx_t depth) const noexcept {
	switch(depth){ 
		case 1: 
			return(full_simplex_r< 1 >(cn));
			break; 
		case 2: 
			return(full_simplex_r< 2 >(cn));
			break; 
		case 3: 
			return(full_simplex_r< 3 >(cn));
			break; 
		case 4: 
			return(full_simplex_r< 4 >(cn));
			break; 
		case 5: 
			return(full_simplex_r< 5 >(cn));
			break; 
		case 6: 
			return(full_simplex_r< 6 >(cn));
			break; 
		case 7: 
			return(full_simplex_r< 7 >(cn));
			break; 
		case 8: 
			return(full_simplex_r< 8 >(cn));
			break; 
		case 9: 
			return(full_simplex_r< 9 >(cn));
			break; 
		default:
			return(full_simplex_r< 0 >(cn));
			break;
	}
}

// Serialize the simplex tree
inline vector< vector< idx_t > > SimplexTree::serialize() const{
  using simplex_t = vector< idx_t >;
  vector< simplex_t > minimal;
	traverse(st::maximal< true >(this, root.get()), [&minimal](node_ptr cn, idx_t depth, simplex_t sigma){
		minimal.push_back(sigma);
		return true; 
	});
	// auto dfs_traversal = dfs< false >(this);
	// for (auto& node: dfs_traversal){
	// 	auto [sigma, d] = node; 
	// 	if (sigma != nullptr && sigma != root.get() && sigma->children.empty()){
  //     vector< node_ptr > coface_roots = locate_cofaces(sigma);
  //     if (coface_roots.size() == 1 && coface_roots.at(0) == sigma){
  //       minimal.push_back(full_simplex(sigma));
  //     }
  //   }
	// }
	// traverse_dfs(roo)
	// std::for_each(begin_dfs(root.get()), end_dfs(), [this, &minimal](node_ptr sigma){
	//   if (sigma != nullptr && sigma != root.get() && sigma->children.empty()){
	//     vector< node_ptr > coface_roots = locate_cofaces(sigma);
	//     if (coface_roots.size() == 1 && coface_roots.at(0) == sigma){
	//       minimal.push_back(full_simplex(sigma));
	//     }
	//   }
	// });
  return(minimal);
}

// Deserialization (C++ side) 
inline void SimplexTree::deserialize(vector< vector< idx_t > > simplices){
  using simplex_t = vector< idx_t >;
  std::for_each(begin(simplices), end(simplices), [this](const simplex_t& sigma){
    insert_simplex(sigma);
  });
}

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

// template <typename Lambda> 
// inline void SimplexTree::traverse_cofaces(node_ptr s, Lambda f) const{
//   if (s != nullptr){
//     for (auto& co: expand_subtrees( locate_cofaces(s) )){ 
// 			f(co); 
// 		}
//   }
// }

// template <typename Lambda> 
// inline void SimplexTree::traverse_link(node_ptr s, Lambda f) const{
//   if (s != nullptr){
//     for (auto& li: link(s)){ f(li); }
//   }
// }

// // Applies the lambda function 'f' to every simplex in the k-skeleton given by 'k'
// template <typename Lambda> 
// inline void SimplexTree::traverse_skeleton(node_ptr s, Lambda f, size_t k) const{
//   const auto valid_eval = [k](const node_ptr cn, const size_t d){ return d <= k + 1; };
//   const auto valid_children = [k](const node_ptr cn, const size_t d){ return d < k + 1; };
//   traverse_bfs_if(s, f, valid_eval, valid_children);
// }

// // Applies the lambda function 'f' to the maximal faces of the k-skeleton given by 'k'
// template <typename Lambda> 
// inline void SimplexTree::traverse_max_skeleton(node_ptr s, Lambda f, size_t k) const{
//   const auto valid_eval = [k](const node_ptr cn, const size_t d){ return d == k + 1; };
//   const auto valid_children = [k](const node_ptr cn, const size_t d){ return d < k + 1; };
//   traverse_bfs_if(s, f, valid_eval, valid_children);
// }

// Intermediate struct to enable faster filtration building
struct weighted_simplex {
  node_ptr np;
  double diam; 
  double face_diam;
  size_t dim; 
};

// Given a dimension k and set of weighted edges (u,v) representing the weights of the ordered edges in the trie, 
// constructs a std::function object which accepts as an argument some weight value 'epsilon' and 
// returns the simplex tree object.  
inline void SimplexTree::rips(vector< double > weights, const size_t k){
  // if (weights.size() != this->n_simplexes.at(1)){ stop("Must have one weight per edge."); return; }
  
  // 1. Perform k-expansion.
  expansion(k);

  // 2. Assign edge weights to a map; relies on precondition that 
  size_t i = 0;
  std::map< node_ptr, double > simplex_weights;
	traverse(st::max_skeleton(this, root.get(), 1), [&simplex_weights, &weights, &i](node_ptr cn, idx_t d){
    simplex_weights.emplace(cn, weights.at(i++)); 
		return true; 
  });

  // 2. Define way getting 'weight' of any simplex.
  std::function< double(node_ptr, idx_t) > sigma_weight; 
  sigma_weight = [this, &simplex_weights, &sigma_weight](node_ptr sigma, const idx_t d) -> double {
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
          traverse_facets(sigma, [&max_weight, &sigma_weight](node_ptr cn, idx_t d){
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
  traverse(st::k_skeleton(this, root.get(), k), [&w_simplices, &sigma_weight](node_ptr cn, size_t d){
    double w = sigma_weight(cn, d); 
    double pw = d == 0 ? -std::numeric_limits< double >::infinity() : sigma_weight(cn->parent, d-1);
    w_simplices.push_back(weighted_simplex{cn, w, pw, d});
		return true; 
  });
  
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
    if (sigma.np->parent == root.get() || sigma.np == root.get()){
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
// auto valid_eval = [&has_cycle](const node_ptr cn, const size_t d){ return (!has_cycle && d == 2); };
// auto valid_children = [&has_cycle](const node_ptr cn, const size_t d){ return (!has_cycle && d < 2); };
// traverse_dfs_if(root.get(), [this, &ds, &has_cycle, &index_of](node_ptr sigma, const size_t d){
//   vector< idx_t > si = full_simplex(sigma);
//   idx_t i1 = index_of(si.at(0)), i2 = index_of(si.at(1)); 
//   if (ds.Find(i1) == ds.Find(i2)){ has_cycle = true; }
//   ds.Union(i1, i2);
// }, valid_eval, valid_children);



// Link of a simplex
// inline vector< node_ptr > SimplexTree::link(node_ptr sigma) const{
//   vector< idx_t > s = full_simplex(sigma);
//   vector< node_ptr > links; 
// 	traverse(link< true >(this, sigma), [&links](const node_ptr cn, const idx_t depth, simplex_t tau){
// 		links.push_back(cn);
// 		return true; 
// 	});
// 	// auto dfs_traversal = dfs< true >(this); 
// 	// for (auto& tau: dfs_traversal){
// 	// 	auto [ np, d, t ] = tau; 
//   //   if (empty_intersection(t, s)){
//   //     vector< idx_t > potential_link;
//   //     std::set_union(s.begin(), s.end(), t.begin(), t.end(), std::back_inserter(potential_link)); 
//   //     node_ptr link_node = find_node(potential_link); 
//   //     if (link_node != nullptr){ links.insert(np); }
//   //   }
// 	// }
//   // vector< node_ptr > link_res(links.begin(), links.end()); 
//   return(links);
// }

// Applies condition DFS, evaluating the first predicate to determine whether to call 
// 'f' on a given element, and applying the second predicate to determine if the children
// of a given element should be added. 
// template <typename Lambda, typename P1, typename P2> 
// inline void SimplexTree::traverse_dfs_if(node_ptr s, Lambda f, P1 p1, P2 p2) const{
//   using d_node = std::pair< node_ptr, idx_t >; 
  
//   // Prepare to iteratively do DFS 
//   d_node current = std::make_pair(s, depth(s)); 
//   std::stack< d_node > node_stack; 
//   node_stack.push(current);
  
//   // Also track depth
//   while (!node_stack.empty()){
//     d_node cn = node_stack.top();
//     if (p1(cn.first, cn.second)){ 
//       f(cn.first, cn.second); // apply function
//     }
//     node_stack.pop();
//     if (p2(cn.first, cn.second) && !cn.first->children.empty()){
//       std::for_each(cn.first->children.rbegin(), cn.first->children.rend(), [&node_stack, &cn](auto& ch){
//         node_stack.push(std::make_pair(ch.get(), cn.second+1));
//       });
//     }
//   }
// }

// template <typename Lambda> 
// inline void SimplexTree::traverse_dfs(node_ptr s, Lambda f) const{
//   using d_node = std::pair< node_ptr, idx_t >; // node ptr + depth marker
  
//   // Prepare to iteratively do DFS 
//   d_node current = std::make_pair(s, depth(s)); 
//   std::stack< d_node > node_stack; 
//   node_stack.push(current);
  
//   // Also track depth
//   while (!node_stack.empty()){
//     current = node_stack.top();
//     f(current.first, current.second); // apply function
//     node_stack.pop();
//     if (!current.first->children.empty()){
//       std::for_each(current.first->children.rbegin(), current.first->children.rend(), [&node_stack, &current](auto& ch){
//         node_stack.push(std::make_pair(ch.get(), current.second+1));
//       });
//     }
//   }
// }

// Applies condition BFS, evaluating the first predicate to determine whether to call 
// 'f' on a given element, and applying the second predicate to determine if the children
// of a given element should be added. 
// template <typename Lambda, typename P1, typename P2> 
// inline void SimplexTree::traverse_bfs_if(node_ptr s, Lambda f, P1 p1, P2 p2) const{
//   using d_node = std::pair< node_ptr, idx_t >; // node ptr + depth marker
  
//   // Prepare to iteratively do BFS 
//   d_node current = std::make_pair(s, depth(s));
//   std::queue< d_node > node_queue; 
//   node_queue.push(current);
  
//   // Also track depth
//   while(!node_queue.empty()){
//     d_node cn = node_queue.front();
//     if (p2(cn.first, cn.second) &&  !cn.first->children.empty()){
//       for (auto nv = cn.first->children.begin(); nv != cn.first->children.end(); ++nv){ 
//         node_queue.emplace(std::make_pair((*nv).get(), cn.second+1)); 
//       }
//     }
    
//     if (p1(cn.first, cn.second)){
//       f(cn.first, cn.second);
//     }
//     node_queue.pop(); 
//   }
// }

// template < typename Lambda > 
// inline void SimplexTree::traverse_bfs(node_ptr s, Lambda f) const{
//   using d_node = std::pair< node_ptr, idx_t >; // node ptr + depth marker
  
//   // Prepare to iteratively do BFS 
//   d_node current = std::make_pair(s, depth(s));
//   std::queue< d_node > node_queue; 
//   node_queue.push(current);
  
//   // Also track depth
//   while(!node_queue.empty()){
//     current = node_queue.front();
//     if (current.first != nullptr && !current.first->children.empty()){
//       for (auto nv = current.first->children.begin(); nv != current.first->children.end(); ++nv){ 
//         node_queue.emplace(std::make_pair((*nv).get(), current.second+1)); 
//       }
//     }
//     f(current.first, current.second);
//     node_queue.pop(); 
//   }
// }

#endif
