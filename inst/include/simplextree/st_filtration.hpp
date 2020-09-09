#ifndef ST_FILTRATION_HPP_
#define ST_FILTRATION_HPP_

#include "simplextree.h"

// Intermediate struct to enable faster filtration building
struct weighted_simplex {
  node_ptr np;
  size_t depth; 
  double weight; 
};

// Indexed simplex 
struct indexed_simplex {
  int parent_idx;   // index of its parent simplex in tree
  idx_t label;      // last(sigma) 
  double value;     // diameter/weight of the simplex 
};

// Lexicographical comparison for simplices
struct s_lex_less {
	bool operator()(const simplex_t& s1, const simplex_t &s2) const {
		if (s1.size() == s2.size()){
			return std::lexicographical_compare(begin(s1), end(s1), begin(s2), end(s2));
		}
    return s1.size() < s2.size(); 
	}
}; 

// Weighted simplex lexicographically-refined comparison
struct ws_lex_less {
	SimplexTree* st; 
	explicit ws_lex_less(SimplexTree* st_) : st(st_){ }
	bool operator()(const weighted_simplex& s1, const weighted_simplex& s2) const {
		if (s1.weight == s2.weight){
			if (s1.depth == s2.depth){
				auto s1_simplex = st->full_simplex(s1.np, s1.depth);
				auto s2_simplex = st->full_simplex(s2.np, s2.depth);
				return s_lex_less()(s1_simplex, s2_simplex);	
			}
			return s1.depth < s2.depth;
		}
		return s1.weight < s2.weight;
	}
}; 

// Computes the maximum weight of a given simplex given 
// inline double max_weight_sparse(simplex_t sigma, const vector< double >& D, const size_t n) noexcept {
// 	switch(sigma.size()){
// 		case 0: case 1: { return(0.0); }
// 		case 2: { return(D[to_natural_2(sigma[0], sigma[1], n)]); }
// 		default: {
// 			double weight = 0.0; 
// 			for_each_combination(sigma.begin(), sigma.begin()+2, sigma.end(), [&D, &weight, n](auto& it1, auto& it2){
// 				const size_t idx = to_natural_2(*it1, *std::next(it1), n); 
// 				if (D[idx] > weight){ weight = D[idx]; }
// 				return false; 
// 			});
// 			return(weight);
// 		}
// 	};
// }

// Computes the maximum weight of a simplex 'sigma' by finding the highest 
// weighted edge representing a face of sigma in the distance vector 'D' 
// Note: sigma must constain 0-based contiguous indices here to use D 
// inline double max_weight_dense(simplex_t sigma, const vector< double >& D, const size_t n) noexcept {
// 	switch(sigma.size()){
// 		case 0: case 1: { return(0.0); }
// 		case 2: { return(D[to_natural_2(sigma[0], sigma[1], n)]); }
// 		default: {
// 			double weight = 0.0; 
// 			for_each_combination(sigma.begin(), sigma.begin()+2, sigma.end(), [&D, &weight, n](auto& it1, auto& it2){
// 				const size_t idx = to_natural_2(*it1, *std::next(it1), n); 
// 				if (D[idx] > weight){ weight = D[idx]; }
// 				return false; 
// 			});
// 			return(weight);
// 		}
// 	};
// }

// Use mixin to define a filtration 
class Filtration : public SimplexTree {
public: 
  // Filtration-specific fields 
  vector< bool > included;
  vector< indexed_simplex > fc; 
  
  // Constructor
  Filtration() : SimplexTree() {}
  
  // copy over the simplex tree
  void initialize(const SimplexTree& sc){
    auto max_tr = st::maximal< true >(&sc, sc.root.get());
    traverse(max_tr, [this](node_ptr cn, idx_t depth, simplex_t sigma){
      this->insert_it(begin(sigma), end(sigma), root.get(), 0);
      return true; 
    });
    id_policy = sc.id_policy;
  }
  
  // Filtration building methods
  void flag_filtration(const vector< double >&, const bool fixed=false);

  // Iterating through the filtration
  template < typename Lambda >
  void traverse_filtration(size_t, size_t, Lambda&&);
  
  // Changing the state of complex 
  void threshold_value(double);
  void threshold_index(size_t);
  
  // Querying information about the filtration
  vector< simplex_t > simplices() const;
  vector< double > weights() const;
  vector< idx_t > dimensions() const; 
  size_t current_index() const;
  double current_value() const;
  
  // Internal helpers
  vector< size_t > simplex_idx(const size_t) const;
  
  template< typename Lambda >
  void apply_idx(size_t, Lambda&&) const;
    
  // template < typename Iter > 
  // simplex_t expand_simplex(Iter, Iter) const;
  simplex_t expand_simplex(simplex_t) const; 
  
}; // End filtration 

// Given sorted vector 'ref', matches elements of 'x' with 'ref', returning
// a vector of the same length as 'x' with the matched indices. 
template < typename T >
inline vector< T > match(const vector< T >& x, const vector< T >& ref){
  vector< T > indices;
  indices.reserve(x.size());
  for(auto& elem: x){
    auto idx = std::distance(begin(ref), std::lower_bound(begin(ref), end(ref), elem));
    indices.push_back(idx);
  }
  return(indices);
}

struct sorted_edges {
  using it_t = vector< size_t >::iterator;
  vector< size_t > keys;
  const vector< double >& values;
  const vector< size_t > vertices; 

  sorted_edges(Filtration* st, const vector< double >& weights) : values(weights), vertices(st->get_vertices()) {
    const size_t n = vertices.size(); 
    auto edge_traversal = st::k_simplices< true >(dynamic_cast< SimplexTree* >(st), st->root.get(), 1);
    st::traverse(edge_traversal, [this, n](node_ptr np, idx_t depth, simplex_t edge){
      auto eid = match(edge, vertices);
      keys.push_back(to_natural_2(eid[0], eid[1], n)); 
      return true; 
    });
    if (!std::is_sorted(keys.begin(), keys.end())){ throw std::invalid_argument("keys not ordered."); }
  }
  
  // Given a simplex 'sigma' whose values are vertex ids in 'vertices', calculates the maximum weight
  double max_weight(simplex_t sigma){
    auto v_ids = match(sigma, vertices);
    double weight = 0.0; 
		for_each_combination(v_ids.begin(), v_ids.begin()+2, v_ids.end(), [this, &weight](it_t it1, it_t it2){
			const size_t idx = to_natural_2(*it1, *std::next(it1), vertices.size()); 
			const auto key_idx = std::lower_bound(keys.begin(), keys.end(), idx);
			const double ew = values[std::distance(keys.begin(), key_idx)];
			if (ew > weight){ weight = ew; }
			return false; 
		});
		return(weight);
  }
};


// Given a dimension k and set of weighted edges (u,v) representing the weights of the ordered edges in the trie, 
// constructs a std::function object which accepts as an argument some weight value 'epsilon' and 
// returns the simplex tree object.  
inline void Filtration::flag_filtration(const vector< double >& D, const bool fixed){
  if (this->tree_max_depth <= 1){ return; }
  if (D.size() != this->n_simplexes.at(1)){ throw std::invalid_argument("Must have one weight per edge."); }
  
  // 0. Create the sorted map between the edges and their weights
  sorted_edges se = sorted_edges(this, D);

  // 1. Calculate simplex weights from the edge weights 
  const size_t ns = std::accumulate(begin(this->n_simplexes), end(this->n_simplexes), 0, std::plus< size_t >());
  std::vector< weighted_simplex > w_simplices;
  w_simplices.reserve(ns);

	size_t i = 0; 
	traverse(st::level_order< true >(this), [&w_simplices, &D, &i, &se](node_ptr cn, idx_t d, simplex_t sigma){
    double c_weight = d == 1 ? 0.0 : (d == 2 ? D.at(i++) : se.max_weight(sigma));
    weighted_simplex ws = { cn, d, c_weight };
    w_simplices.push_back(ws); 
		return true; 
  });

  // 2. Sort simplices by weight to create the filtration
  std::sort(begin(w_simplices), end(w_simplices), ws_lex_less(this));
    
  // 3. Index the simplices
  fc.clear();
  fc.reserve(w_simplices.size());
  i = 0;  
  for (auto sigma_it = begin(w_simplices); sigma_it != end(w_simplices); ++sigma_it){
    auto& sigma = *sigma_it; 
    indexed_simplex tau;
    tau.label = sigma.np->label;
    tau.value = sigma.weight; 
    
    // Find the index of sigma's parent, or 0 if sigma if the parent is the empty face. 
    if (sigma.np->parent == this->root.get() || sigma.np == this->root.get()){
      tau.parent_idx = -1;
    } else {
      const auto sp = sigma.np->parent;
      auto p_it = std::find_if(begin(w_simplices), sigma_it, [&sp](const weighted_simplex& si){
        return(si.np == sp);
      });
      if (p_it == sigma_it){ throw std::range_error("sigma detected itself as the parent!"); }
      tau.parent_idx = std::distance(begin(w_simplices), p_it);
    }
    fc.push_back(tau);
  }
  
  // Set state of the filtration to the max
  included = vector< bool >(fc.size(), true);
}

// Accepts lambda that takes iterator start and end
template< typename Lambda >
inline void Filtration::apply_idx(size_t idx, Lambda&& f) const {
  if (idx >= fc.size()){ throw std::out_of_range("Bad simplex index"); }
  
  SmallVector< size_t >::allocator_type::arena_type arena;
  SmallVector< size_t > indices{ arena };
  indices.push_back(idx);
  while (fc[idx].parent_idx != -1){
    idx = fc[idx].parent_idx; 
    indices.push_back(idx);
  }
  
  std::for_each(indices.rbegin(), indices.rend(), [&f](size_t index){
    f(index);
  });
} 
 
// Returns the indices of where the labels that make up the simplex 
// at index 'idx' are in the filtration in ascending order.
inline vector< size_t > Filtration::simplex_idx(size_t idx) const {
  if (idx >= fc.size()){ throw std::out_of_range("Bad simplex index"); }
  vector< size_t > expanded = { idx };
  while (fc[idx].parent_idx != -1){
    idx = fc[idx].parent_idx; 
    expanded.push_back(idx);
  }
  std::reverse(begin(expanded), end(expanded));
  return(expanded);
}

// Given a vector of indices, expand the indices to form the simplex
// template < typename Iter > 
// inline vector< idx_t > Filtration::expand_simplex(Iter s, Iter e) const {
//   using index_t = typename std::iterator_traits< Iter >::value_type;
//   static_assert(std::is_integral< index_t >::value, "Integral type required for indexing.");
//   simplex_t sigma(std::distance(s, e));
//   std::transform(s, e, begin(sigma), [this](auto i){ std::cout << i << std::endl; return(fc.at(i).label); });
//   return(sigma);
// }

// inline void Filtration::get_simplex(const size_t idx, SmallVector< size_t >& out) const {
//   out.resize(0);
//   out.push_back(fc[idx].label);
//   size_t pidx = fc[idx].parent_idx;
//   while (pidx != -1){
//     out.push_back(fc[pidx].label);
//     pidx = fc[pidx].parent_idx; 
//   }
//   std::reverse(out.begin(), out.end());
// }
// template< typename Iter, typename OutputIt >
// inline void Filtration::expand_it(Iter s, Iter e, OutputIt out) const {
//   for(; s != e; ++s){
//     *out = fc[*s].label;
//   }
// }

inline simplex_t Filtration::expand_simplex(vector< size_t > indices) const {
  std::transform(begin(indices), end(indices), begin(indices), [this](size_t i){ return(fc.at(i).label); });
  return(indices);
}


template < typename Lambda >
inline void Filtration::traverse_filtration(size_t a, size_t b, Lambda&& f){
  if (b > fc.size()) { b = fc.size(); };
  if (a == b){ return; }
  
  SmallVector< size_t >::allocator_type::arena_type arena;
  SmallVector< size_t > expanded{ arena };
  expanded.reserve(tree_max_depth);

  const auto apply_f = [this, &expanded, &f](const size_t i){
    apply_idx(i, [this, &expanded](size_t index){ expanded.push_back(fc[index].label); });
    f(i, expanded.begin(), expanded.end());
    expanded.resize(0);
  };
  
  if (a < b){
    for (size_t i = a; i < b; ++i){ apply_f(i); }
    // for (size_t i = a; i < b; ++i){ f(i, expand_simplex(simplex_idx(i))); }
  }
  if (a > b){
    int i = a >= fc.size() ? fc.size() - 1 : a; // i possibly negative!
    //for (; i >= b && i >= 0; --i){ f(i, expand_simplex(simplex_idx(i))); }
    for (; i >= int(b) && i >= 0; --i){ apply_f(i); }
  }
  return;
}


// Get index corresponding to eps (inclusive)
inline void Filtration::threshold_value(double eps){
  using IS = indexed_simplex; 
  const auto eps_it = std::lower_bound(begin(fc), end(fc), eps, [](const IS s, double val) -> bool {
    return(s.value <= val); 
  });
  const size_t eps_idx = std::distance( begin(fc), eps_it );
  threshold_index(eps_idx);
}

inline void Filtration::threshold_index(size_t req_index){
  const size_t c_idx = current_index();
  const bool is_increasing = c_idx < req_index;
  using it_t = SmallVector< size_t >::iterator;
  traverse_filtration(c_idx, req_index, [this, is_increasing](const size_t i, it_t s, it_t e){
    included.at(i) = is_increasing; 
    is_increasing ? insert_it(s,e,root.get(),0) : remove(find_it(s,e,root.get()));
  });
}

// Returns the current index in the filtration
inline size_t Filtration::current_index() const {
  if (included.size() == 0){ return 0; }
  const size_t current_idx = std::distance(
    begin(included), std::find(begin(included), end(included), false)
  );
  return(current_idx);
}

inline double Filtration::current_value() const {
  const double INF = std::numeric_limits< double >::infinity(); 
  if (included.size() == 0){ return -INF; }
  const size_t c_idx = current_index();
  return(c_idx == fc.size() ? INF : fc[c_idx].value);
}

// Returns the simplices in the filtration in a list
inline vector< vector< idx_t > > Filtration::simplices() const {
  const size_t n = current_index();
  vector< vector< idx_t > > simplices(n);
  for (size_t i = 0; i < n; ++i){
    simplices[i] = expand_simplex(simplex_idx(i)); 
  }
  return simplices;
}

// Retrieves the filtration weights
inline vector< double > Filtration::weights() const{
  const size_t n = current_index();
  vector< double > weights = vector< double >(n);
  for (size_t i = 0; i < n; ++i){
    weights[i] = fc[i].value; 
  }
  return weights;
}

inline vector< idx_t > Filtration::dimensions() const{
  const size_t n = current_index();
  vector< idx_t > dims = vector< idx_t >(n);
  for (size_t i = 0; i < n; ++i){
    dims[i] = simplex_idx(i).size()-1;
  }
  return dims;
}

#endif 

