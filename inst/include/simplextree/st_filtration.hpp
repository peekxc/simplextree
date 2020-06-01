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

// Computes the maximum weight of a simplex 'sigma' by finding the highest 
// weighted edge representing a face of sigma in the distance vector 'D' 
inline double max_weight(simplex_t sigma, const vector< double >& D, const size_t n) noexcept {
	switch(sigma.size()){
		case 0: case 1: { return(0.0); }
		case 2: { return(D[to_natural_2(sigma[0], sigma[1], n)]); }
		default: {
			double weight = 0.0; 
			for_each_combination(sigma.begin(), sigma.begin()+2, sigma.end(), [&D, &weight, n](auto& it1, auto& it2){
				const size_t idx = to_natural_2(*it1, *std::next(it1), n); 
				if (D[idx] > weight){ weight = D[idx]; }
				return false; 
			});
			return(weight);
		}
	};
}

// Use mixin to define a filtration 
class Filtration : public SimplexTree {
public: 
  // Filtration-specific fields 
  vector< bool > included;
  vector< indexed_simplex > fc; 
  
  // Constructor
  Filtration() : SimplexTree() {}
  
  // copy over the simplex tree
  void initialize(const SimplexTree& st){
    deserialize(st.serialize());
    id_policy = st.id_policy;
  }
  
  // Filtration building methods
  void flag(const vector< double >&, const bool fixed=false);

  // Changing the state of complex 
  void threshold_value(double);
  void threshold_index(size_t);
  
  // Querying information about the filtration
  vector< vector< idx_t > > simplices() const;
  vector< double > weights() const;
  size_t current_index() const;
  double current_value() const;
  
  // Internal helpers
  vector< size_t > simplex_idx(const size_t) const;
  
  // template < typename Iter > 
  // simplex_t expand_simplex(Iter, Iter) const;
  simplex_t expand_simplex(simplex_t) const; 
};

// use PIMPL to handle filtration 
// class Filtration {
//   private:
//   	SimplexTree* st { nullptr }; // back pointer
//   public: 
//     Filtration(SimplexTree* st_ptr) : SimplexTree(st_ptr) { }
//   	~Filtration() = default;
// };

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


// Given a dimension k and set of weighted edges (u,v) representing the weights of the ordered edges in the trie, 
// constructs a std::function object which accepts as an argument some weight value 'epsilon' and 
// returns the simplex tree object.  
inline void Filtration::flag(const vector< double >& D, const bool fixed){
  if (D.size() != this->n_simplexes.at(1)){ throw std::invalid_argument("Must have one weight per edge."); }
  const size_t n = this->n_simplexes.at(0);
  
  // 0. Get sorted vertices to match against
  const auto v = this->get_vertices();
  
  // 1. Calculate simplex weights from the distances 
  const size_t ns = std::accumulate(begin(this->n_simplexes), end(this->n_simplexes), 1, std::multiplies< size_t >());
  std::vector< weighted_simplex > w_simplices;
  w_simplices.reserve(ns);
	
	size_t i = 0; 
	traverse(st::level_order< true >(this), [&w_simplices, &D, &i, &v, n](node_ptr cn, idx_t d, simplex_t sigma){
    const double NEG_INF = -std::numeric_limits< double >::infinity(); 
	  auto sigma_idx = match(sigma, v);
    double c_weight = d == 1 ? NEG_INF : (d == 2 ? D.at(i++) : max_weight(sigma_idx, D, n));
    weighted_simplex ws = { cn, d, c_weight };
    // std::cout << "depth: " << d << ", weight: " << c_weight << std::endl; 
    w_simplices.push_back(ws); 
		return true; 
  });

  // 2. Sort simplices by weight to create the filtration
  std::sort(begin(w_simplices), end(w_simplices), ws_lex_less(this));
  
  size_t ii = 0; 
  for (auto& ws: w_simplices){
    std::cout << ii << ": label=" << ws.np->label << ", depth=" << ws.depth << ", weight=" << ws.weight << std::endl;
    ii++;
  }
    
  // 3. Index the simplices
  fc.clear();
  fc.reserve(w_simplices.size());
  i = 0;  
  for (auto sigma_it = begin(w_simplices); sigma_it != end(w_simplices); ++sigma_it){
    auto& sigma = *sigma_it; 
    indexed_simplex tau;
    tau.label = sigma.np->label;
    tau.value = sigma.weight; 
    // std::cout << "Iter: " << i++ << ", weight: " << tau.value << ", label: " << tau.label << std::endl; 
    // Find the index of sigma's parent, or 0 if sigma if the parent is the empty face. 
    if (sigma.np->parent == this->root.get() || sigma.np == this->root.get()){
      tau.parent_idx = -1;
    } else {
      // Search for the lower bound on where sigma's parents weight is
      // auto lb = std::lower_bound(begin(w_simplices), sigma_it, sigma.face_diam, 
      //   [](const weighted_simplex& si, const double& val) -> bool {
      //   return(si.diam < val);
      // });
      const auto sp = sigma.np->parent;
      
      auto p_it = std::find_if(begin(w_simplices), sigma_it, [&sp](const weighted_simplex& si){
        return(si.np == sp);
      });
      if (p_it == sigma_it){ throw std::range_error("sigma detected itself as the parent!"); }
      tau.parent_idx = std::distance(begin(w_simplices), p_it);
      // std::cout << "parent idx: " << *tau.parent_idx << std::endl; 
    }
    fc.push_back(tau);
  }
  
  size_t jj = 0; 
  for (auto& is: fc){
    std::cout << jj << ": label=" << is.label << ", value=" << is.value << ", parent index=" << is.parent_idx << std::endl;
    jj++;
  }
  
  
  // Set state of the filtration to the max
  included = vector< bool >(fc.size(), true);
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
  // std::cout << std::endl; 
  std::reverse(begin(expanded), end(expanded));
  // std::cout << idx << " expanded: "; for (auto el: expanded){ std::cout << el << ","; }; std::cout << std::flush << std::endl; 
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

inline simplex_t Filtration::expand_simplex(vector< size_t > indices) const {
  std::transform(begin(indices), end(indices), begin(indices), [this](auto i){ return(fc.at(i).label); });
  return(indices);
}

// Get index corresponding to eps (inclusive)
inline void Filtration::threshold_value(double eps){
  using IS = indexed_simplex; 
  const auto eps_it = std::lower_bound(begin(fc), end(fc), eps, [](const IS s, double val) -> bool {
    return(s.value <= val); 
  });
  const size_t eps_idx = std::distance( begin(fc), (eps_it-1) );
  threshold_index(eps_idx);
}

inline void Filtration::threshold_index(size_t req_index){
  if (req_index > fc.size()) { req_index = fc.size(); };
  
  // Get current index; if matches requested index, end
  auto current_idx = current_index();
  if (current_idx == req_index){ return; }
 
  // Insert simplices as needed
  if (current_idx < req_index){
    for (size_t i = current_idx; i < req_index; ++i){
      std::cout << i << "," << std::flush << std::endl; 
      included.at(i) = true; 
      insert_simplex(expand_simplex(simplex_idx(i)));
    }
    return; 
  } 
  
  // Or remove simplices as needed
  if (current_idx > req_index){
    int i = current_idx >= fc.size() ? fc.size() - 1 : current_idx; // possibly negative
    for (; i >= req_index && i > 0; --i){
      std::cout << i << "," << std::flush << std::endl; 
      included.at(i) = false; 
      remove_simplex(expand_simplex(simplex_idx(i)));
    }
    return;
  }
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
  const size_t n = fc.size();
  vector< vector< idx_t > > simplices(n);
  for (size_t i = 0; i < n; ++i){
    simplices[i] = expand_simplex(simplex_idx(i)); 
  }
  return simplices;
}

// Retrieves the filtration weights
inline vector< double > Filtration::weights() const{
  const size_t n = fc.size();
  vector< double > weights = vector< double >(n);
  for (size_t i = 0; i < n; ++i){
    weights[i] = fc[i].value; 
  }
  return weights;
}

#endif 