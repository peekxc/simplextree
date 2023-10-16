#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]  
#include "simplextree.h"

// Given list of integer vectors
// [[Rcpp::export]]
bool nfold_intersection(vector< vector< int > > x, const size_t n){
  using it_t = vector< int >::iterator;
  auto ranges = vector< std::pair< it_t, it_t > >();
  std::transform(begin(x), end(x), std::back_inserter(ranges), [](vector< int >& cs){
    return std::make_pair(cs.begin(), cs.end());
  });
  bool is_connected = n_intersects(ranges, n);
  return(is_connected);
}

// Computes the nerve up to dimension k
// [[Rcpp::export]]
void nerve_expand(SEXP stx, vector< size_t > ids, vector< vector< int > > cover, const size_t k, const size_t threshold){
  const size_t n_sets = cover.size(); 
  if (ids.size() != n_sets){ stop("Invalid id/cover combination."); }
  
  // Extract the simplex tree 
  SimplexTree& st = *(Rcpp::XPtr< SimplexTree >(stx));  
  
  // Inserts vertices 
  std::array< idx_t, 1 > v; 
  for (auto v_id: ids){
    v[0] = v_id;
    st.insert_it(begin(v), end(v), st.root.get(), 0);
  }
  
  // Extract the range pairs
  using it_t = vector< int >::iterator; 
  using range_t = std::pair< it_t, it_t >;
  auto ranges = map< size_t, range_t >();
  size_t i = 0; 
  for (auto& c_set: cover){ ranges.emplace(ids[i++], std::make_pair(begin(c_set), end(c_set))); }

  // First insert all the edges w/ a common intersection
  using it_t2 = typename vector< size_t >::iterator; 
  for_each_combination(begin(ids), begin(ids)+2, end(ids), [&st, &ranges, threshold](it_t2 b, it_t2 e){
    auto edge = std::make_pair(*b, *std::next(b, 1));
    vector< range_t > sets = { ranges[edge.first], ranges[edge.second] };
    bool valid_edge = n_intersects(sets, threshold);
    if (valid_edge){
      st.insert_it(b, e, st.root.get(), 0);
    }
    return false; // always continue 
  });
  
  // Then perform the conditional k-expansion
  st.expansion_f(k, [&](node_ptr parent, idx_t depth, idx_t label){
    
    // Collect simplex to test
    auto k_simplex = st.full_simplex(parent, depth);
    k_simplex.push_back(label);
    
    // Extract current set of ranges
    auto current_sets = vector< std::pair< it_t, it_t > >();
    for (auto c_label: k_simplex){
      auto it = ranges.find(c_label);
      if (it != ranges.end()){ current_sets.push_back(it->second); }
    }
    
    // Test that their intersection is at least 'threshold'; if so, insert 
    if ((current_sets.size() == k_simplex.size()) && n_intersects(current_sets, threshold)){
      std::array< idx_t, 1 > int_label = { label };
      st.insert_it(begin(int_label), end(int_label), parent, depth);
    }
  });
  
  return; // Return nothing
}

// [[Rcpp::export]]
void nerve_expand_f(SEXP stx, vector< size_t > ids, Function include_f, const size_t k){
  
  // Extract the simplex tree 
  SimplexTree& st = *(Rcpp::XPtr< SimplexTree >(stx));  
  
  // Inserts vertices 
  std::array< idx_t, 1 > v; 
  for (auto v_id: ids){
    v[0] = v_id;
    st.insert_it(begin(v), end(v), st.root.get(), 0);
  }

  // First insert all the edges w/ a common intersection
  using it_t = vector< size_t >::iterator;
  for_each_combination(begin(ids), begin(ids)+2, end(ids), [&st, &include_f](it_t b, it_t e){
    IntegerVector edge = IntegerVector(b, e);
    LogicalVector valid_check = include_f(edge);  // This is needed to coerce correctly to bool
    bool valid_edge = is_true(all(valid_check));
    if (valid_edge){ st.insert_it(b, e, st.root.get(), 0); }
    return false; // always continue
  });

  // Then perform the conditional k-expansion
  st.expansion_f(k, [&](node_ptr parent, idx_t depth, idx_t label){

    // Collect simplex to test
    auto k_simplex = st.full_simplex(parent, depth);
    k_simplex.push_back(label);

    LogicalVector valid_check = include_f(k_simplex);  // This is needed to coerce correctly to bool
    bool valid = is_true(all(valid_check));

    // Test that their intersection is at least 'threshold'; if so, insert
    if (valid){
      std::array< idx_t, 1 > int_label = { label };
      st.insert_it(begin(int_label), end(int_label), parent, depth);
    }
  });
  
  return; // Return nothing
}

/*** R
library(simplextree)

# st <- simplex_tree(combn(4,2))
# simplextree:::nerve_expand(st$as_XPtr(), 2)

# set.seed(1234)
# st <- simplex_tree(as.list(seq(3)))
# simplextree:::nerve_expand(st$as_XPtr(), ids = st$vertices, cover = list(c(1,2,3), c(3, 4, 5, 6), c(3, 6, 7)), k = 2, threshold = 1)
# 
# 
# set.seed(1234)
# alphabet <- seq(50)
# cover <- lapply(seq(15), function(i){ 
#   set_size <- as.integer(runif(n = 1, min = 1, max = 15))
#   sample(alphabet, size = set_size, replace = FALSE) 
# })
# st <- simplex_tree(as.list(seq(length(cover))))
# simplextree:::nerve_expand(st$as_XPtr(), ids = st$vertices, cover = cover, k = 3, threshold = 1)


# simplextree:::nerve_comb(st$as_XPtr(), ids = seq(3), cover = list(c(1,2,3), c(3, 4, 5, 6), c(3, 6, 7)), k = 10, threshold = 1)
st
# 
# st <- simplex_tree()
# simplextree:::nerve_comb(st$as_XPtr(), list(c(1), c(3, 4, 5, 6)), seq(2), k = 5, n = 1)
# st
# 
# st <- simplex_tree()
# simplextree:::nerve_comb(st$as_XPtr(), list(c(1, 2), c(3, 4, 5, 6, 1), c(3, 1, 2)), seq(3), k = 5, n = 1)
# st
# 
# st <- simplex_tree()
# simplextree:::nerve_comb(st$as_XPtr(), list(c(1), c(3, 4, 5, 6, 1), c(3, 1, 2)), seq(3), k = 2, n = 1)
# st
# 
# set.seed(1234)
# cover <- lapply(seq(15), function(i){ sample(seq(as.integer(runif(1, min = 1, max = 30)))) })
# st <- simplex_tree()
# simplextree:::nerve_comb(st$as_XPtr(), cover, seq(length(cover)), k = 1, threshold = 1)
# 
# 
# all(combn(length(cover), 3, function(idx){ length(Reduce(function(x, y){ intersect(x, y) }, x = cover[idx])) > 1 }))
# 
# combn(length(cover), 2, function(idx){
#   simplextree:::nfold_intersection(cover[idx], 2) 
# })
# 
# 
# 
# simplextree:::nerve_comb(st$as_XPtr(), cover, seq(length(cover)), k = 1, n = 1)
# 
# 
# simplextree:::nfold_intersection(list(c(1,2,3), c(3, 4, 5, 6)), 1)

## Test if there is at least a single element in common
# nerve_cpp(st$as_XPtr(), list(c(1,2,3), c(3, 4, 5, 6)), 1)
# nerve_cpp(st$as_XPtr(), list(sample(seq(100)), sample(seq(150)), sample(seq(200))), 101)
# nerve_cpp(st$as_XPtr(), list(c(1,2,3), c(3, 4, 5, 6)), 1)
# nerve_cpp(st$as_XPtr(), list(c(2,3,1), c(5, 3, 4, 6), c(4, 6, 8, 3)), 1)
*/
