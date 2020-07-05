#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]  
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
// vector< size_t > nerve_comb(SEXP st, vector< vector< int > > cover, vector< size_t > ids, const size_t k, const size_t threshold) {
//   if (k == 0){ return vector< size_t >(); }
//   if (cover.size() < ids.size()){ stop("Must supply at least one id for each cover element"); }
//   Rcpp::XPtr< SimplexTree > stp(st);
//   using it_t = vector< int >::iterator; 
//   
//   // Extract the range pairs
//   auto ranges = vector< std::pair< it_t, it_t > >();
//   std::transform(begin(cover), end(cover), std::back_inserter(ranges), [](auto cs){
//     return std::make_pair(begin(cs), end(cs));
//   });
//   
//   // Get number of sets 
//   const size_t n_sets = cover.size(); 
//   
//   // Lambda to handle k-simplex in the nerve
//   const auto do_something = [](auto k_simplex){
//     for (auto i: k_simplex){ std::cout << i << ","; }
//     std::cout << std::endl;
//   };
//   
//   // Setup 
//   auto proper_faces = vector< size_t >();
//   auto current_faces = vector< size_t >();
//   
//   // Start comparing combinations
//   bool continue_nerve = true; 
//   for (size_t ki = 1; ki <= k && continue_nerve; ++ki){
//     
//     const size_t d = ki+1; // length of ki-simplex to test
//     const size_t fd = d-1; // length of (ki-1)-face to test
//     
//     // complete the loop if the previous combinations yielded no intersections
//     continue_nerve = (d < n_sets) && ( ki == 1 || proper_faces.size() > 0);
//     if (continue_nerve){
//       std::cout << "checking n: " << n_sets << " k: " << d << " combs" << std::endl;
//       
//       for_each_combination_idx(n_sets, d, [&current_faces, &proper_faces, &ranges, d, fd, n_sets, threshold](auto indices){
// 
//         // Extract current subset of ranges
//         auto current_sets = vector< std::pair< it_t, it_t > >();
//         std::transform(begin(indices), end(indices), std::back_inserter(current_sets), [&ranges](size_t idx){ return ranges.at(idx); });
// 
//         // Check all previous faces included
//         bool faces_included = true;
//         if (d > 2){
//           //std::cout << "fd: " << fd << ",";
//           for_each_combination(begin(indices), begin(indices)+fd, end(indices), [&proper_faces, &faces_included, fd, n_sets](auto fb, auto fe){
//             to_natural(fb, fe, n_sets, fd, [&proper_faces, &faces_included](const size_t idx){
//               auto it = std::lower_bound(begin(proper_faces), end(proper_faces), idx);
//               faces_included &= (it != end(proper_faces)) && (*it == idx);
//             });
//             return !faces_included; // if any face doesn't exist, stop testing combinations
//           });
//         } 
//         
//         // If all faces included, test the actual intersection; if all connected, add to the vector
//         if (faces_included && n_intersects(current_sets, threshold)){
//           const size_t N = BinomialCoefficient(n_sets, d);
//           // std::cout << "N: " << N << "d: " << d << std::endl;
//           to_natural(begin(indices), end(indices), n_sets, d, [N, &indices, &current_faces](const size_t idx){
//             if (idx >= N){
//               std::cout << "bad: " << idx << std::endl;
//               for (auto el: indices){ std::cout << el << ","; }
//               std::cout << std::endl;
//             }
//             current_faces.push_back(idx);
//           });
//         }
//         
//       }); // for_each_combination_idx
//       
//       // Do something for each of the k-faces
//       const auto [min, max] = std::minmax_element(begin(current_faces), end(current_faces));
//       std::cout << *min << "," << *max << std::endl; 
//       // to_subscript(begin(current_faces), end(current_faces), n_sets, d, do_something);
//     } // if (continue_nerve)
//     
//     if (ki == k){
//       return(current_faces);
//     }
//     // Prepare for next k 
//     proper_faces = std::move(current_faces);
//     current_faces.clear();
//     
//   } // for (size_t ki = 1; ki <= k; ++ki)
//   return vector< size_t >();
// }

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
    bool valid_edge = include_f(edge);
    if (valid_edge){ st.insert_it(b, e, st.root.get(), 0); }
    return false; // always continue 
  });
  
  // Then perform the conditional k-expansion
  st.expansion_f(k, [&](node_ptr parent, idx_t depth, idx_t label){
    
    // Collect simplex to test
    auto k_simplex = st.full_simplex(parent, depth);
    k_simplex.push_back(label);
    
    bool valid = include_f(k_simplex);
    
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
