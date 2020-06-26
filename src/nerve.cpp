#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]  
#include "simplextree.h"

// Template specialization for vector output.
// From: http://ptspts.blogspot.com/2013/12/how-to-implement-in-place-set.html
// template< typename Iter, typename T >
// static void IntersectionUpdate(std::vector<T>* bc, Iter a, const Iter a_end) {
//   auto b = bc->begin();
//   auto b_high = bc->end();
//   while (a != a_end && b != b_high) {
//     if (*a < *b) { ++a; } 
//     else if (*a > *b) {
//       std::iter_swap(b, --b_high);  // Works even if swapping with itself.
//     } else {  // Elements are equal, keep them in the intersection.
//       ++a; ++b;
//     }
//   }
//   bc->erase(b, bc->end());  // Remove remaining elements in bc but not in ac.
// }

template < typename T >
bool intervals_disjoint(vector< pair< T, T > > intervals){
  if (intervals.size() <= 1){ return(true); }
  
  auto interval_ids = vector< std::pair< T, T > >();
  T i = 0; 
  for (auto& interval: intervals){
    interval_ids.push_back( std::make_pair(i, interval.first) );
    interval_ids.push_back( std::make_pair(i, interval.second) );
    ++i;
  }
  
  // Sort by the values of value types  
  std::stable_sort(begin(interval_ids), end(interval_ids), [](auto& p1, auto& p2){
    return p1.second < p2.second; 
  });

  // Search for adjacent indices of the same id
  auto adj = std::adjacent_find(begin(interval_ids), end(interval_ids), [](auto& p1, auto& p2){
    return p1.first == p2.first;
  });
  
  // If any adjacent elements were found, then there's at least one interval disjoint from the rest
  if (adj != end(interval_ids)){
    bool is_disjoint = true; 
    for (size_t i = 0; i < (interval_ids.size())-1; ++i){
      auto& p1 = interval_ids.at(i);
      auto& p2 = interval_ids.at(i+1);
      is_disjoint = is_disjoint && (p1.first == p2.first ? p1.second <= p2.second : p1.second < p2.second);
    }
    return(is_disjoint);
  } 
  return(false);
}

template <class Iter, class Incr>
auto safe_advance(const Iter& curr, const Iter& end, Incr n) {
  size_t remaining(std::distance(curr, end));
  if (remaining < n) { n = remaining; }
  return std::next(curr, n);
}

// Checks if there are at least n elements common to all given sorted ranges
template < typename Iter >
bool n_intersects_sorted(vector< pair< Iter, Iter > > ranges, const size_t n){
  using T = typename Iter::value_type; 
  if (n == 0){ return(true); }
  if (ranges.size() <= 1){ return(false); }
  
  // Sort by size, then fold a set_intersection
  std::sort(begin(ranges), end(ranges), [](auto& p1, auto& p2){
    return std::distance(p1.first, p1.second) < std::distance(p2.first, p2.second);
  });
 
  // Fold a sorted intersection 
  auto common = vector< T >();
  std::set_intersection(ranges[0].first, ranges[0].second, ranges[1].first, ranges[1].second, std::back_inserter(common));
  for (size_t i = 2; i < ranges.size(); ++i){
    auto aux = vector< T >();
    std::set_intersection(begin(common), end(common), ranges[i].first, ranges[i].second, std::back_inserter(aux));
    if (aux.size() < n){
      return(false);
    }
    common.resize(aux.size());
    std::move(begin(aux), end(aux), begin(common));
  }
  return(common.size() >= n);
}

// Checks if there are at least n elements common to all given unsorted ranges
template < typename Iter >
bool n_intersects_unsorted(vector< pair< Iter, Iter > > ranges, const size_t n, const size_t k=32){
  using T = typename Iter::value_type; 
  
  // Use lambda to determine when finished
  const auto finished = [&ranges](){
    size_t n_finished = std::accumulate(begin(ranges), end(ranges), (size_t) 0, [](size_t val, auto& p){ return p.first == p.second; });
    return(n_finished >= 2);
  };
  
  auto mv = std::set< T >();          // the values mutally common to all the ranges
  auto vc = std::map< T, size_t >();  // values counts; i.e. values -> the # times they've been seen
  const auto m = ranges.size(); 
  while(!finished()){
    Rcpp::checkUserInterrupt();
    
    // Extract range with minimum first element
    auto rng_it = std::min_element(begin(ranges), end(ranges), [](auto& p1, auto &p2){
      if (p1.first == p1.second){ return false; }
      if (p2.first == p2.second){ return true; }
      return(*p1.first < *p2.first);
    });
    if (rng_it->first == rng_it->second){ return(false); }
    const size_t min_idx = std::distance(begin(ranges), rng_it); 
    auto rng = ranges[min_idx];
    
    // Partial sort up to the end
    std::nth_element(rng.first, safe_advance(rng.first, rng.second, k), rng.second);

    // Go through the current range 
    const auto e = safe_advance(rng.first, rng.second, k);
    for (auto b = rng.first; b != e; ++b){
      // auto value_it = std::lower_bound(begin(vc), end(vc), *b, [](auto& p, auto v){ return p.first < v; });
      auto value_it = vc.find(*b);
      if (value_it != end(vc) && ++(value_it->second) == m){
        mv.insert(value_it->first); // Increase number of things found in mutual
      } else {
        vc.emplace_hint(value_it, *b, 1);
      }
    }
    
    // If the common intersection grows to threshold, we're done
    if (mv.size() >= n){ return true;  }
    
    // Always replace beginning range with it's next beginning
    ranges[min_idx].first = e;
  }
  return(false);
}


// Tests a set of ranges to see if they all have at least n elements in their intersection.
template < typename Iter >
bool n_intersects(const vector< pair< Iter, Iter > >& ranges, const size_t n){
  using T = typename Iter::value_type; 
  
  // Do linear O(n) scan to determine if everything is sorted
  const bool is_sorted = std::all_of(begin(ranges), end(ranges), [](auto& rng){ return std::is_sorted(rng.first, rng.second); });
  
  // Collect (min,max) of each range
  auto minmaxes = vector< std::pair< T, T > >();
  minmaxes.reserve(ranges.size());
  std::transform(begin(ranges), end(ranges), std::back_inserter(minmaxes), [is_sorted](auto& rng){
    if (is_sorted){ return std::make_pair(*rng.first, *rng.second); }
    else {
      auto mm = std::minmax_element(rng.first, rng.second);
      return std::make_pair(*mm.first, *mm.second);
    }
  });
  
  // Do initial check to see if they are all disjoint
  if (intervals_disjoint(minmaxes)){
    return(false);
  }

  // Use the appropriate intersection check
  return is_sorted ? n_intersects_sorted(ranges, n) : n_intersects_unsorted(ranges, n);
}
      
      
      // // Add those less than the pivot to the id counts, checking for intersections
      // // The pivot element represents elements already checked.
      // // const size_t id_size = std::distance(begin(id_counts), end(id_counts));
      // for (auto it = b; it != b + inc; ++it){
      //   auto id_it = std::lower_bound(begin(id_counts), end(id_counts), *it, [](auto& p, const auto val){
      //     return(p.first < val);
      //   });
      //   if (id_it == end(id_counts)){
      //     std::cout << "emplacing " << *it << std::endl;
      //     id_counts.emplace_back(*it, 1); // element not found
      //   } else if (id_it->second+1 >= n_rngs){
      //     return true; // found an element in all the ranges!
      //   } else {
      //     std::cout << "found " << *it << " with count " << id_it->second << std::endl;
      //     ++id_it->second;
      //   }
      // }
    // }
    
    // const auto finished = [&offsets, &rng_sizes, n_rngs](){
    //   size_t n_finished = 0;
    //   for (size_t i = 0; i < n_rngs; ++i){
    //     n_finished += size_t(offsets[i] >= rng_sizes[i]);
    //   }
    //   std::cout << "n_finished: " << n_finished << std::endl;
    //   return n_finished >= 2;
    // };

    // Scan each range in sequence, adding ids to the multiset
  //   const T t_inf = std::numeric_limits< T >::infinity();
  //   auto id_counts = vector< pair< T, size_t > >();
  //   const auto first_less = [](auto& p1, auto& p2){ return p1.first < p2.first; };
  //   while(!finished()){
  //     Rcpp::checkUserInterrupt();
  // 
  //     std::cout << "getting min" << std::endl;
  //     size_t j = 0, i = 0;
  //     T min_elem = t_inf;
  //     for (auto& rng: ranges){
  //       T c_elem = offsets[j] >= rng_sizes[j] ? t_inf : *std::next(rng.first, offsets[j]);
  //       if (c_elem < min_elem){
  //         min_elem = c_elem;
  //         i = j;
  //       }
  //       ++j;
  //     }
  // 
  //     // Get subrange begin + end
  //     auto b = ranges[i].first + offsets[i];
  //     auto e = ranges[i].second;
  // 
  //     // move partially sorted elements around offsetted pivot
  //     auto inc = std::distance(b, e) < k ? std::distance(b, e) : k;
  //     std::nth_element(b, b + inc, e);
  // 
  //     // Add those less than the pivot to the id counts, checking for intersections
  //     // The pivot element represents elements already checked.
  //     // const size_t id_size = std::distance(begin(id_counts), end(id_counts));
  //     for (auto it = b; it != b + inc; ++it){
  //       auto id_it = std::lower_bound(begin(id_counts), end(id_counts), *it, [](auto& p, const auto val){
  //         return(p.first < val);
  //       });
  //       if (id_it == end(id_counts)){
  //         std::cout << "emplacing " << *it << std::endl;
  //         id_counts.emplace_back(*it, 1); // element not found
  //       } else if (id_it->second+1 >= n_rngs){
  //         return true; // found an element in all the ranges!
  //       } else {
  //         std::cout << "found " << *it << " with count " << id_it->second << std::endl;
  //         ++id_it->second;
  //       }
  //     }
  // 
  //     // Keep the vector sorted
  //     std::sort(begin(id_counts), end(id_counts), first_less);
  //     // std::sort(begin(id_counts) + id_size, end(id_counts), first_less);
  //     // std::inplace_merge(begin(id_counts), begin(id_counts) + id_size, end(id_counts));
  // 
  //     // otherwise increase the offset range
  //     offsets[i] = std::min(offsets[i]+k, rng_sizes[i]);
  //     Rcout << "offsets @ i" << offsets[i] << std::endl;
  //   }
  // }
  // return false;
// }

// [[Rcpp::export]]
bool nerve_cpp(SEXP st, vector< vector< int > > cover, const size_t n) {
  Rcpp::XPtr< SimplexTree > stp(st);
  using it_t = vector< int >::iterator; 
  
  auto ranges = vector< std::pair< it_t, it_t > >();
  std::transform(begin(cover), end(cover), std::back_inserter(ranges), [](auto& cs){
    return std::make_pair(begin(cs), end(cs));
  });
  
  bool is_connected = n_intersects(ranges, n);  
  return(is_connected);
}


/*** R
library(simplextree)
st <- simplex_tree()

## Test if there is at least a single element in common
nerve_cpp(st$as_XPtr(), list(c(1,2,3), c(3, 4, 5, 6)), 1)
# nerve_cpp(st$as_XPtr(), list(sample(seq(100)), sample(seq(150)), sample(seq(200))), 101)
# nerve_cpp(st$as_XPtr(), list(c(1,2,3), c(3, 4, 5, 6)), 1)
# nerve_cpp(st$as_XPtr(), list(c(2,3,1), c(5, 3, 4, 6), c(4, 6, 8, 3)), 1)
*/
