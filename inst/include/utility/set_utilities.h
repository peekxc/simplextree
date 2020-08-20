// nerve_utility.cpp
// Contains utility functions related to the nerve construction

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <iostream>

using std::vector; 
using std::begin; 
using std::end;
using std::pair;
using std::unordered_map;
using std::unordered_set;


// Safely advance iterator to prevent passing the end
template <class Iter, class Incr>
Iter safe_advance(const Iter& curr, const Iter& end, Incr n) {
  size_t remaining(std::distance(curr, end));
  if (remaining < size_t(n)) { n = remaining; }
  return std::next(curr, n);
}

// Moves through two ordered sets, returning a boolean indicating if they are disjoint
// returns true on first element found in both sets, otherwise iterates through both. 
template <typename Iter>
bool disjoint_sorted(Iter it_a, const Iter a_end, Iter it_b, const Iter b_end) {
  while (it_a != a_end && it_b != b_end) {
    switch (*it_a == *it_b ? 0 : *it_a < *it_b ? -1 : 1) {
    case 0:
      return false;
    case -1:
      it_a = std::lower_bound(++it_a, a_end, *it_b);
      break;
    case 1:
      it_b = std::lower_bound(++it_b, b_end, *it_a);
      break;
    }
  }
  return true;
}

// Given a set of intervals as pairs, checks if any of them are disjoint from the others
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
  using rng_t = std::pair< T, T >; 
  std::stable_sort(begin(interval_ids), end(interval_ids), [](const rng_t& p1, const rng_t& p2){
    return p1.second < p2.second; 
  });
  
  // Check if any adjacent values are equal
  auto adj_it = std::adjacent_find(begin(interval_ids), end(interval_ids), [](const rng_t& p1, const rng_t& p2){
    return(p1.second == p2.second);
  });
  if (adj_it != end(interval_ids)){ return(false); }
  
  // Check if sequence is strictly increasing
  auto sinc_it = std::adjacent_find(begin(interval_ids), end(interval_ids), [](const rng_t& p1, const rng_t& p2){
    return(p1.first > p2.first);
  });
  if (sinc_it != end(interval_ids)){ return(false); }
  
  // Otherwise they are disjoint
  return(true);
}

// Given two random-access iterator ranges, (a1, a2), (b1, b2), return a boolean indicating 
// whether or not they have a non-empty intersection. Does not assumes either is sorted.
template <typename Iter>
bool intersects_nonempty(Iter a1, Iter a2, Iter b1, Iter b2){
  using it_cat = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<it_cat, std::random_access_iterator_tag>::value, "Iterator type must be random-access."); 
  using T = typename std::iterator_traits<Iter>::value_type;
    
  // Either empty == they do not have an intersection
  const size_t a_sz = std::distance(a1, a2); 
  const size_t b_sz = std::distance(b1, b2); 
  if (a_sz == 0 || b_sz == 0) { return false; }
  
  // a is much smaller than b => partial sort b, then do binary search on b for each element of a
  if (a_sz * 100 < b_sz) {
    vector<T> b_sort(b_sz);
    std::partial_sort_copy(b1, b2, begin(b_sort), end(b_sort)); // partial-sorted elements of y copied to y_sort
    while (a1 != a2){
      if (std::binary_search(begin(b_sort), end(b_sort), T(*a1))) { return(true); }
      ++a1;
    }
    return(false);
  } else if (b_sz * 100 < a_sz) { // Opposite case
    vector<T> a_sort(a_sz);
    std::partial_sort_copy(a1, a2, begin(a_sort), end(a_sort)); // partial-sorted elements of y copied to y_sort
    while (b1 != b2){
      if (std::binary_search(begin(a_sort), end(a_sort), T(*b1))) { return(true); }
      ++b1;
    }
    return(false);
  }
  
  // Otherwise, sort both, then use lower_bound type approach to potentially skip massive sections.
  vector<T> a_sort(a_sz), b_sort(b_sz);
  std::partial_sort_copy(a1, a2, begin(a_sort), end(a_sort)); // partial-sorted elements of y copied to y_sort
  std::partial_sort_copy(b1, b2, begin(b_sort), end(b_sort)); // partial-sorted elements of y copied to y_sort
  return !disjoint_sorted(a_sort, b_sort);
}

// Checks if there are at least n elements common to all given sorted ranges
template < typename Iter >
bool n_intersects_sorted(vector< pair< Iter, Iter > > ranges, const size_t n){
  using T = typename Iter::value_type; 
  if (n == 0){ return(true); }
  if (ranges.size() <= 1){ return(false); }
  
  // Sort by size, then fold a set_intersection
  using rng_t = pair< Iter, Iter >;
  std::sort(begin(ranges), end(ranges), [](rng_t& p1, rng_t& p2){
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
  using rng_t = pair< Iter, Iter >;
  const auto finished = [&ranges](){
    size_t n_finished = std::accumulate(begin(ranges), end(ranges), (size_t) 0, [](size_t val, rng_t& p){ return p.first == p.second; });
    return(n_finished >= 2);
  };
  
  auto mv = std::set< T >();          // the values mutally common to all the ranges
  auto vc = std::map< T, size_t >();  // values counts; i.e. values -> the # times they've been seen
  const auto m = ranges.size(); 
  while(!finished()){

    // Partial sort to put minimum element at the beginning of each range
    for (auto& rng: ranges){
      std::nth_element(rng.first, safe_advance(rng.first, rng.second, 1), rng.second);
    }
    
    // Extract range with minimum first element
    auto rng_it = std::min_element(begin(ranges), end(ranges), [](rng_t& p1, rng_t &p2){
      if (p1.first == p1.second){ return false; }
      if (p2.first == p2.second){ return true; }
      return(*p1.first < *p2.first);
    });
    
    // If the range with the minimum element is an empty range, no candidate ranges available
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
  using rng_t = pair< Iter, Iter >;
  
  // Check if any of the ranges don't even have n elements
  const bool too_small = std::any_of(begin(ranges), end(ranges), [n](const rng_t& rng){ return size_t(std::distance(rng.first, rng.second)) < size_t(n); });
  if (too_small){ return(false); }
  
  // Do linear O(n) scan to determine if everything is sorted
  const bool is_sorted = std::all_of(begin(ranges), end(ranges), [](const rng_t& rng){ return std::is_sorted(rng.first, rng.second); });
  
  // Collect (min,max) of each range
  auto minmaxes = vector< std::pair< T, T > >();
  minmaxes.reserve(ranges.size());
  std::transform(begin(ranges), end(ranges), std::back_inserter(minmaxes), [is_sorted](const rng_t& rng){
    if (is_sorted){ 
      auto min = *rng.first;
      auto max = std::distance(rng.first, rng.second) == 1 ? *rng.first : *std::prev(rng.second);
      return std::make_pair(min,max); 
    } else {
      auto mm = std::minmax_element(rng.first, rng.second);
      return std::make_pair(*mm.first, *mm.second);
    }
  });
  
  // Do initial check to see if they are all disjoint
  if (intervals_disjoint(minmaxes)){ 
    return(false); 
  }

  // Use the appropriate intersection check
  if (is_sorted && n == 1 && ranges.size() == 2){
    bool disjoint = disjoint_sorted(ranges[0].first, ranges[0].second, ranges[1].first, ranges[1].second);
    return(!disjoint);
  } else {
    return is_sorted ? n_intersects_sorted(ranges, n) : n_intersects_unsorted(ranges, n);
  }
}
