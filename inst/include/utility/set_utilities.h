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


// Moves through two ordered sets, returning a boolean indicating if they are disjoint
// returns true on first element found in both sets, otherwise iterates through both. 
template <typename SetA, typename SetB>
bool disjoint_sorted(const SetA &a, const SetB &b) {
  auto it_a = a.begin();
  auto it_b = b.begin();
  while (it_a != a.end() && it_b != b.end()) {
    switch (*it_a == *it_b ? 0 : *it_a < *it_b ? -1 : 1) {
    case 0:
      return false;
    case -1:
      it_a = std::lower_bound(++it_a, a.end(), *it_b);
      break;
    case 1:
      it_b = std::lower_bound(++it_b, b.end(), *it_a);
      break;
    }
  }
  return true;
}

// Given two random-access iterator ranges, (a1, a2), (b1, b2), return a boolean indicating 
// whether or not they have a non-empty intersection. Does not assumes either is sorted.
template <typename Iter>
bool nonempty(Iter a1, Iter a2, Iter b1, Iter b2){
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


template < typename Iter >
bool nfold_nonempty_sorted(vector< pair< Iter, Iter > > ranges){
  using T = typename Iter::value_type; 
  const size_t n_rngs = ranges.size(); 
  
  // Predicate to indicate if all ranges have been exhausted
  const auto finished = [&ranges](){
    bool all_ended = std::all_of(begin(ranges), end(ranges), [](auto& rng){
      return rng.first == rng.second;
    });
    return all_ended;
  };
  
  // Advance the iterator pointing to the smallest element 
  std::multiset< T > ids; 
  while(!finished()){
    auto smallest_elem = std::min_element(begin(ranges), end(ranges), [](auto& rng1, auto& rng2){
      if (rng1.first == rng1.second){ return false; }
      if (rng2.first == rng2.second){ return true; }
      return(*rng1.first < *rng2.first);
    });
    const size_t idx = std::distance(begin(ranges), smallest_elem);
    
    // Get the current valid id with the smallest value
    T id = *ranges[idx].first; 
    ids.insert(id);
    
    // If the id just inserted has appeared at least 'n_rngs' times, then we've found a nonempty intersection
    if (ids.count(id) >= n_rngs){
      return(true);
    }
    
    // otherwise increment the smallest iterator 
    std::advance(ranges[idx].first, 1);
  }
  return(false);
}

// Higher order generalization of nonempty. Checks to see if any number of ranges all have a nonempty intersection.
// The ranges are assumed not to be sorted apriori.
template < typename Iter >
bool n_nonempty_unsorted(vector< pair< Iter, Iter > > ranges, const size_t k = 32){
  using T = typename Iter::value_type; 
  const size_t n_rngs = ranges.size();  
  vector< pair< Iter, Iter > > current_range;
  
  // todo: do minmax on each range, stuff into pairs w/ indices, sort, then check adjacency 
  // to see if there are any disjoint. they are in non-overlapping ranges
  // Also: remove ranges that *are* adjacent and not-enclosed by other ranges
  
  // Do minmax on each range 
  vector< std::pair< T, T > > minmaxes_idx;
  T i = 0; 
  for (auto& rng: ranges){
    auto mm = std::minmax_element(rng.first, rng.second);
    minmaxes_idx.push_back( std::make_pair(i, *mm.first) );
    minmaxes_idx.push_back( std::make_pair(i, *mm.second) );
    ++i;
  }
  
  // Sort by the values of extrema  
  std::stable_sort(begin(minmaxes_idx), end(minmaxes_idx), [](auto& p1, auto& p2){
    return p1.second < p2.second; 
  });

  // Search for adjacent indices
  auto adj = std::adjacent_find(begin(minmaxes_idx), end(minmaxes_idx), [](auto& p1, auto& p2){
    return p1.first == p2.first;
  });
  
  // If any adjacent elements were found, then there's at least one interval disjoint from the rest
  if (adj != end(minmaxes_idx)){
    bool is_disjoint = std::is_sorted(begin(minmaxes_idx), end(minmaxes_idx), [](auto& p1, auto& p2){
      return p1.first == p2.first ? p1.second <= p2.second : p1.second < p2.second;
    });
    if (is_disjoint){ return(false); }
  } 
  
  // Otherwise, sort by size, then fold a set_intersection 

  // Get range sizes 
  vector< size_t > rng_sizes(n_rngs);
  std::transform(begin(ranges), end(ranges), begin(rng_sizes), [](auto& rng){
    return std::distance(rng.first, rng.second);
  });
  
  // Create a set of initial offsets giving beginnings
  vector< size_t > offsets(n_rngs, 0);
  
  // Use lambda to determine when finished
  const auto finished = [&offsets, &rng_sizes, n_rngs](){
    size_t n_finished = 0; 
    for (size_t i = 0; i < n_rngs; ++i){
      n_finished += size_t(offsets[i] >= rng_sizes[i]);
    }
    std::cout << "n_finished: " << n_finished << std::endl;
    return n_finished >= 2; 
  };
  
  // Scan each range in sequence, adding ids to the multiset
  const T t_inf = std::numeric_limits< T >::infinity();
  auto id_counts = vector< pair< T, size_t > >(); 
  const auto first_less = [](auto& p1, auto& p2){ return p1.first < p2.first; };
  while(!finished()){
    
    std::cout << "getting min" << std::endl;
    size_t j = 0, i = 0; 
    T min_elem = t_inf;
    for (auto& rng: ranges){
      T c_elem = offsets[j] >= rng_sizes[j] ? t_inf : *std::next(rng.first, offsets[j]);
      if (c_elem < min_elem){
        min_elem = c_elem;
        i = j;
      }
      ++j;
    }
    
    // Get subrange begin + end 
    auto b = ranges[i].first + offsets[i];
    auto e = ranges[i].second;

    // move partially sorted elements around offsetted pivot 
    auto inc = std::distance(b, e) < k ? std::distance(b, e) : k;
    std::nth_element(b, b + inc, e);
    
    // Add those less than the pivot to the id counts, checking for intersections
    // The pivot element represents elements already checked.  
    // const size_t id_size = std::distance(begin(id_counts), end(id_counts));
    for (auto it = b; it != b + inc; ++it){
      auto id_it = std::lower_bound(begin(id_counts), end(id_counts), *it, [](auto& p, const auto val){
        return(p.first < val);
      });
      if (id_it == end(id_counts)){
        std::cout << "emplacing " << *it << std::endl;
        id_counts.emplace_back(*it, 1); // element not found 
      } else if (id_it->second+1 >= n_rngs){
        return true; // found an element in all the ranges!
      } else {
        std::cout << "found " << *it << " with count " << id_it->second << std::endl;
        ++id_it->second;
      }
    }
    
    // Keep the vector sorted
    std::sort(begin(id_counts), end(id_counts), first_less);
    // std::sort(begin(id_counts) + id_size, end(id_counts), first_less);
    // std::inplace_merge(begin(id_counts), begin(id_counts) + id_size, end(id_counts));

    // otherwise increase the offset range 
    offsets[i] = std::min(offsets[i]+k, rng_sizes[i]);
  }
  return false; 
}
    
  //   for (size_t i = 0; i < n_rngs; ++i){
  //     auto b = ranges[i].first;
  //     auto e = ranges[i].second; 
  //     
  //     // partial sort around a pivot 
  //     std::nth_element(b, b + std::min(offsets[i], rng_sizes[i]), e);
  //     
  //     // Add those less than the pivot to the multiset, checking for intersections
  //     for (auto it = b; it != b + offsets[i]; ++it){
  //       const size_t m = ids.count(*it);
  //       if ((m+1) >= n_rngs){
  //         return true; // found an element in all the ranges!
  //       }
  //     }
  //     // otherwise increase the offset range 
  //     offsets[i] += k;
  //   }
  // }
// }

// Higher order generalization of nonempty. Checks to see if any number of ranges all have a nonempty intersection.
template <typename Iter>
bool nfold_nonempty(vector< std::pair< Iter, Iter > > ranges){
  using it_cat = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<it_cat, std::random_access_iterator_tag>::value, "Iterator type must be random-access.");
  using T = typename std::iterator_traits<Iter>::value_type;

  // Either empty == they do not have an intersection
  const size_t n_rng = ranges.size();
  vector< size_t > rng_sizes = vector< size_t >(n_rng);
  std::transform(begin(ranges), end(ranges), begin(rng_sizes), [](const std::pair<Iter, Iter> rng){
    return std::distance(rng.first, rng.second);
  });
  bool any_empty = std::any_of(begin(rng_sizes), end(rng_sizes), [](const size_t sz){ return sz == 0; });
  if (any_empty){ return(false); };

  // Use multiset to track number of ids
  std::multiset<T> ids; 
  const auto insert_rng = [&ids](Iter a, Iter b){ 
    std::for_each(a, b, [&ids](T elem){ ids.insert(elem); }); 
  };
  for (auto rng: ranges){ insert_rng(rng.first, rng.second); };
  
  // If any ids appeared k times, there's a nonempty intersection between them all 
  bool nonempty_int = std::any_of(ids.begin(), ids.end(), [&ids, &n_rng](const T id){
    return (ids.count(id) == n_rng);
  });
  return nonempty_int;
}

template <typename Iter>
auto intersection(Iter a1, Iter a2, Iter b1, Iter b2) 
  -> vector< typename std::iterator_traits<Iter>::value_type > {
  using T = typename std::iterator_traits<Iter>::value_type;
  std::unordered_set< T > c;
  for (; a1 != a2; ++a1) {
    if (std::find(b1, b2, *a1) != b2) { 
      c.insert(*a1); 
    }
  }
  vector< T > res(begin(c), end(c));
  return(res);
}

template <typename Iter>
vector< int > nfold_intersection(vector< std::pair< Iter, Iter > > ranges){
  using it_cat = typename std::iterator_traits<Iter>::iterator_category;
  static_assert(std::is_same<it_cat, std::random_access_iterator_tag>::value, "Iterator type must be random-access.");
  using T = typename std::iterator_traits<Iter>::value_type;
  vector< T > res;
  const size_t n_pairs = ranges.size()-1; 
  for (size_t i = 0; i < n_pairs; ++i){
    const std::pair< Iter, Iter > rng1 = ranges[i];
    const std::pair< Iter, Iter > rng2 = ranges[i+1];
    vector< T > c_int = intersection(rng1.first, rng1.second, rng2.first, rng2.second);
    res.insert(res.end(), begin(c_int), end(c_int));
    if (res.size() == 0){
      return(res);
    }
  }
  return(res);
}
