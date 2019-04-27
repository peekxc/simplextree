#ifndef UTILITIES_H_
#define UTILITIES_H_

// Aliases
using std::vector;
using std::size_t;
using std::begin; 
using std::end;
template <typename T> using s_ptr = std::shared_ptr<T>; // Shared pointer
template <typename T> using u_ptr = std::unique_ptr<T>; // Unique pointer

// To use alloca portably
#include <cstdlib> // alloca
#ifdef __GNUC__
/* Includes GCC, clang and Intel compilers */
# undef alloca
# define alloca(x) __builtin_alloca((x))
#elif defined(__sun) || defined(_AIX)
/* this is necessary (and sufficient) for Solaris 10 and AIX 6: */
# include <alloca.h>
#endif

// Identity functor. Invoke using identity() for functional algorithms requiring unary functions
struct identity {
  template<typename U>
  constexpr auto operator()(U&& v) const noexcept
    -> decltype(std::forward<U>(v)) {
      return std::forward<U>(v);
    }
};

// Cantor pairing function. Takes as input two unsigned integral types (a, b), and uniquely 
// maps (a, b) to a number c, where c is possibly a different integral type 
template <typename T1, typename T2> 
inline T2 cantor_pair(T1 a, T1 b){
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  return static_cast<T2>(a (a + b) * (a + b + 1) / 2 + a);
}

// Szudziks pairing function. Takes as input two unsigned integral types (a, b), and uniquely 
// maps (a, b) to a number c, where c is possibly a different integral type 
template <typename T1, typename T2> 
inline T2 szudzik_pair(T1 a, T1 b){
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  return static_cast<T2>(a >= b ? a * a + a + b : a + b * b);
}


// Given two indices 'from' and 'to', each of which are between [0, 1, ..., N-1], this 
// function returns the corresponding 'flat' 0-based index the entry would appear in the  
// lower triangular portion of an (N x N) column-major matrix.
template <typename T> 
inline T index_lt(T from, T to, const T N){
  if (from < to){ std::swap(from, to); }
  return((N)*(to)-(to)*(to+1)/2+from-to-1);
}

// Given a flat 0-based index 'k' and the size 'N', returns the 0-based column index 'k' would represent 
// if 'k' indexed into the lower triangular portion of an (N x N) column-major matrix.
template <typename T> 
inline T index_to(T k, T n){
  return n - 2 - std::floor(std::sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5); // expects 0-based, returns 0-based 
}

// Given a flat 0-based index 'k', the size 'N', and the column index 'i', returns the 0-based row index 'k' 
// would represent if 'k' indexed into the lower triangular portion of an (N x N) column-major matrix.
template <typename T> 
inline T index_from(T k, T n, T i){
  return k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2; // expects 0-based, returns 0-based 
}

// Creates a vector with the range [i, j]
template <typename T>
vector<T> seq_ij(const T i, const T j){
  static_assert(std::is_integral<T>::value, "Integral-type required as a range storage type.");
  size_t sz = i < j ? (j - i)+1 : (i - j)+1;
  vector<T> rng = vector<T>(sz);
  std::iota(rng.begin(), rng.end(), static_cast<T>(i));
  return(rng);
}

// Compile-time enumeration that allows generic weak order
enum Ordering { EQ = 0, LT = -1, GT = 1 };

// Comparison function 
template <typename A, typename B>
inline Ordering compare(const A &a, const B &b) {
  return a == b ? EQ : a < b ? LT : GT;
}

// Given two sorted or partitioned containers (a, b), traverses each container, looking to see 
// if they are disjoint. If one equivalent element is found that lies in both containers, the 
// function ends returning false. This should be faster than set_difference, set_intersection, etc.
template <typename SetA, typename SetB>
bool is_disjoint(const SetA &a, const SetB &b) {
  auto it_a = a.begin();
  auto it_b = b.begin();
  while (it_a != a.end() && it_b != b.end()) {
    switch (compare(*it_a, *it_b)) {
    case EQ:
      return false;
    case LT:
      it_a = std::lower_bound(++it_a, a.end(), *it_b);
      break;
    case GT:
      it_b = std::lower_bound(++it_b, b.end(), *it_a);
      break;
    }
  }
  return true;
}


// Applies the lambda function f to all pairwise combinations in the range [first, last)
template<typename Iter, typename Lambda>
inline void combine_pairwise(Iter first, Iter last, Lambda f) {
  for(; first != last; ++first){
    for(Iter next = std::next(first); next != last; ++next){
      f(*first, *next);
    }
  }
}

// Implements a generic n-vector cartesian product for a vector of vectors using an iterative design pattern. 
// Individual items are put into a fixed-sized vector and given to a passed in Lambda function. Limited to products 
// of vectors having the same type. Original design based on: 
// https://stackoverflow.com/questions/18732974/c-dynamic-number-of-nested-for-loops-without-recursion/30808351
template <typename T, typename Func> 
inline void CartesianProduct(const vector< vector<T> >& elems, Func&& f) {
  
  // Initialize the slots to hold the current iteration value for each depth
  const size_t depth = elems.size();
  size_t* slots = (size_t*) alloca(sizeof(size_t) * depth);
  for (size_t i = 0; i < depth; i++) { slots[i] = 0; }
  
  // Extract the sizes of each vector in the product
  vector<size_t> max = vector<size_t>(depth);
  std::transform(elems.begin(), elems.end(), max.begin(), [](const vector<T>& lst){ return(lst.size()); });
  vector<T> current_element(depth);
  
  size_t index = 0, i = 0;
  while (true) {
    
    // Fill the element and apply the lambda 
    i = 0; 
    std::transform(elems.begin(), elems.end(), current_element.begin(), [&i, &slots](const vector<T>& v){
      return(v.at(slots[i++]));
    });
    f(current_element);
    
    // Increment
    slots[0]++;
    
    // Carry
    while (slots[index] == max.at(index)) {
      if (index == depth - 1) { return; } // Overflow, we're done
      slots[index++] = 0;
      slots[index]++;
    }
    index = 0;
  }
}

// Alternative ordering of the cartesian product 
template <typename T, typename Lambda> 
inline void CartesianProductOrdered(const vector< vector<T> >& elems, Lambda&& f, const vector< size_t > order) {
  // Initialize the slots to hold the current iteration value for each depth
  const size_t depth = elems.size();
  size_t* slots = (size_t*) alloca(sizeof(size_t) * depth);
  for (size_t i = 0; i < depth; i++) { slots[i] = 0; }
  
  // Extract the sizes of each vector in the product
  vector<size_t> max = vector<size_t>(depth);
  std::transform(elems.begin(), elems.end(), max.begin(), [](const vector<T>& lst){ return(lst.size()); });
  vector<T> current_element(depth);
  
  size_t index = order.at(0), i = 0, o_idx = 0;
  while (true) {
    
    // Fill the element and apply the lambda 
    i = 0; 
    std::transform(elems.begin(), elems.end(), current_element.begin(), [&i, &slots](const vector<T>& v){
      return(v.at(slots[i++]));
    });
    f(current_element);
    
    // Increment + Carry
    slots[order.at(0)]++;
    while (slots[index] == max.at(index)) {
      if (index == order.back()) { return; } // Overflow, we're done
      slots[index] = 0;
      o_idx++;
      index = order.at(o_idx);
      slots[index]++;
    }
    o_idx = 0; 
    index = order.at(o_idx);
  }
}

// Generates the set of all (n choose k) combinations, passing each combination 
// to the lambda function f
template <typename Lambda>
bool apply_combinations(const size_t n, const size_t k, Lambda f){
  vector< bool > v(n);
  vector< size_t > idx(k);
  std::fill(v.end() - k, v.end(), true);
  bool all_eval = true; 
  do {
    for (size_t i = 0, cc = 0; i < n; ++i) { 
      if (v[i]) { idx[cc++] = i; }
    }
    all_eval = bool(f(idx));
  } while (all_eval && std::next_permutation(v.begin(), v.end()));
  return(all_eval);
}

// Given two random-access iterator ranges, (a1, a2), (b1, b2), return a boolean indicating 
// whether or not they have a non-empty intersection. Does not assumes either is sorted.
template <typename Iter>
bool nonempty_intersection(Iter a1, Iter a2, Iter b1, Iter b2){
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
    }
    return(false);
  } else if (b_sz * 100 < a_sz) { // Opposite case
    vector<T> a_sort(a_sz);
    std::partial_sort_copy(a1, a2, begin(a_sort), end(a_sort)); // partial-sorted elements of y copied to y_sort
    while (a1 != a2){
      if (std::binary_search(begin(a_sort), end(a_sort), T(*b1))) { return(true); }
    }
    return(false);
  }
  
  // Otherwise, sort both, then use lower_bound type approach to potentially skip massive sections
  vector<T> a_sort(a_sz), b_sort(b_sz);
  std::partial_sort_copy(a1, a2, begin(a_sort), end(a_sort)); // partial-sorted elements of y copied to y_sort
  std::partial_sort_copy(b1, b2, begin(b_sort), end(b_sort)); // partial-sorted elements of y copied to y_sort
  return !is_disjoint(a_sort, b_sort);
}


template<typename ForwardIterator>
inline auto unique_indices(ForwardIterator first, ForwardIterator last) 
  -> decltype(std::map< typename std::iterator_traits<ForwardIterator>::value_type, size_t>()){
  using value_type = typename std::iterator_traits<ForwardIterator>::value_type;
  std::map< value_type, size_t > val_idx_map; // map between values (keys) and indices (size_t)
  for(size_t i = 0; first != last; ++i, ++first){
    auto it = val_idx_map.find(*first);
    if (it == val_idx_map.end()) { // value doesn't exist
      val_idx_map.emplace(*first, i);
    }
  }
  return val_idx_map;
}

template <typename T> 
vector<T> merge_vectors(const vector< vector<T>* >& vec){
  size_t total_vec_size = 0;
  for (auto& v: vec){ total_vec_size += (*v).size(); }
  
  // Concatenate all the vectors into a new vector
  std::vector< T > final_res = std::vector< T >();
  final_res.reserve(total_vec_size);
  for (auto& v: vec){ final_res.insert(final_res.end(), v->begin(), v->end()); }
  return(final_res);
}

// Erase-remove idiom
template <typename T, typename Iter> 
void erase_remove(vector<T>& v, Iter b, Iter e, T val){
  v.erase(std::remove(b, e, val), v.end());
}

// Erase-partition idiom
template <typename T, typename Iter, typename Lambda> 
void erase_partition(vector<T>& v, Iter b, Iter e, Lambda f){
  v.erase(std::partition(b, e, f), v.end());
}


// Combinations chooser
template <typename T>
size_t choose(T n, T k) {
  static_assert(std::is_integral<T>::value, "Integral-type required.");
  if (k > n) { return 0; }
  if (k * 2 > n) { k = n-k; }
  if (k == 0) { return 1; }
  
  size_t result = n;
  for(T i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}


#endif 