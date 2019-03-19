// utility_rcpp.h
// Utility functions using Rcpp. 
#ifndef UTIL_RCPP_H
#define UTIL_RCPP_H

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]   

#include <memory> //shared_ptr

// Choice of base types to use for indexing/storage
using uint = unsigned int;                      // unsigned int
using uint8_t = uint_fast8_t;                   // 8+ bit signed integer type
using sint8 = int_fast8_t;                      // 8+ bit signed integer type
using uidx_t = std::size_t;                     // signed index type 
using sidx_t = std::ptrdiff_t;                  // unsigned index type 
template <typename T> 
using s_ptr = std::shared_ptr<T>;               // Shared pointer
template <typename T> 
using u_ptr = std::unique_ptr<T>;               // Unique pointer

template <typename T>
using enable_int = typename std::enable_if<std::is_integral<T>::value>;

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

// Nice simple interface from: https://stackoverflow.com/questions/17973442/index-element-from-list-in-rcpp
template <typename WHAT>
class ListOf : public List {
public:
  template <typename T>
  ListOf( const T& x) : List(x){}
  WHAT operator[](int i){ return as<WHAT>( ( (List*)this)->operator[]( i) ) ; }
  WHAT at(int i){ return as<WHAT>( ( (List*)this)->at(i) ) ; }
};

// Given a flat 0-based index 'k' and the size 'N', returns the 0-based column index 'k' would represent 
// if 'k' indexed into the lower triangular portion of an (N x N) column-major matrix.
#ifndef INDEX_TO
#define INDEX_TO(k, n) n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5) // expects 0-based, returns 0-based
#endif

// Given a flat 0-based index 'k', the size 'N', and the column index 'i', returns the 0-based row index 'k' 
// would represent if 'k' indexed into the lower triangular portion of an (N x N) column-major matrix.
#ifndef INDEX_FROM
#define INDEX_FROM(k, n, i) k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 // expects 0-based, returns 0-based
#endif

// namespace util {

  sidx_t index_lower_triangular(sidx_t from, sidx_t to, const sidx_t N);

  template <typename T>
  std::vector<T> seq_ij(const T i, const T j);
  
  template<typename ForwardIterator>
  std::map<int, int> get_unique_indices(ForwardIterator first, ForwardIterator last); 
  
  template <typename T> 
  std::vector<T> merge_vectors(const std::vector< std::vector<T>* >& vec);
    
  // Rcpp-specific utilities 
  IntegerMatrix make_cartesian_product(const List& vecs);
  bool any_is_in(const IntegerVector& x, const IntegerVector& y);
  IntegerMatrix rbindlist_int(std::list<IntegerVector>& lst);
  List resize_list(const List& x, int n);
  template <typename T> IntegerVector to_ivec(std::vector<T> vec);
// }

// utility_rcpp.hpp
// Contains pure-template implementations suitable for being inlined.  

// Applies the function Func to all pairwise combinations in the range [first, last)
template<typename Iter, typename Func>
inline void combine_pairwise(Iter first, Iter last, Func func)
{
  for(; first != last; ++first){
    for(Iter next = std::next(first); next != last; ++next){
      func(*first, *next);
    }
  }
}
// Implements a generic n-vector cartesian product for a vector of vectors using an iterative design pattern. 
// Individual items are put into a fixed-sized vector and given to a passed in Lambda function. Limited to products 
// of vectors having the same type. Original design based on: 
// https://stackoverflow.com/questions/18732974/c-dynamic-number-of-nested-for-loops-without-recursion/30808351
template <typename T, typename Func> 
inline void CartesianProduct(const std::vector< std::vector<T> >& elems, Func&& f) {
  
  // Initialize the slots to hold the current iteration value for each depth
  const std::size_t depth = elems.size();
  std::size_t* slots = (std::size_t*) alloca(sizeof(std::size_t) * depth);
  for (std::size_t i = 0; i < depth; i++) { slots[i] = 0; }
  
  // Extract the sizes of each vector in the product
  std::vector<std::size_t> max = std::vector<std::size_t>(depth);
  std::transform(elems.begin(), elems.end(), max.begin(), [](const std::vector<T>& lst){ return(lst.size()); });
  std::vector<T> current_element(depth);
  
  std::size_t index = 0, i = 0;
  while (true) {
    
    // Fill the element and apply the lambda 
    i = 0; 
    std::transform(elems.begin(), elems.end(), current_element.begin(), [&i, &slots](const std::vector<T>& v){
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
template <typename T, typename Func> 
inline void CartesianProductOrdered(const std::vector< std::vector<T> >& elems, Func&& f, const std::vector< std::size_t > order) {
  // Initialize the slots to hold the current iteration value for each depth
  const std::size_t depth = elems.size();
  std::size_t* slots = (std::size_t*) alloca(sizeof(std::size_t) * depth);
  for (std::size_t i = 0; i < depth; i++) { slots[i] = 0; }
  
  // Extract the sizes of each vector in the product
  std::vector<std::size_t> max = std::vector<std::size_t>(depth);
  std::transform(elems.begin(), elems.end(), max.begin(), [](const std::vector<T>& lst){ return(lst.size()); });
  std::vector<T> current_element(depth);
  
  std::size_t index = order.at(0), i = 0, o_idx = 0;
  while (true) {
    
    // Fill the element and apply the lambda 
    i = 0; 
    std::transform(elems.begin(), elems.end(), current_element.begin(), [&i, &slots](const std::vector<T>& v){
      return(v.at(slots[i++]));
    });
    f(current_element);
    
    // Increment
    slots[order.at(0)]++;
    
    // Carry
    while (slots[index] == max.at(index)) {
      if (index == order.back()) { return; } // Overflow, we're done
      slots[ index ] = 0;
      o_idx++;
      index = order.at(o_idx);
      slots[index]++;
    }
    o_idx = 0; 
    index = order.at(o_idx);
  }
}

#endif /* UTIL_RCPP_H */

