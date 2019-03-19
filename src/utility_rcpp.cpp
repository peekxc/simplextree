// utility_rcpp.cpp
// Utility functions using Rcpp. 

#include "utility_rcpp.h"

// namespace util {

  // Given two indices 'from' and 'to', each of which are between [0, 1, ..., N-1], this 
  // function returns the corresponding 'flat' 0-based index the entry would appear in the  
  // lower triangular portion of an (N x N) column-major matrix.
  sidx_t index_lower_triangular(sidx_t from, sidx_t to, const sidx_t N){
    if (from < to){ std::swap(from, to); }
    return((N)*(to) - (to)*(to+1)/2 + (from) - (to) - (1));
  }
  
  // Creates a vector with the range [i, j]
  template <typename T>
  std::vector<T> seq_ij(const T i, const T j){
    static_assert(std::is_integral<T>::value, "Integral-type required as a range storage type.");
    // static_assert(std::is_integral<C>::value, "Integral-type required as a range type.");
    std::size_t sz = std::abs(j - i)+1;
    std::vector<T> rng = std::vector<T>(sz);
    std::iota(rng.begin(), rng.end(), static_cast<T>(i));
    return(rng);
  }
  
  template <typename T> 
  std::vector<T> merge_vectors(const std::vector< std::vector<T>* >& vec) {
    std::size_t total_vec_size = 0;
    std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>* v){ total_vec_size += v->size(); });
    std::vector< T > final_res = std::vector< T >();
    final_res.reserve(total_vec_size);
    // std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>* v){ std::copy(v->begin(), v->end(), std::back_inserter(final_res)); });
    std::for_each(vec.begin(), vec.end(), [&](const std::vector<T>* v){
      final_res.insert(final_res.end(), v->begin(), v->end());
    });
    return(final_res);
  }

  
  template<typename ForwardIterator>
  inline std::map<int, int> get_unique_indices(ForwardIterator first, ForwardIterator last){
    std::map<int, int> pt_to_unique_idx;
    for(std::size_t i = 0; first != last; ++i, ++first){
      auto it = pt_to_unique_idx.find(*first);
      if (it == pt_to_unique_idx.end()) { // value doesn't exist
        pt_to_unique_idx.emplace(*first, i);
      }
    }
    return pt_to_unique_idx;
  }
  
  template <typename T> 
  IntegerVector to_ivec(std::vector<T> vec){
    static_assert(std::is_integral<T>::value, "T must be integral type");
    IntegerVector res = IntegerVector(vec.size());
    std::size_t i = 0; 
    std::for_each(vec.begin(), vec.end(), [&i, &res](const T val){
      res.at(i++) = static_cast<int>(val);
    });
    return(res);
  }
  
  // Fast partial-sort/binary-search check to see if the intersection between two given vectors has non-zero length
  // Loosely based off of bugged version found here: https://stackoverflow.com/questions/21359432/a-c-version-of-the-in-operator-in-r
  bool any_is_in(const IntegerVector& x, const IntegerVector& y){
    std::vector<int> y_sort(y.size());
    std::partial_sort_copy (y.begin(), y.end(), y_sort.begin(), y_sort.end()); // partial-sorted elements of y copied to y_sort
    const std::size_t nx = x.size();
    for (std::size_t i = 0; i < nx; ++i) {
      if (std::binary_search(y_sort.begin(), y_sort.end(), x[i])) {
        return(true); // end the search
      }
    }
    return(false);
  }
  
  // rbindlist_int: Takes a list of integer vectors and rbind's them together.
  IntegerMatrix rbindlist_int(std::list<IntegerVector>& lst){
    std::size_t n = lst.size();
    if(n == 0) { Rcpp::stop("Invalid sized list."); }
    std::size_t d = lst.front().size();
    Rcpp::IntegerMatrix res = Rcpp::no_init(n, d);
    std::size_t i = 0;
    for (std::list<IntegerVector>::iterator it = lst.begin(); it != lst.end(); ++it, ++i) {
      if (static_cast<uidx_t>((*it).size()) != d) { Rcpp::stop("Invalid integer vector size."); }
      res(i, _) = *it;
    }
    return res;
  }
  
  // Resizes a list 
  List resize_list(const List& x, int n){
    const std::size_t lst_sz = x.size();
    List y(n);
    for(std::size_t i = 0; i < lst_sz; i++) { y[i] = x[i]; }
    return(y);
  }
  
  
  IntegerMatrix make_cartesian_product(const List& vecs){
    std::vector< std::vector<int> > vv(vecs.size());
    std::size_t n_rows = 1;
    uidx_t vsize = vecs.size();
    for (std::size_t i = 0; i < vsize; ++i){
      IntegerVector v = vecs.at(i);
      vv.at(i) = as< std::vector<int> >(v);
      n_rows = n_rows * v.size();
    }
    
    // Copy to result
    std::size_t i = 0;
    IntegerMatrix res = IntegerMatrix(n_rows, vecs.size());
    CartesianProduct(vv, [&i, &res](std::vector<int> idx){
      res(i++, _) = as<IntegerVector>(wrap(idx));
    });
    return(res);
  }

  // template <typename T> to_ivec< enable_int::type >();
  
// } // util namespace 


template IntegerVector to_ivec<uint8_t>(std::vector<uint8_t> vec);
template std::vector<uint8_t> seq_ij<uint8_t>(uint8_t i, uint8_t j);
template std::vector<sidx_t> merge_vectors<sidx_t>(const std::vector< std::vector<sidx_t>* >& vec);
  

