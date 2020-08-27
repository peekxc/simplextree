#ifndef COMBINATIONS_H
#define COMBINATIONS_H

//  (C) Copyright Howard Hinnant 2005-2011.
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).
//
//  See http://www.boost.org/libs/type_traits for most recent version including documentation.

// This code was adapted from Howard Hinnants excellent combinations header:
// https://github.com/HowardHinnant/combinations
// The original copyright is included above. 

#include <iterator>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace detail {
	template<typename T>
	using it_diff_t = typename std::iterator_traits<T>::difference_type;

	// Rotates two discontinuous ranges to put *first2 where *first1 is.
	// Adapted from: https://github.com/HowardHinnant/combinations
	template <class It>
	void rotate_discontinuous(
		It first1, It last1, it_diff_t< It > d1,
	  It first2, It last2, it_diff_t< It > d2) 
	{
		using std::swap;
		if (d1 <= d2){ std::rotate(first2, std::swap_ranges(first1, last1, first2), last2); }
		else {
			It i1 = last1;
			while (first2 != last2){ swap(*--i1, *--last2); }
			std::rotate(first1, i1, last1);
		}
	}

	// Call f() for each combination of the elements [first1, last1) + [first2, last2)
	// swapped/rotated into the range [first1, last1).  
	template < typename It, typename Lambda >
	bool combine_discontinuous(
		It first1, It last1, it_diff_t< It > d1,  
	 	It first2, It last2, it_diff_t< It > d2,
		Lambda&& f, it_diff_t< It > d = 0)
	{
		using D = it_diff_t< It >;
		using std::swap;
		if (d1 == 0 || d2 == 0){ return f(); }
		if (d1 == 1) {
			for (It i2 = first2; i2 != last2; ++i2) {
				if (f()){ return true; }
				swap(*first1, *i2);
			}
		}
		else {
			It f1p = std::next(first1), i2 = first2;
			for (D d22 = d2; i2 != last2; ++i2, --d22){
				if (combine_discontinuous(f1p, last1, d1-1, i2, last2, d22, f, d+1))
					return true;
				swap(*first1, *i2);
			}
		}
		if (f()){ return true; }
		if (d != 0){ rotate_discontinuous(first1, last1, d1, std::next(first2), last2, d2-1); }
		else { rotate_discontinuous(first1, last1, d1, first2, last2, d2); }
		return false;
	}
	
	template < typename Lambda, class... Ts >
	struct NullaryPredicate {
		Lambda f_;
		std::tuple< Ts... > params;
		NullaryPredicate(Lambda& f, Ts ... args) : f_(f), params(std::make_tuple(std::move(args)...)) {};
		bool operator()() { return f_(params); }
	};
}; // namespace detail
	
using namespace detail; 	
	
template < class It, class Function >
Function for_each_combination(It first, It mid, It last, Function&& f) {
	combine_discontinuous(first, mid, std::distance(first, mid),
												mid, last, std::distance(mid, last),
												[&f, &first, &mid]() -> bool { return f(first, mid); });
	return std::move(f);
}

template < class I, class Function >
void for_each_combination_idx(I n, I k, Function&& f) {
  static_assert(std::is_integral<I>::value, "Must be integral type.");
  using It = typename std::vector< I >::iterator;
  std::vector< I > seq_n(n);
  std::iota(begin(seq_n), end(seq_n), 0);
  for_each_combination(begin(seq_n), begin(seq_n)+k, end(seq_n), [&f](It a, It b) -> bool {
    std::vector< I > cc(a, b);
    f(cc);
    return false; 
  });
	return; 
}

#endif  // COMBINATIONS_H

