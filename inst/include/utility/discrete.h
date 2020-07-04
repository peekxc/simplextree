#ifndef DISCRETE_H_
#define DISCRETE_H_

#include "meta_utilities.h" // template-utilities
#include "short_alloc.h"    // stack-based allocation helpers
#include <assert.h>         // assertions
#include <array>

template < class T, std::size_t BufSize = 32 >
using SmallVector = std::vector<T, short_alloc<T, BufSize, alignof(T)>>;

// Integral-type to use in the combinadics
#include <cstdint>
using cc_int_t = uint_fast64_t;

// Szudziks pairing function. Takes as input two unsigned integral types (a, b), and uniquely 
// maps the pair (a, b) to a distinct number c, where c is possibly a different integral type
template < typename T1, typename T2 = T1 > 
constexpr inline T2 szudzik_pair(T1 x, T1 y){
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  T2 a = static_cast< T2 >(x), b = static_cast< T2 >(y);
  return a >= b ? a * a + a + b : a + b * b;
}
template < typename T1, typename T2 = T1 > 
constexpr inline T2 szudzik_unpair(T1 z) {
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  T1 sqrtz = std::floor(std::sqrt(z));
  T1 sqz = sqrtz * sqrtz;
  return ((z - sqz) >= sqrtz) ? 
    std::make_pair(static_cast< T2 >(sqrtz), static_cast< T2 >(z - sqz - sqrtz)) : 
    std::make_pair(static_cast< T2 >(z - sqz), static_cast< T2 >(sqrtz));
}

template < size_t n, size_t k >
constexpr size_t bc_recursive() noexcept {
	if constexpr (k == 0) { return 1; }
 	else { return (n * bc_recursive< n - 1, k - 1>()) / k; }
}

template< size_t n, size_t k >
struct BinomialCoefficientTable {
  size_t combinations[n+1][k];
  constexpr BinomialCoefficientTable() : combinations() {
		auto n_dispatcher = make_index_dispatcher< n+1 >();
		auto k_dispatcher = make_index_dispatcher< k >();
		n_dispatcher([&](auto i) { 
			k_dispatcher([&](auto j){
				combinations[i][j] = bc_recursive< i, j>();
			});
		});
  }
};

static constexpr size_t max_choose = 10;
static constexpr auto BC = BinomialCoefficientTable< max_choose, max_choose >();

// Recursive binomial coefficient dispatcher; should be compiled to support up to max_choose - 1 at compile time
inline size_t BinomialCoefficient(const size_t n, const size_t k) noexcept {
	if (n < k){ return(0); }
	switch(n){
		case 0: return 0;
		case 1: 
			return BC.combinations[n][k];
		case 2: 
			return BC.combinations[n][k];
		case 3: 
			return BC.combinations[n][k];
		case 4: 
			return BC.combinations[n][k];
		case 5: 
			return BC.combinations[n][k];
		case 6: 
			return BC.combinations[n][k];
		case 7: 
			return BC.combinations[n][k];
		case 8: 
			return BC.combinations[n][k];
		case 9: 
			return BC.combinations[n][k];
		default: {
			if (n == k || k == 0){ return 1; }
			return BinomialCoefficient(n-1, k-1) + BinomialCoefficient(n-1, std::min(n-1-k, k));		
		}	
	}
}

// Given a binomial coefficient 'x' representing (n choose 2), finds 'n'
inline size_t inv_choose_2(const size_t x) noexcept {
  const size_t a = floor(sqrt(2*x)), b = ceil(sqrt(2*x)+2);
  SmallVector< size_t >::allocator_type::arena_type arena;
  SmallVector< size_t > rng{ arena };
  rng.resize((b - a) + 1);
  std::iota(begin(rng), end(rng), a);
  auto it = std::find_if(begin(rng), end(rng), [x](size_t n){ return(x == BinomialCoefficient(n, 2)); });
  return it == end(rng) ? 0 : *it;
}

// constexpr implicitly inlined 
constexpr cc_int_t to_natural_2(cc_int_t i, cc_int_t j, cc_int_t n) noexcept {
	return i < j ? cc_int_t(n*i - i*(i+1)/2 + j - i - 1) : cc_int_t(n*j - j*(j+1)/2 + i - j - 1);
}

// 0-based
inline std::array< cc_int_t, 2 > to_subscript_2(const cc_int_t x, const cc_int_t n) noexcept {
	auto i = static_cast< cc_int_t >( (n - 2 - floor(sqrt(-8*x + 4*n*(n-1)-7)/2.0 - 0.5)) );
	auto j = static_cast< cc_int_t >( x + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 );
	return { i, j };
}

// converts to natural number
template < typename Iter > 
inline size_t to_natural_k(Iter s, Iter e, const cc_int_t n, const cc_int_t k) {
  if (n == k){ return(0); }
  
	// Apply the dual index mapping
	const cc_int_t N = BinomialCoefficient(n,k); 
	cc_int_t i = k; 
  const cc_int_t index = std::accumulate(s, e, 0, [n, &i](cc_int_t val, cc_int_t num){ 
	  return val + BinomialCoefficient((n-1) - num, i--); 
	});
  const cc_int_t combinadic = (N-1) - index;
  return(combinadic);
}
// if (combinadic >= N || index > (N-1)){
//   // throw std::out_of_range ("Combinadic mapping failed.");
//   // std::cout << "n: " << n << ", k: " << k << ", N: " << N << ", combinadic: " << combinadic << ", index: " << index << std::endl;
//   // std::for_each(s,e, [](auto el){ std::cout << el << ","; });
//   // std::cout << std::endl;
// }

// 0-based conversion of (n choose k) combinadic subscripts to natural number
template < bool check_inputs = true, typename Iter, typename Lambda > 
inline void to_natural(Iter s, Iter e, const size_t n, const size_t k, Lambda&& f) {
  if constexpr (check_inputs){
    if (k > n) { throw std::out_of_range ("Combinadic out of range (k > n)"); }
    if (std::distance(s,e) % k != 0){ throw std::out_of_range ("Invalid input; not aligned to k."); }
    bool any_overflow = std::any_of(s, e, [n](auto& i){ return(i < 0 || i >= n); });
    if (any_overflow){ throw std::out_of_range ("Invalid combinadic input."); } 
  }
  if (n == k){ while (s != e){ f(0); s += k; } }
  while (s != e){
    switch(k){
      case 2: 
        f(to_natural_2(*s, *std::next(s), n));
        break;
      default: 
        f(to_natural_k(s, s+k, n, k));
        break;
    }
    s += k;
  }
}

// Converts each value between [s,e) to its corresponding combinadic
template< bool check_inputs = true, typename Iter, typename Lambda >
inline void to_subscript(Iter s, Iter e, const size_t n, const size_t k, Lambda&& f) {
  
  // If check inputs true (default), make sure the input is valid. 
  if constexpr (check_inputs){
    if (k > n) { throw std::out_of_range ("Combinadic out of range."); }
    const size_t N = BinomialCoefficient(n, k); 
    bool any_overflow = std::any_of(s, e, [N](auto& i){ return(i < 0 || i >= N); });
    if (any_overflow){ throw std::out_of_range ("Invalid combinadic input."); } 
  }
  
  SmallVector< cc_int_t >::allocator_type::arena_type a;
  SmallVector< cc_int_t > combination{ a };
  combination.resize(k);
  switch(k){
    case 2:{
      std::array< cc_int_t, 2 > cc; 
      std::for_each(s, e, [n, &cc, &f, &combination](auto& i){
        cc = to_subscript_2(i, n);
        std::move(cc.begin(), cc.end(), combination.begin());
        f(combination);
      });
      break;
    }
    default: {
      const size_t N = BinomialCoefficient(n, k);
      std::for_each(s, e, [&](auto m){
        m = ((N-1)-m); 
        auto guess = 0; 
        auto pc = n;
        for (auto j = k; j > 0; --j){
        	cc_int_t value = m + 1;
          for (; value > m; pc = guess, --guess){
            guess = pc - 1; 
            value = BinomialCoefficient(guess, j);
          }
        	m = m - value;
        	combination[k-j] = (n-1)-guess-1;
      	}
        f(combination);
      });
      break; 
    }
  }
}

#endif 