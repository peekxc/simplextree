#ifndef DISCRETE_H_
#define DISCRETE_H_

#include "short_alloc.h"    // stack-based allocation helpers
#include <assert.h>         // assertions
#include <array>
#include <cmath>  // round, sqrt, etc. 

using std::floor; 
using std::sqrt;
using std::round; 

template < class T, std::size_t BufSize = 32 >
using SmallVector = std::vector<T, short_alloc<T, BufSize, alignof(T)>>;

// Integral-type to use in the combinadics
#include <cstdint>
using cc_int_t = uint_fast64_t;

// Szudziks pairing function. Takes as input two unsigned integral types (a, b), and uniquely
// maps the pair (a, b) to a distinct number c, where c is possibly a different integral type
template < typename T1 = uint_fast32_t, typename T2 = uint_fast64_t >
constexpr inline T2 szudzik_pair(T1 x, T1 y){
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  return static_cast< T2 >(x >= y ? x * x + x + y : x + y * y);
}
template < typename T1 = uint_fast32_t, typename T2 = uint_fast64_t >
inline std::pair< T1, T1 > szudzik_unpair(T2 z) {
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  T2 sqrtz = std::floor(std::sqrt(z));
  T2 sqz = sqrtz * sqrtz;
  return ((z - sqz) >= sqrtz) ?
    std::make_pair(static_cast< T1 >(sqrtz), static_cast< T1 >(z - sqz - sqrtz)) :
    std::make_pair(static_cast< T1 >(z - sqz), static_cast< T1 >(sqrtz));
}

// constexpr implicitly inlined 
template < bool i_less_j = false >
constexpr auto to_natural_2(cc_int_t i, cc_int_t j, cc_int_t n) noexcept -> 
  typename std::enable_if< i_less_j == true, cc_int_t >::type  {
  return cc_int_t(n*i - i*(i+1)/2 + j - i - 1);
}

template < bool i_less_j = false >
constexpr auto to_natural_2(cc_int_t i, cc_int_t j, cc_int_t n) noexcept -> 
  typename std::enable_if< i_less_j == false, cc_int_t >::type  {
  return i < j ? cc_int_t(n*i - i*(i+1)/2 + j - i - 1) : cc_int_t(n*j - j*(j+1)/2 + i - j - 1);
}

// 0-based
inline std::array< cc_int_t, 2 > to_subscript_2(const cc_int_t x, const cc_int_t n) noexcept {
	auto i = static_cast< cc_int_t >( (n - 2 - floor(sqrt(-8*x + 4*n*(n-1)-7)/2.0 - 0.5)) );
	auto j = static_cast< cc_int_t >( x + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 );
	return (std::array< cc_int_t, 2 >{ i, j });
}

// static constexpr size_t max_choose = 16;
// static constexpr std::array< size_t, 120 > BC = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,3,6,10,15,21,28,36,45,55,66,78,91,105,4,10,20,35,56,84,120,165,220,286,364,455,5,15,35,70,126,210,330,495,715,1001,1365,6,21,56,126,252,462,792,1287,2002,3003,7,28,84,210,462,924,1716,3003,5005,8,36,120,330,792,1716,3432,6435,9,45,165,495,1287,3003,6435,10,55,220,715,2002,5005,11,66,286,1001,3003,12,78,364,1365,13,91,455,14,105,15 };

template < size_t n, size_t k >
constexpr auto bc_recursive() noexcept -> typename std::enable_if< k == 0, size_t >::type {
	return 1;
}
template < size_t n, size_t k >
constexpr auto bc_recursive() noexcept -> typename std::enable_if< (k > 0), size_t >::type {
  return (n * bc_recursive< n - 1, k - 1>()) / k;
}

// template< size_t n, size_t k >
// struct BinomialCoefficientTable {
//   size_t combinations[n+1][k];
//   constexpr BinomialCoefficientTable() : combinations() {
// 		// auto n_dispatcher = make_index_dispatcher< n+1 >();
// 		// auto k_dispatcher = make_index_dispatcher< k >();
// 		// n_dispatcher([&](auto i) { 
// 		// 	k_dispatcher([&](auto j){
// 		// 		combinations[i][j] = bc_recursive< i, j>();
// 		// 	});
// 		// });
// 		combinations[0][0] = bc_recursive< 0, 0 >();
//     combinations[0][1] = bc_recursive< 0, 1 >();
//     combinations[0][2] = bc_recursive< 0, 2 >();
//     combinations[1][0] = bc_recursive< 0, 0 >();
//     combinations[1][1] = bc_recursive< 0, 1 >();
//     combinations[1][2] = bc_recursive< 0, 2 >();
//     combinations[2][0] = bc_recursive< 0, 0 >();
//     combinations[2][1] = bc_recursive< 0, 1 >();
//     combinations[2][2] = bc_recursive< 0, 2 >();
//   }
// };
// static constexpr auto BC = BinomialCoefficientTable< max_choose, max_choose >();

// Recursive binomial coefficient dispatcher; should be compiled to support up to max_choose - 1 at compile time
	// return BinomialCoefficient(n-1, k-1) + BinomialCoefficient(n-1, std::min(n-1-k, k));	
inline size_t binomial_coeff_(const double n, size_t k) noexcept {
  double bc = n;
  for (size_t i = 2; i <= k; ++i){ bc *= (n+1-i)/i; }
  return(static_cast< size_t >(std::round(bc)));
}

static constexpr size_t max_comb = 16; 
static constexpr std::array< size_t, 120 > BC = { 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,3,6,10,15,21,28,36,45,55,66,78,91,105,120,4,10,20,35,56,84,120,165,220,286,364,455,560,5,15,35,70,126,210,330,495,715,1001,1365,1820,6,21,56,126,252,462,792,1287,2002,3003,4368,7,28,84,210,462,924,1716,3003,5005,8008,8,36,120,330,792,1716,3432,6435,11440,9,45,165,495,1287,3003,6435,12870,10,55,220,715,2002,5005,11440,11,66,286,1001,3003,8008,12,78,364,1365,4368,13,91,455,1820,14,105,560,15,120,16 };
  
inline size_t BinomialCoefficient(const size_t n, const size_t k){
  if (k == 0 || n == k){ return 1; }
  if (n < k){ return 0; }
  return n < max_comb ? BC[to_natural_2< true >(k-1,n-1,max_comb)] : binomial_coeff_(n,std::min(k,n-k));
}

// Given a binomial coefficient 'x' representing (n choose 2), finds 'n'
inline size_t inv_choose_2(const size_t x) noexcept {
  const size_t a = floor(sqrt(2*x));
  const size_t b = ceil(sqrt(2*x)+2);
  SmallVector< size_t >::allocator_type::arena_type arena;
  SmallVector< size_t > rng{ arena };
  rng.resize((b - a) + 1);
  std::iota(begin(rng), end(rng), a);
  auto it = std::find_if(begin(rng), end(rng), [x](size_t n){ return(x == BinomialCoefficient(n, 2)); });
  return it == end(rng) ? 0 : *it;
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
  // if constexpr (check_inputs){
  //   if (k > n) { throw std::out_of_range ("Combinadic out of range (k > n)"); }
  //   if (std::distance(s,e) % k != 0){ throw std::out_of_range ("Invalid input; not aligned to k."); }
  //   bool any_overflow = std::any_of(s, e, [n](auto& i){ return(i < 0 || i >= n); });
  //   if (any_overflow){ throw std::out_of_range ("Invalid combinadic input."); } 
  // }
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
  // if constexpr (check_inputs){
  //   if (k > n) { throw std::out_of_range ("Combinadic out of range."); }
  //   const size_t N = BinomialCoefficient(n, k); 
  //   bool any_overflow = std::any_of(s, e, [N](auto& i){ return(i < 0 || i >= N); });
  //   if (any_overflow){ throw std::out_of_range ("Invalid combinadic input."); } 
  // }
  
  SmallVector< cc_int_t >::allocator_type::arena_type a;
  SmallVector< cc_int_t > combination{ a };
  combination.resize(k);
  switch(k){
    case 2:{
      using id_t = typename Iter::value_type; 
      std::array< cc_int_t, 2 > cc; 
      std::for_each(s, e, [n, &cc, &f, &combination](id_t i){
        cc = to_subscript_2(i, n);
        std::move(cc.begin(), cc.end(), combination.begin());
        f(combination);
      });
      break;
    }
    default: {
      using id_t = typename Iter::value_type; 
      const size_t N = BinomialCoefficient(n, k);
      std::for_each(s, e, [&](id_t m){
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
