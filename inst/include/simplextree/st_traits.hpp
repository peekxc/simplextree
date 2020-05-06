// st_traits.h
// Type traits and other meta-programming necessities for the simplex tree implementation
#ifndef ST_TRAITS_H
#define ST_TRAITS_H

#include <type_traits>
#include <memory> 

using std::is_pointer;
using std::is_pointer_v;
using std::is_same; 
using std::is_same_v;
using std::is_base_of; 
using std::enable_if;
using std::decay; 
using std::decay_t;
using std::conjunction; 
using std::conjunction_v; 
using std::disjunction;
using std::disjunction_v;
using std::remove_pointer;

namespace detail {
  template <template <class...> class Trait, class Enabler, class... Args>
  struct is_detected : std::false_type{};

  template <template <class...> class Trait, class... Args>
  struct is_detected<Trait, std::void_t<Trait<Args...>>, Args...> : std::true_type{};

	template <typename T>
	class remove_pointer_ {
			template <typename U=T>
			static auto test(int) -> std::remove_reference<decltype(*std::declval<U>())>;
			static auto test(...) -> std::remove_cv<T>;
		public:
			using type = typename decltype(test(0))::type;
	};
}

template <typename T> using s_ptr = std::shared_ptr<T>; // Shared pointer
template <typename T> using u_ptr = std::unique_ptr<T>; // Unique pointer
template <typename T> using w_ptr = std::weak_ptr<T>;   // Weak pointer

template <typename S, typename T> 
using similar = is_same< decay_t<S>, decay_t<T> >;

template <template <class...> class Trait, class... Args>
using is_detected = typename detail::is_detected<Trait, void, Args...>::type;

template <typename T>
using dereference_t = decltype(std::declval<T>().operator*());

template <typename T> 
using is_dereferenceable = is_detected< dereference_t, T>;

template <typename T>
using remove_pointer_t = typename detail::remove_pointer_<T>::type;

template <typename T> 
using is_ptr = disjunction< is_pointer<T>,  is_dereferenceable<T> >; 

template <typename T>
inline constexpr bool is_ptr_v = is_ptr< T >::value;

template <typename T, typename U = std::decay_t<T> >
struct less_ptr {
  bool operator() (const T& lhs, const T& rhs) const { 
		if constexpr (is_ptr_v< U >){ return (*lhs) < (*rhs); } 
		else { return lhs < rhs; }
	}
};

// So good: https://blog.tartanllama.xyz/exploding-tuples-fold-expressions/
template <std::size_t... Idx>
auto make_index_dispatcher(std::index_sequence<Idx...>) {
	return [] (auto&& f) { (f(std::integral_constant<std::size_t,Idx>{}), ...); };
};

template <std::size_t N>
auto make_index_dispatcher() {
	return make_index_dispatcher(std::make_index_sequence<N>{}); 
};

// Szudziks pairing function. Takes as input two unsigned integral types (a, b), and uniquely 
// maps (a, b) to a number c, where c is possibly a different integral type 
template <typename T1, typename T2> 
constexpr inline T2 szudzik_pair(T1 a, T1 b){
  static_assert(std::is_integral<T1>::value, "Integral-type required as a range storage type.");
  static_assert(std::is_unsigned<T1>::value, "Integral-type required as a range storage type.");
  return static_cast<T2>(a >= b ? a * a + a + b : a + b * b);
}
  

#endif 