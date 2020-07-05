// Type traits and other meta-programming necessities for the simplex tree implementation
#ifndef META_UTILITIES_H_
#define META_UTILITIES_H_

#include <type_traits>
#include <memory> 

  // detail namespace contains internal template boilerplate, including several implementations of the detection idiom 
namespace detail {
//   template <typename T> using s_ptr = std::shared_ptr<T>; // Shared pointer
//   template <typename T> using u_ptr = std::unique_ptr<T>; // Unique pointer
//   template <typename T> using w_ptr = std::weak_ptr<T>;   // Weak pointer
//   
//   // Detection Idiom
//   template <template <class...> class Trait, class Enabler, class... Args>
//   struct is_detected_trait : std::false_type{};
// 
//   template <template <class...> class Trait, class... Args>
//   struct is_detected_trait<Trait, std::void_t<Trait<Args...>>, Args...> : std::true_type{};
// 
//   template <template <class...> class Trait, class... Args>
//   using is_detected = typename is_detected_trait<Trait, void, Args...>::type;
//   
// 	template <typename T>
// 	class remove_pointer_ {
// 			template <typename U=T>
// 			static auto test(int) -> std::remove_reference<decltype(*std::declval<U>())>;
// 			static auto test(...) -> std::remove_cv<T>;
// 		public:
// 			using type = typename decltype(test(0))::type;
// 	};
//   
  // template <typename S, typename T> 
  // using similar = is_same< decay_t<S>, decay_t<T> >;
  // 
  // template <typename T>
  // using dereference_t = decltype(std::declval<T>().operator*());
  // 
  // template <typename T> 
  // using is_dereferenceable = is_detected< dereference_t, T>;
  // 
  // template <typename T> 
  // using is_ptr = disjunction< is_pointer<T>,  is_dereferenceable< T > >; 
  // 
  // template <typename T>
  // inline constexpr bool is_ptr_v = is_ptr< T >::value;
  // 
  // template <typename T, typename U = std::decay_t<T> >
  // struct less_ptr {
  //   bool operator() (const T& lhs, const T& rhs) const { 
  // 		if constexpr (is_ptr_v< U >){ return (*lhs) < (*rhs); } 
  // 		else { return lhs < rhs; }
  // 	}
  // };
  // template <typename T>
  // struct has_dereference {
  //   template <typename U> 
  // 	static constexpr auto test_dereference(int) -> decltype(std::declval<U>().operator*(), bool()) { 
  // 		return true;
  //   }
  //   template <typename U>
  //   static constexpr bool test_dereference(...) {
  //     return false;
  //   }
  //   static constexpr bool value = test_dereference<T>(int());
  // };
}; // end namespace detail

// Compile-time index dispatchers from amazing post: https://blog.tartanllama.xyz/exploding-tuples-fold-expressions/
// template <std::size_t... Idx>
// constexpr auto make_index_dispatcher(std::index_sequence<Idx...>) {
// 	return [] (auto&& f) { (f(std::integral_constant<std::size_t,Idx>{}), ...); };
// };
// 
// template <std::size_t N>
// constexpr auto make_index_dispatcher() {
// 	return make_index_dispatcher(std::make_index_sequence< N >{}); 
// };
// 

#endif 