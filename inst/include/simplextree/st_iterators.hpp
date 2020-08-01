#ifndef ST_ITERS_H
#define ST_ITERS_H

#include "simplextree.h"
#include <stack>
#include <queue>
#include <functional>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <iostream>
#include <type_traits>

using std::get; 
using simplex_t = SimplexTree::simplex_t;
using node_ptr = SimplexTree::node_ptr;
using node_uptr = SimplexTree::node_uptr;


// #include <string_view>
// 
// template <typename T>
// constexpr std::string_view 
// type_name() {
//     std::string_view name, prefix, suffix;
// #ifdef __clang__
//     name = __PRETTY_FUNCTION__;
//     prefix = "std::string_view type_name() [T = ";
//     suffix = "]";
// #elif defined(__GNUC__)
//     name = __PRETTY_FUNCTION__;
//     prefix = "constexpr std::string_view type_name() [with T = ";
//     suffix = "; std::string_view = std::basic_string_view<char>]";
// #elif defined(_MSC_VER)
//     name = __FUNCSIG__;
//     prefix = "class std::basic_string_view<char,struct std::char_traits<char> > __cdecl type_name<";
//     suffix = ">(void)";
// #endif
//     name.remove_prefix(prefix.size());
//     name.remove_suffix(suffix.size());
//     return name;
// }

// Simplextree namespace
namespace st {

// detail namespace contains internal template boilerplate, including several implementations of the detection idiom 
  namespace detail {
    template <typename T>
    struct has_dereference {
      template <typename U> 
    	static constexpr auto test_dereference(int) -> decltype(std::declval<U>().operator*(), bool()) { 
    		return true;
      }
      template <typename U>
      static constexpr bool test_dereference(...) {
        return false;
      }
      static constexpr bool value = test_dereference<T>(int());
    };
    
    // Detection idom applied to member function operator==
    template <typename T>
    struct has_equality {
      template <typename U> 
    	static constexpr auto test_equality(int) -> decltype(std::declval<U>().operator==(), bool()) { 
    		return true;
      }
      template <typename U>
      static constexpr bool test_equality(...) {
        return false;
      }
      static constexpr bool value = test_equality<T>(int());
    };
    
    // Detection idiom applied to member function operator!=
    template <typename T>
    struct has_not_equality {
      template <typename U> 
    	static constexpr auto test_not_equality(int) -> decltype(std::declval<U>().operator!=(), bool()) { 
    		return true;
      }
      template <typename U>
      static constexpr bool test_not_equality(...) {
        return false;
      }
      static constexpr bool value = test_not_equality<T>(int());
    };
    
    // Detection idiom applied to member function operator++
    template <typename T>
    struct has_increment {
      template <typename U> 
    	static constexpr auto test_increment(int) -> decltype(std::declval<U>().operator++(), bool()) { 
    		return true;
      }
      template <typename U>
      static constexpr bool test_increment(...) {
        return false;
      }
      static constexpr bool value = test_increment<T>(int());
    };
    
    template <typename T>
    struct has_begin {
      template <typename U> 
    	static constexpr auto test_begin(int) -> decltype(std::declval<U>().begin(), bool()) { return true; }
      template <typename U>
      static constexpr bool test_begin(...) { return false; }
      static constexpr bool value = test_begin<T>(int());
    };
    
    template <typename T>
    struct has_end {
      template <typename U> 
    	static constexpr auto test_end(int) -> decltype(std::declval<U>().end(), bool()) { return true; }
      template <typename U>
      static constexpr bool test_end(...) { return false; }
      static constexpr bool value = test_end<T>(int());
    };
    
    template <typename T>
    struct has_update_simplex {
      template <typename U> 
    	static constexpr auto test_update_simplex(int) -> decltype(std::declval<U>().has_update_simplex(), bool()) { return true; }
      template <typename U>
      static constexpr bool test_update_simplex(...) { return false; }
      static constexpr bool value = test_update_simplex<T>(int());
    };
  }; // end namespace detail

  template < bool ts, template< bool > class Derived > 
  struct TraversalInterface {
  	using d_type = Derived< ts >;
  	using d_node = std::tuple< node_ptr, idx_t >;
  	using t_node = typename std::conditional< ts, std::tuple< node_ptr, idx_t, simplex_t >, d_node >::type;
  	using pair_pred_t = delegate< bool (t_node&) >;
  
  	public: 
  	  static const bool is_tracking = ts;
  		static const idx_t NP = 0;
    	static const idx_t DEPTH = 1;
  		static const idx_t LABELS = 2;
  		node_ptr init;
  		const SimplexTree* st;
  		pair_pred_t p1 = [](t_node& cn){ return true; };
  		pair_pred_t p2 = [](t_node& cn){ return true; };
  
  		TraversalInterface() = default;
  		TraversalInterface(const SimplexTree* st_) : st(st_){
  			init = nullptr;
  		}
  		TraversalInterface(const SimplexTree* st_, node_ptr start) : st(st_){
  			init = start;
  		};
  		template< typename P1, typename P2 >
  		TraversalInterface(const SimplexTree* st_, node_ptr start, P1 pred1, P2 pred2) : init(start), st(st_), p1(pred1), p2(pred2){
  			init = start;
  		};
  
  		struct iterator {
  			using difference_type = std::ptrdiff_t;
  			using value_type = t_node; 
  			using pointer = t_node*; 
  			using reference = t_node&;
  			using iterator_category = std::forward_iterator_tag;
  			
  			static constexpr d_node sentinel() { return std::make_tuple(nullptr, 0); };
  
  			using d_iter = typename d_type::iterator;
  			std::reference_wrapper< d_type > info;
  			d_node current; 
  			simplex_t labels; 
  			t_node output; 
  
  			// Constructors 
  			iterator(d_type& dd) : info(dd){
					labels = simplex_t();
					labels.reserve(dd.st->tree_max_depth);
  			};
  
  			constexpr d_type& base() const { return(this->info.get()); };
  			constexpr const SimplexTree& trie() const { return(*(this->info.get().st)); };
  			
  			template < bool T = ts >
        typename std::enable_if< T == true, t_node& >::type
  			current_t_node(){
  			  output = std::tuple_cat(current, std::make_tuple(labels)); 
  				return(output);
  			}
  			template < bool T = ts>
        typename std::enable_if< T == false, t_node& >::type
  			current_t_node(){ return(current); }
  			
  			// Operators
  			bool operator==(const d_iter& t) const { return (get< 0 >(t.current) == get< 0 >(current)); }
  			bool operator!=(const d_iter& t) const { return !(*this == t); }
  			t_node& operator*() { return current_t_node(); };
  		};
  		
  		// methods within TraversalInterface can use template to access members of Derived
  		// auto begin() { return static_cast< d_type* >(this)->begin(); };
  		// auto end() { return static_cast< d_type* >(this)->end(); };
  		
  		static constexpr bool derived_has_begin = detail::has_begin< Derived< ts > >::value;
  	  static constexpr bool derived_has_end = detail::has_end< Derived< ts > >::value;
  		
  	// 	template < typename std::enable_if< true, int >::type = 0 >
  	// 	auto begin(){
  	// 		return static_cast< d_type* >(this)->begin();
  	// 	};
  	//   template < typename std::enable_if< false, int >::type = 0 >
  	// 	auto begin(){
  	// 		return static_cast< d_type* >(this)->begin();
  	// 	};
  		  		
  		// template < typename std::enable_if< !detail::has_begin< d_type >::value, int >::type = 0 >
  		// auto begin() {
  		// 	return static_cast< d_type* >(this)->begin();
  		// };
  		// -> std::enable_if< false, decltype(static_cast< d_type* >(this)->begin()) >
  		// 
  		// auto end() -> std::enable_if< detail::has_begin< d_type >::value, decltype(d_type::iterator(static_cast< d_type& >(*this), nullptr)) >{
  		// 	return d_type::iterator(static_cast< d_type& >(*this), nullptr);
  		// };
  		// auto end() -> std::enable_if< !detail::has_begin< d_type >::value, decltype(static_cast< d_type* >(this)->begin()) >{
  		// 	return static_cast< d_type* >(this)->end();
  		// };
  		
  		// Tag dispatching
  		// using d_iter = typename Derived< ts >::iterator;
  		// template<typename T> struct TD;
  		
  // 		template<typename T> struct TD;
  // 		
  //     auto begin(std::true_type) {
  //       return static_cast< d_type* >(this)->begin();
  //     };
  //     
  //     auto begin(std::false_type) {
  //       using d_iter = typename Derived< ts >::iterator;
  // 			return d_iter(static_cast< d_type& >(*this), nullptr);
  //     };
  //     
  //     auto begin() {
  //       auto b = begin(derived_has_begin); 
  //       std::cout << type_name<decltype(b)>() << std::endl;
  //       return b; 
  //     };
  // 		
  // 		// Tag dispatching
  //     auto end(std::true_type) {
  //       return static_cast< d_type* >(this)->end();
  //     };
  //     
  //     auto end(std::false_type) {
  //       using d_iter = typename d_type::iterator;
  // 			return d_iter(static_cast< d_type& >(*this), nullptr);
  //     };
  //     
  //     auto end() {
  //       return end(derived_has_end); 
  //     };

  		// auto end() {
  		// 	if constexpr(derived_has_end){
  		// 		return static_cast< d_type* >(this)->end();
  		// 	} else {
  		// 		using d_iter = typename d_type::iterator;
  		// 		return d_iter(static_cast< d_type& >(*this), nullptr);
  		// 	}
  		// };
  }; // end Traversal Interface
  
  // https://eli.thegreenplace.net/2011/05/17/the-curiously-recurring-template-pattern-in-c
  
  template < class T, std::size_t BufSize = 64 >
  using small_vec = std::vector<T, short_alloc<T, BufSize, alignof(T)>>;
    		
  // Preorder traversal iterator
  template < bool ts = false > 
  struct preorder : TraversalInterface< ts, preorder > {
  	using B = TraversalInterface< ts, preorder >;
  	using B::init;
  	using B::st; 
  	using B::p1;
  	using B::p2;
  	using B::NP;
  	using B::DEPTH;
  	using B::LABELS;
  
  	// Constructors
  	preorder(const SimplexTree* st_) : TraversalInterface<  ts, preorder >(st_, st_->root.get()) {};
  	preorder(const SimplexTree* st_, node_ptr start) : TraversalInterface<  ts, preorder >(st_, start) {};
  	template< typename P1, typename P2 >
  	preorder(const SimplexTree* st_, node_ptr start, P1 pred1, P2 pred2) : TraversalInterface<  ts, preorder >(st_, start, pred1, pred2){};
  
  	void reset(node_ptr start){
  		init = start; 
  	}
  
  	// Iterator type
  	struct iterator : public TraversalInterface< ts, preorder >::iterator {
  		using Bit = typename B::iterator;
  		using Bit::current; 
  		using Bit::labels;
  		using Bit::info;
  		using Bit::base;
  		using Bit::trie;
  		using Bit::current_t_node;
  		using Bit::sentinel;
  
  		// DFS specific data structures
  		using dn_t = typename B::d_node;
  		std::stack< typename B::d_node > node_stack;
  		
    	// 	small_vec< dn_t >::allocator_type::arena_type a;
    	//   small_vec< dn_t > v = {a};
    	// 	std::stack< dn_t, small_vec< dn_t > > node_stack = {v};

  		// Iterator constructor
  		iterator(preorder& dd, node_ptr cn = nullptr) : TraversalInterface< ts, preorder >::iterator(dd){
  			const idx_t d = dd.st->depth(cn);
  			current = std::make_tuple(cn, d);
  		  labels = dd.st->full_simplex(cn, d);
  		}
  		
  		// Increment operator is all that is needed by default
  		auto operator++() -> decltype(*this){
  			do {
  				if (get< NP >(current) != nullptr && base().p2(current_t_node())) {
  					const auto& ch = get< NP >(current)->children;
  					for (auto cn_up = ch.rbegin(); cn_up != ch.rend(); ++cn_up){
  						node_stack.push(std::make_tuple((*cn_up).get(), get< DEPTH >(current)+1)); 
  					}
  				}
  				if (node_stack.empty()){ 
  					current = sentinel();
  				} else { 
  					current = node_stack.top(); 
  					node_stack.pop(); 
  				}
  				update_simplex();
  			} while (!base().p1(current_t_node()) && get< NP >(current) != nullptr);
  			return(*this);
  		}
  		
  		template < bool T = ts >
      auto update_simplex() -> typename std::enable_if< T == true, void >::type {
        // std::cout << "update_simplex: " << base().init << std::endl;
    	  if (get< NP >(current) != nullptr && get< DEPTH >(current) > 0){
  				labels.resize(get< DEPTH >(current));
  				labels.at(get< DEPTH >(current)-1) = get< NP >(current)->label;
  			}
    	}
  		template < bool T = ts >
      auto update_simplex() -> typename std::enable_if< T == false, void >::type { 
        return;
    	}
  	};
  
  	// Only start at a non-root node, else return the end
  	auto begin() -> decltype(iterator(*this, init)) {
  		if (init == st->root.get()){
  			return st->n_simplexes.empty() ? iterator(*this, nullptr) : ++iterator(*this, st->root.get());
  		} else {
  			return iterator(*this, init);
  		}
  	};
  	auto end() -> decltype(iterator(*this, nullptr)){
  	  return iterator(*this, nullptr);
    }
  	
  };
  
  // level order traversal iterator
  template < bool ts = false > 
  struct level_order : TraversalInterface< ts, level_order > {
  	using B = TraversalInterface< ts, level_order >;
  	using B::init;
  	using B::st; 
  	using B::p1;
  	using B::p2;
  	using B::NP;
  	using B::DEPTH;
  	using B::LABELS;
  	using d_node = typename TraversalInterface< ts, level_order >::d_node;
  
  	// Constructors
  	level_order(const SimplexTree* st_) : TraversalInterface<  ts, level_order >(st_, st_->root.get()) {};
  	level_order(const SimplexTree* st_, node_ptr start) : TraversalInterface<  ts, level_order >(st_, start) {};
  	template< typename P1, typename P2 >
  	level_order(const SimplexTree* st_, node_ptr start, P1 pred1, P2 pred2) : TraversalInterface<  ts, level_order >(st_, start, pred1, pred2){};
  
  	void reset(node_ptr start){ init = start;  }
  
  	// Iterator type
  	struct iterator : public TraversalInterface< ts, level_order >::iterator {
  		using Bit = typename B::iterator;
  		using Bit::current; 
  		using Bit::labels;
  		using Bit::info;
  		using Bit::base;
  		using Bit::trie;
  		using Bit::current_t_node;
  		using Bit::sentinel;
  
  		// BFS specific data structures
  		std::queue< typename B::d_node > node_queue;  
  
  		// Iterator constructor
  		iterator(level_order& dd, node_ptr cn = nullptr) : TraversalInterface< ts, level_order >::iterator(dd){
  			current = std::make_tuple(cn, dd.st->depth(cn));
  		  update_simplex();
  		}
  		
  		// Increment operator is all that is needed by default
  		auto operator++() -> decltype(*this){
  			do {
  				if (get< NP >(current) != nullptr && base().p2(current_t_node())) {
  					const auto& ch = get< NP >(current)->children;
  					for (auto cn_up = ch.begin(); cn_up != ch.end(); ++cn_up){
  						node_queue.emplace(std::make_tuple((*cn_up).get(), get< DEPTH >(current)+1)); 
  					}
  				}
  				if (node_queue.empty()){ 
  					current = sentinel();
  				} else { 
  					current = node_queue.front(); 
  					node_queue.pop(); 
  				}
  				update_simplex();
  			} while (!base().p1(current_t_node()) && get< NP >(current) != nullptr);
  			return(*this);
  		}
  		
  		template < bool T = ts >
      typename std::enable_if< T == true, void >::type
    	update_simplex(){ labels = trie().full_simplex(get< NP >(current), get< DEPTH >(current)); }
  		
    	template < bool T = ts >
      typename std::enable_if< T == false, void >::type
    	update_simplex(){ return; };
  	};// iterator 
  
  	// Only start at a non-root node, else return the end
  	auto begin() -> decltype(iterator(*this, init)){
  		if (init == st->root.get()){
  			return st->n_simplexes.empty() ? iterator(*this, nullptr) : ++iterator(*this, st->root.get()); 
  		} else {
  			return iterator(*this, init);
  		}
  	};
  	auto end() -> decltype(iterator(*this, nullptr)) { return iterator(*this, nullptr); }
  };
  
  
  // Coface-roots search
  template < bool ts = false >
  struct coface_roots : TraversalInterface< ts, coface_roots > {
  	using B = TraversalInterface< ts, coface_roots >;
  	using B::init;
  	using B::st; 
  	using B::p1;
  	using B::p2;
  	using B::NP;
  	using B::DEPTH;
  	using B::LABELS;
  	using d_node = typename TraversalInterface< ts, coface_roots >::d_node;
  
  	// Constructors
  	coface_roots() : TraversalInterface< ts, coface_roots >() {};
  	coface_roots(const SimplexTree* st_, node_ptr start = nullptr) : TraversalInterface<  ts, coface_roots >(st_, start) {};
  	template< typename P1, typename P2 >
  	coface_roots(const SimplexTree* st_, node_ptr start, P1 pred1, P2 pred2) : TraversalInterface<  ts, coface_roots >(st_, start, pred1, pred2){};
  
  	// Iterator type
  	struct iterator : public TraversalInterface< ts, coface_roots >::iterator {
  		using Bit = typename B::iterator;
  		using Bit::current; 
  		using Bit::labels;
  		using Bit::info;
  		using Bit::base;
  		using Bit::trie;
  		using Bit::current_t_node;
  		using Bit::sentinel;
  
  		// coface-roots-specific structures 
  		simplex_t start_coface_s; 
  		size_t c_level_key = 0; // the current level_map key
  		size_t c_level_idx = 0; // the current level_map index 
  
  		// Iterator constructor
  		iterator(coface_roots& dd, node_ptr cn) : TraversalInterface< ts, coface_roots >::iterator(dd){
  			if (cn == dd.st->root.get()){ throw std::invalid_argument("Invalid given coface."); };
  			const size_t c_depth = dd.st->depth(cn); 
  			start_coface_s = dd.st->full_simplex(cn, c_depth);
  			current = std::make_tuple(cn, c_depth);
  			update_simplex();
  			++get< DEPTH >(current); // start at next depth
  		}
  		
  		// Finds the next coface of a given face starting at some offset + key, or nullptr if there is none
  		std::pair< node_ptr, bool > next_coface(simplex_t face, size_t offset, idx_t depth){
  		  auto& st = trie();
  		  bool has_cousins = st.cousins_exist(base().init->label, depth);
  		  
  		  // If the key doesn't exist or the cousins are empty, return the end
  		  if (!has_cousins || offset >= st.cousins(base().init->label, depth).size()){
  		    return std::make_pair(nullptr, false);
  		  }
  		  
  		  // Else the cousins exist; see if any of the are a coface
  		  const auto& c_cousins = st.cousins(base().init->label, depth); 
  		  auto it = c_cousins.cbegin(); //+ offset
  		  std::advance(it, offset);
  		  const auto coface_it = std::find_if(it, c_cousins.cend(), [&st, &face, depth](const node_ptr np){
  		    return st.is_face(face, st.full_simplex(np, depth));
  		  });
  		  
  		  // If it exists, return it, otherwise return sentinel
  		  return coface_it != c_cousins.end() ? std::make_pair(*coface_it, true) : std::make_pair(nullptr, false);
  		}
  		
  		// Increment operator is all that is needed by default
  		auto operator++() -> decltype(*this){
  			// If root was given, end the traversal immediately.
  			if (get< NP >(current) == trie().root.get() || get< NP >(current) == nullptr){ 
  				current = sentinel();
  				return(*this);
  			}
  			
			  // While the coface doesn't exist at the current depth, increment depth
			  std::pair< node_ptr, bool > coface = next_coface(start_coface_s, c_level_idx, get< DEPTH >(current));
			  while (coface.second == false && get< DEPTH >(current) <= trie().tree_max_depth){
			    c_level_idx = 0; 
			    get< DEPTH >(current)++;
			    coface = next_coface(start_coface_s, c_level_idx, get< DEPTH >(current));
			  }
			  
			  // If it doesn't exist, we're done
			  if (coface.second == false){
			    current = sentinel();
			  } else {
			    get< NP >(current) = coface.first;
			    c_level_idx++;
			  }
			  update_simplex();
			  return *this; 
  		}; // operator++
  
  		template < bool T = ts >
      auto update_simplex() -> typename std::enable_if< T == true, void >::type { 
        labels = trie().full_simplex(get< NP >(current), get< DEPTH >(current)); 
      }
  		
    	template < bool T = ts >
      auto update_simplex() -> typename std::enable_if< T == false, void >::type { 
        return; 
      };
    	
  	}; // iterator 
  
  	// Only start at a non-root node, else return the end
  	auto begin() -> decltype(iterator(*this, nullptr)) {
  		if (init == st->root.get() || init == nullptr){ return iterator(*this, nullptr); } 
  		else { return iterator(*this, init); }
  	};
  	auto end() -> decltype(iterator(*this, nullptr)) { return iterator(*this, nullptr); };
  
  }; // end coface_roots 
  
  // ---- coface iterator ------
  template < bool ts = false > 
  struct cofaces : TraversalInterface< ts, cofaces > {
  	using B = TraversalInterface< ts, cofaces >;
  	using B::init;
  	using B::st; 
  	using B::p1;
  	using B::p2;
  	using B::NP;
  	using B::DEPTH;
  	using B::LABELS;
  	using d_node = typename TraversalInterface< ts, cofaces >::d_node;
  
  	cofaces(const SimplexTree* st, node_ptr start) : TraversalInterface< ts, cofaces >(st, start){ }
  
  	struct iterator : public TraversalInterface< ts, cofaces >::iterator {
  		using Bit = typename B::iterator;
  		using Bit::current; 
  		using Bit::labels;
  		using Bit::info;
  		using Bit::base;
  		using Bit::trie;
  		using Bit::current_t_node;
  		using Bit::sentinel;
  
  		using preorder_it = decltype(std::declval< preorder< ts > >().begin());
  		using coface_root_it = decltype(std::declval< coface_roots< false > >().begin());
  
  		// coface-roots-specific structures 
  		coface_roots< false > roots; 
  		coface_root_it c_root;
  		preorder< ts > subtree; 
  		preorder_it c_node; 
  
  		// Iterator constructor
  		iterator(cofaces& dd, node_ptr cn) : TraversalInterface< ts, cofaces >::iterator(dd),
  		 roots(coface_roots< false >(dd.st, cn)), c_root(roots, cn), subtree(preorder< ts >(dd.st, cn)), c_node(subtree.begin()) {
  			 current = std::make_tuple(cn, dd.st->depth(cn));
  		   update_simplex();
  		}
  		 
  		// Increment operator is all that is needed by default
  		auto operator++() -> decltype(*this) {
  			// Need to increment iterator by one and return the result 
  			// Logically what needs to happen is: 
  			// - if next node is end of subtree 
  			//   - ... then if next root is end of coface_roots, we're finished return end 
  			//   - else next root is not end of coface_roots, advance root and reset subtree
  			// else while next node is not at end of subtree 
  			// - advance one in subtree, report it back
  			if (get< NP >(*c_root) == trie().root.get()){ ++c_root; }
  			if (std::next(c_node) == subtree.end()){
  				if (c_root == roots.end()){
  					current = sentinel();
  				} else {
  					++c_root;
  					subtree.reset(get< NP >(*c_root));
  					c_node = subtree.begin();
  					current = std::make_tuple(get< NP >(*c_node), get< DEPTH >(*c_node));
  				}
  			} else {
  				++c_node;
  			  current = std::make_tuple(get< NP >(*c_node), get< DEPTH >(*c_node));
  			}
  			update_simplex();
  			return(*this);
  		};
  		
  		template < bool T = ts >
      typename std::enable_if< T == true, void >::type
    	update_simplex(){ labels = get< LABELS >(*c_node); }
  		
    	template < bool T = ts >
      typename std::enable_if< T == false, void >::type
    	update_simplex(){ return; };
    	
  	}; // iterator
  	
  	// Only start at a non-root node, else return the end
  	auto begin() -> decltype(iterator(*this, init))  { return iterator(*this, init); };
  	auto end() -> decltype(iterator(*this, nullptr)) { return iterator(*this, nullptr); };
  }; // cofaces
  
  // ---- expansion iterator ------
  // template < bool ts = false > 
  // struct expansion : TraversalInterface< ts, expansion > {
  // 	using B = TraversalInterface< ts, expansion >;
  // 	using B::init, B::st, B::p1, B::p2, B::NP, B::DEPTH, B::LABELS;
  // 	using d_node = typename TraversalInterface< ts, expansion >::d_node;
  // 
  // 	expansion(const SimplexTree* st, node_ptr start) : TraversalInterface< ts, expansion >(st, start){ }
  // 
  // 	struct iterator : public TraversalInterface< ts, expansion >::iterator {
  // 		using Bit = typename B::iterator;
  // 		using Bit::current, Bit::labels, Bit::info, Bit::update_simplex, Bit::base, Bit::trie, Bit::current_t_node, Bit::sentinel;
  // 
  // 		// expansion-specific structures 
  // 		const size_t k;
  // 		std::reference_wrapper< node_set_t > c_set; 
  // 	  node_set_t::iterator c_sib;
  // 	  node_set_t::iterator c_node;
  //     vector< node_ptr > intersection; 
  // 		
  // 		// Iterator constructor
  // 		iterator(expansion& dd, node_set_t& exp_set, size_t dim) : TraversalInterface< ts, expansion >::iterator(dd), k(dim) {
  // 			// current = std::make_tuple(cn, dd.st->depth(cn));
  // 		  // update_simplex();
  // 		  // node_set_t& c_set, const idx_t k
  // 		  c_set = std::ref(exp_set);
  // 		  c_node = begin(c_set);
  // 		  c_sib = std::next(c_set.begin(), 1);
  // 		   
  // 		}
  // 		
  // 		node_ptr next_top_node() {
  // 		  if (c_node == end(c_node)){ return(nullptr); }
  // 		  node_ptr top_v = find_vertex((*c_node)->label);
  // 		  while (top_v != nullptr && c_node != end(c_node) && top_v->children.empty()){
  // 		    std::advance(c_node,1);
  // 		    if (c_node == end(c_node)){
  // 		      top_v <- nullptr; 
  // 		    } else {
  // 		      top_v = find_vertex((*c_node)->label);
  // 		    }
  // 		  }
  // 		  return(top_v);
  // 		}
  // 		
  // 		// Increment operator is all that is needed by default
  // 		auto& operator++(){
  // 			// Postcondition: iterator is incremented by one, returns the result
  // 	    node_ptr top_v = next_top_node();
  // 		  node_ptr cn = (*c_node).get(); 
  // 		  
  // 		  if (cn != nullptr && top_v != nullptr){
  // 		    vector< node_ptr > sib_ptrs; 
  //   			std::transform(c_sib, end(c_set), std::back_inserter(sib_ptrs), [](const auto& n){
  //   				return n.get();
  //   			});
  // 
  //         // Get the intersection
  //         intersection.clear();
  //         std::set_intersection(
  //           begin(sib_ptrs), end(sib_ptrs),
  //           begin(top_v->children), end(top_v->children), 
  //           std::back_inserter(intersection), 
  //           [](auto& sib_n, auto& child_n) -> bool {
  //             return node_label(sib_n) < node_label(child_n); 
  //           }
  //         );
  //         
  //         // Insert and recursively expand 
  //         if (intersection.size() > 0){
  //           std::array< idx_t, 1 > int_label;
  //           for (auto& int_node: intersection){
  //             int_label[0] = int_node->label;
  //             auto np = find_it(begin(int_label), end(int_label), cn);
  //             if (np == nullptr){ // not found, need to insert it 
  //               current = { cn, depth(cn) };
  //               
  //             }
  //           }
  //           
  //           expand(cn->children, k-1); // recurse
  //         }
  //         if (siblings != end(c_set)){ ++siblings; }
  //         
  // 		  }
  // 		  
  // 			update_simplex();
  // 			return(*this);
  // 		};
  // 		
  // 		// Doesn't work with regular update method, so override
  //   	constexpr void update_simplex(){
  //   		if constexpr (ts){
  //   		  labels = get< LABELS >(*c_node);
  //   		}
  //   	};
  // 	}; // iterator
  // 	
  // 	// Only start at a non-root node, else return the end
  // 	auto begin(){ return iterator(*this, init); };
  // 	auto end(){ return iterator(*this, nullptr); };
  // }; // expansion
  
  
  // K-skeleton iterator
  template < bool ts = false > 
  struct k_skeleton : preorder< ts > {
    using P = preorder< ts >;
    using B = TraversalInterface< ts, preorder >;
  	using B::init;
  	using t_node = typename B::t_node;
  	using iterator_t = typename P::iterator;
  	
  	k_skeleton(const SimplexTree* st, node_ptr start, const size_t k) 
  	  : preorder< ts >(
  	      st, start, 
          [k](t_node& cn) -> bool { return get< 1 >(cn) <= (k+1); },
          [k](t_node& cn) -> bool { return get< 1 >(cn) <= k; }
        )
  	{};
  	// Only start at a non-root node, else return the end
  	auto begin() -> iterator_t {
  	  return static_cast< P& >(*this).begin();
  	  //return static_cast< P& >(*this).begin();
      // return iterator_t(static_cast< P& >(*this), (node_ptr) init);
  	};
  	auto end() -> iterator_t {
  	  return static_cast< P& >(*this).end();
  	  // return iterator_t(static_cast< P& >(*this), nullptr);
  	}
  };
  
  template < bool ts = false > 
  struct k_simplices : preorder< ts > {
    using P = preorder< ts >;
    using B = TraversalInterface< ts, preorder >;
  	using B::init;
  	using t_node = typename B::t_node;
  	using iterator_t = typename P::iterator;
  	
  	//using d_node = tuple< node_ptr, idx_t >;
  	// const static auto valid_eval = [](const size_t k) { 
  	// 	return([k](t_node& cn) -> bool { return get< 1 >(cn) == (k+1); });
  	// };
  	// const static auto valid_children = [](const size_t k) { 
  	// 	return([k](t_node& cn) -> bool{ return get< 1 >(cn) <= k; });
  	// };
  	k_simplices(const SimplexTree* st, node_ptr start, const size_t k) 
  	  : preorder< ts >(
  	      st, start, 
          [k](t_node& cn) -> bool { return get< 1 >(cn) == (k+1); }, 
          [k](t_node& cn) -> bool { return get< 1 >(cn) <= k; }
        )
  	{};
  	// Only start at a non-root node, else return the end
  	auto begin() -> iterator_t {
  	  return static_cast< P& >(*this).begin();
  	  // return iterator_t(static_cast< P& >(*this), (node_ptr) init);
  	};
  	auto end() -> iterator_t {
  	  return static_cast< P& >(*this).end();
  	  //return iterator_t(static_cast< P& >(*this), nullptr);
  	}
  };
  
  
  template < bool ts = false > 
  struct maximal : preorder< ts > {
    using P = preorder< ts >;
    using B = TraversalInterface< ts, preorder >;
  	using B::init;
  	using t_node = typename B::t_node;
  	using iterator_t = typename P::iterator;
  	// Check if a given node has no children
  	// const static bool has_no_children(const t_node& cn){ return(get< 0 >(cn)->children.empty()); }
  	// 
  	// // Check cn is only coface 
  	// const static bool is_only_root(const SimplexTree* st, const t_node& cn){
  	// 	auto cr = coface_roots< ts >(st, get< 0 >(cn));
  	// 	return(std::next(cr.begin()) == cr.end());
  	// }
  	// 
  	// A valid member in a preorder-traversal is just 
  	// const static auto valid_eval = [](const SimplexTree* st) { 
  	// 	return([&st](t_node& cn) -> bool { 
  	// 		return get<0>(cn) == nullptr ? false : has_no_children(cn) && is_only_root(st, cn);
  	// 	});
  	// };
  	// const static auto always_true = [](t_node& cn) -> bool { return true; };
  
  	// Constructor is all that is needed
  	maximal(const SimplexTree* st, node_ptr start) : 
  		preorder< ts >(st, start, [st](t_node& cn) -> bool { 
    		node_ptr np = get< 0 >(cn);
    		if (np != nullptr && np != st->root.get()){
    		  auto cr = coface_roots< false >(st, np);
    		  return(np->children.empty() && std::next(cr.begin()) == cr.end());
    		}
    		return false;
  		}, [](t_node& cn) -> bool { return true; })
  	{};
  	auto begin() -> iterator_t {
  	  return static_cast< P& >(*this).begin();
  	  // return iterator_t(static_cast< P& >(*this), (node_ptr) init);
  	};
  	auto end() -> iterator_t {
  	  return static_cast< P& >(*this).end();
  	  // return iterator_t(static_cast< P& >(*this), nullptr);
  	}
  
  };
  
  // Checks for empty intersection
  inline bool empty_intersection(const vector< idx_t > x, const vector< idx_t > y){
    vector<idx_t>::const_iterator i = x.begin(), j = y.begin();
    while (i != x.end() && j != y.end()){
      if (*i<*j) ++i; else if(*j<*i) ++j; else return false;
    }
    return true;
  }
  
  template< typename tnode_t >
  std::function< bool(tnode_t&) > link_condition(const SimplexTree* st, const node_ptr s_np){
    const simplex_t s = st->full_simplex(s_np);
		return([st, s](tnode_t& cn) -> bool {
			bool is_link = false;
			const simplex_t t = st->full_simplex(get< 0 >(cn));
			bool is_disjoint = empty_intersection(t, s);
			if (is_disjoint){
				vector< idx_t > pot_link;
      	std::set_union(s.begin(), s.end(), t.begin(), t.end(), std::back_inserter(pot_link));
      	if (st->find_it(begin(pot_link), end(pot_link), st->root.get()) != nullptr){
      	  is_link = true;
      	}
			}
			return(is_link);
		});
  }
  
  template < bool ts = false > 
  struct link : preorder< ts > {
  	using P = preorder< ts >;
    using B = TraversalInterface< ts, preorder >;
  	using B::init;
  	using t_node = typename B::t_node;
  	using iterator_t = typename P::iterator;
  
  //   template < bool T = ts >
  // 	static auto get_simplex2(const SimplexTree* st, t_node& cn) -> std::enable_if< T == true, simplex_t > {
  // 		return(get< 2 >(cn));
  // 	}
  //   template < bool T = ts >
  // 	static auto get_simplex2(const SimplexTree* st, t_node& cn) -> std::enable_if< T == false, simplex_t > {
  // 		return(st->full_simplex(get< 0 >(cn)));
  // 	}
  
  	// A valid member in a preorder-traversal is just 
  // 	const auto valid_eval = [](const SimplexTree* st, const node_ptr s_np) {
  // 		const simplex_t s = st->full_simplex(s_np);
  // 		return([st, s](t_node& cn) -> bool {
  // 			bool is_link = false;
  // 			const simplex_t t = st->full_simplex(get< 0 >(cn));
  // 			bool is_disjoint = empty_intersection(t, s);
  // 			if (is_disjoint){
  // 				vector< idx_t > pot_link;
  //       	std::set_union(s.begin(), s.end(), t.begin(), t.end(), std::back_inserter(pot_link));
  //       	if (st->find_it(begin(pot_link), end(pot_link), st->root.get()) != nullptr){
  //       	  is_link = true;
  //       	}
  // 			}
  // 			return(is_link);
  // 		});
  // 	};
  // 	static const auto always_true = [](t_node& cn) -> bool { return true; };
  
  
  // [st, &start](t_node& cn) -> bool {
  // 		  bool is_link = false;
  // 		  std::cout << "testing: \n";        
  // 			std::cout << "1: " << start << std::endl;
  // 		  st->print_simplex(std::cout, get<0>(cn), true);
  //       st->print_simplex(std::cout, start, true);
  // 			// std::cout << "2: " << start << std::endl;
  // 		  const simplex_t s = st->full_simplex(start);
  // 			const simplex_t t = st->full_simplex(get< 0 >(cn));
  // 			// std::cout << "3: " << start << std::endl;
  // 			bool is_disjoint = empty_intersection(t, s);
  // 			if (is_disjoint){
  // 				vector< idx_t > pot_link;
  //       	std::set_union(s.begin(), s.end(), t.begin(), t.end(), std::back_inserter(pot_link));
  //       	if (st->find_it(pot_link.begin(), pot_link.end(), st->root.get()) != nullptr){
  //       	  is_link = true;
  //       	}
  // 			}
  // 			// std::cout << "4: " << start << std::endl;
  // 			std::cout << "is_link: " << is_link << std::endl;        
  // 			return(is_link);
  // 		}
  
  	// Constructor is all that is needed
  	link(const SimplexTree* st, const node_ptr start) : 
  		preorder< ts >(
  		    st, st->root.get(), 
  		    link_condition< t_node >(st, start),
  		    [](const t_node& cn) -> bool { return true; }
  		)
  	{};
  	auto begin() -> iterator_t {
  	  return static_cast< P& >(*this).begin();// maybe change this to not be P&
  	  // return iterator_t(static_cast< P& >(*this), (node_ptr) init);
  	};
  	auto end() -> iterator_t {
  	  return static_cast< P& >(*this).end();
  	  // return iterator_t(static_cast< P& >(*this), nullptr);
  	}
  }; // link
  
  template< typename tnode_t >
  std::function<bool(tnode_t&)> face_condition(const SimplexTree* st, node_ptr start){
    const simplex_t sigma = st->full_simplex(start);
    return [st, sigma](tnode_t& cn) -> bool { 
      node_ptr np = get<0>(cn);
      if (np == nullptr || np == st->root.get()){ return false; }
      simplex_t tau = st->full_simplex(np, get<1>(cn));

      return std::includes(sigma.begin(), sigma.end(), tau.begin(), tau.end());
		};
  };
  
  template< typename tnode_t >
  std::function<bool(tnode_t&)> face_condition2(const SimplexTree* st, node_ptr start){
    const size_t d = st->depth(start);
    return [d](tnode_t& cn) -> bool { return get< 1 >(cn) <= d; };
  };
  
  template < bool ts = false > 
  struct faces : level_order< ts > {
  	using P = level_order< ts >;
    using B = TraversalInterface< ts, level_order >;
  	using B::init;
  	using t_node = typename B::t_node;
  	using iterator_t = typename P::iterator;
  
  	// A valid member in a preorder-traversal is just 
  	// const static auto valid_eval = [](const SimplexTree* st, node_ptr start) { 
  	// 	simplex_t sigma = st->full_simplex(start);
  	// 	return([sigma](t_node& cn) -> bool { 
  	// 	  return SimplexTree::is_face(get< 2 >(cn), sigma);
  	// 	});
  	// };
  	// const static auto valid_children = [](const SimplexTree* st, node_ptr start) { 
  	//   idx_t k = st->depth(start);
  	// 	return([k](t_node& cn) -> bool{ return get< 1 >(cn) <= k; });
  	// };
  	// Constructor is all that is needed
  	faces(const SimplexTree* st, node_ptr start) : 
  		level_order< ts >(
  		  st, st->root.get(),
  		  // [](t_node& cn) -> bool { return true; }, 
  		  face_condition< t_node >(st, start),
  		  face_condition2< t_node >(st, start)
  		  // [&st, &start](t_node& cn) ->bool { return get< 1 >(cn) <= st->depth(start); }
        // face_condition< t_node >(st, start),
    		// [&st, &start](t_node& cn) ->bool { return get< 1 >(cn) <= st->depth(start); }
  	  )
  	{};
  	auto begin() -> iterator_t {
  	  return static_cast< P& >(*this).begin();
  	  // return iterator_t(static_cast< P& >(*this), (node_ptr) init);
  	};
  	auto end() -> iterator_t {
  	  return static_cast< P& >(*this).end();
  	  // return iterator_t(static_cast< P& >(*this), nullptr);
  	}
  };
  
  
  template <class T>
  auto get_node_ptr(T& cn) -> node_ptr { return std::get< 0 >(cn); }
  
  template <class T>
  auto get_depth(T& cn) -> idx_t { return std::get< 1 >(cn); }
  
  template <class T>
  auto get_simplex(T& cn) -> simplex_t { return std::get< 2 >(cn); }
  
  // Generic traversal function which unpacks the tuple and allows for early termination of the iterable
  template <class Iterable, class Lambda> 
  auto traverse(Iterable traversal, Lambda f) -> typename std::enable_if< Iterable::is_tracking, void >::type {
  	for (auto& cn: traversal){ 
  		bool should_continue = f(get_node_ptr(cn), get_depth(cn), get_simplex(cn));
  		if (!should_continue){ break; }
  	}
  }
  template <class Iterable, class Lambda> 
  auto traverse(Iterable traversal, Lambda f) -> typename std::enable_if< !Iterable::is_tracking, void >::type {
  	for (auto& cn: traversal){ 
  		bool should_continue = f(get_node_ptr(cn), get_depth(cn));
  		if (!should_continue){ break; }
  	}
  }
  

  // template < class Iterable, typename Lambda > 
  // void traverse_node_pairs(Iterable traversal, Lambda&& f) -> typename std::enable_if< Iterable::is_tracking, void >::type {
  //   for (auto& cn: traversal){ 
  //     f(std::make_pair(get_node_ptr(cn), get_depth(cn)));
  //   }
  // }
  // template < class Iterable, typename Lambda > 
  // void traverse_node_pairs(Iterable traversal, Lambda&& f) -> typename std::enable_if< !Iterable::is_tracking, void >::type {
  //   for (auto& cn: traversal){ 
  //     f(std::make_pair(get_node_ptr(cn), get_depth(cn)));
  //   }
  // }
  // 
  // template < class Iterable, typename Lambda > 
  // void traverse_simplices(Iterable traversal, Lambda&& f) -> typename std::enable_if< Iterable::is_tracking, void >::type {
  //   for (auto& cn: traversal){ 
  //     f(get_simplex(cn));
  //   }
  // }
  // template < class Iterable, typename Lambda > 
  // void traverse_simplices(Iterable traversal, Lambda&& f) -> typename std::enable_if< !Iterable::is_tracking, void >::type {
  //   for (auto& cn: traversal){ 
  //     f(get_simplex(cn));
  //   }
  // }
  
   // Generic traversal function which unpacks the tuple and allows for early termination of the iterable
  // template <class Iterable, class Lambda> 
  // decltype(auto) generate(Iterable traversal, Lambda f){
  //   using V = typename Iterable::value_type;
  //   using T = decltype( std::declval< Lambda >()(V) );
  //   vector< T > result; 
  // 	for (auto& cn: traversal){ 
  // 		bool should_continue = std::apply(f, cn);
  // 		if (!should_continue){ break; }
  // 	}
  // }
  
  
  
}; // end namespace st 


#endif 


