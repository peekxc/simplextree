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

using std::get; 
using simplex_t = SimplexTree::simplex_t;
using node_ptr = SimplexTree::node_ptr;
using node_uptr = SimplexTree::node_uptr;

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
  using namespace detail; 
    

  template <  bool ts, template< bool > class Derived > 
  struct TraversalInterface {
  	using d_type = Derived< ts >;
  	using d_node = std::tuple< node_ptr, idx_t >;
  	using t_node = typename std::conditional< ts, std::tuple< node_ptr, idx_t, simplex_t >, d_node >::type;
  	using pair_pred_t = delegate< bool (t_node&) >;
  	// using pair_pred_t = std::function< bool (t_node&) >;
  
  	public: 
  		// enum : idx_t { np = 0, depth = 1, simplex = 2 };
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
  				if constexpr (ts){	
  					labels = simplex_t();
  					labels.reserve(dd.st->tree_max_depth);
  				}
  			};
  
  			constexpr d_type& base() const {
  				return(this->info.get());
  			};
  			constexpr t_node& current_t_node(){
  				if constexpr (ts){ 
  					output = std::tuple_cat(current, std::make_tuple(labels)); 
  					return(output);
  				} 
  				else { return(current); }
  			}
  
  			// Operators
  			auto operator==(const d_iter& t) const {
  				if constexpr (has_equality< d_iter >::value){
  					return static_cast< d_iter* >(this)->operator==(t);
  				} else {
  					return (get< 0 >(t.current) == get< 0 >(current)); // && (get< 1 >(t.current) == get< 1 >(current)
  				}
  			}
  			auto operator!=(const d_iter& t) const {
  				if constexpr (has_not_equality< d_iter >::value){
  					return static_cast< d_iter* >(this)->operator!=(t);
  				} else {
  					return !(*this == t);
  				}
  			}
  			auto& operator*() {
  				if constexpr (has_dereference< d_iter >::value){
  					return static_cast< d_iter* >(this)->operator*();
  				} else {
  					return current_t_node();
  				}
  			};
  
  			// Updates the current simplex labels, if tracking
  			constexpr void update_simplex() noexcept {
  				if constexpr(has_update_simplex< d_type >::value){
  					return static_cast< d_type* >(this)->update_simplex();
  				} else {
  					if constexpr (ts){
  						if (get< NP >(current) != nullptr && get< DEPTH >(current) > 0){
  							labels.resize(get< DEPTH >(current));
  							labels[get< DEPTH >(current)-1] = get< NP >(current)->label;
  						}
  					}
  				}
  			}; // update_simplex
  		};
  		
  		// methods within TraversalInterface can use template to access members of Derived
  		auto begin() {
  			if constexpr(has_begin< d_type >::value){
  				return static_cast< d_type* >(this)->begin();
  			} else {
  				return d_type::iterator(static_cast< d_type& >(*this), init);
  			}
  		};
  		auto end() {
  			if constexpr(has_end< d_type >::value){
  				return static_cast< d_type* >(this)->end();
  			} else {
  				using d_iter = typename d_type::iterator;
  				return d_iter(static_cast< d_type& >(*this), nullptr);
  			}
  		};
  };
  // https://eli.thegreenplace.net/2011/05/17/the-curiously-recurring-template-pattern-in-c
  
  
  // Preorder traversal iterator
  template < bool ts = false > 
  struct preorder : TraversalInterface< ts, preorder > {
  	using B = TraversalInterface< ts, preorder >;
  	using B::init, B::st, B::p1, B::p2, B::NP, B::DEPTH, B::LABELS;
  
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
  		using Bit::current, Bit::labels, Bit::info, Bit::update_simplex, Bit::base, Bit::current_t_node, Bit::sentinel;
  
  		// DFS specific data structures
  		std::stack< typename B::d_node > node_stack;  
  
  		// Iterator constructor
  		iterator(preorder& dd, node_ptr cn = nullptr) : TraversalInterface< ts, preorder >::iterator(dd){
  			current = std::make_tuple(cn, dd.st->depth(cn));
  		}
  		
  		// Increment operator is all that is needed by default
  		auto& operator++(){
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
  	};
  
  	// Only start at a non-root node, else return the end
  	auto begin(){
  		if (init == st->root.get()){
  			return st->n_simplexes.empty() ? iterator(*this, nullptr) : ++iterator(*this, st->root.get()); 
  		} else {
  			return iterator(*this, init);
  		}
  	};
  };
  
  // level order traversal iterator
  template < bool ts = false > 
  struct level_order : TraversalInterface< ts, level_order > {
  	using B = TraversalInterface< ts, level_order >;
  	using B::init, B::st, B::p1, B::p2, B::NP, B::DEPTH, B::LABELS;
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
  		using Bit::current, Bit::labels, Bit::info, Bit::update_simplex, Bit::base, Bit::current_t_node, Bit::sentinel;
  
  		// BFS specific data structures
  		std::queue< typename B::d_node > node_queue;  
  
  		// Iterator constructor
  		iterator(level_order& dd, node_ptr cn = nullptr) : TraversalInterface< ts, level_order >::iterator(dd){
  			current = std::make_tuple(cn, dd.st->depth(cn));
  		}
  		
  		// Increment operator is all that is needed by default
  		auto& operator++(){
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
  		
  		// Doesn't work with regular update method, so override
    	constexpr void update_simplex(){
    		if constexpr (ts){
    			labels = base().st->full_simplex(get< NP >(current), get< DEPTH >(current));
    		}
    	};
  	};// iterator 
  
  	// Only start at a non-root node, else return the end
  	auto begin(){
  		if (init == st->root.get()){
  			return st->n_simplexes.empty() ? iterator(*this, nullptr) : ++iterator(*this, st->root.get()); 
  		} else {
  			return iterator(*this, init);
  		}
  	};
  };
  
  
  // Coface-roots search
  template < bool ts = false >
  struct coface_roots : TraversalInterface< ts, coface_roots > {
  	using B = TraversalInterface< ts, coface_roots >;
  	using B::init, B::st, B::p1, B::p2, B::NP, B::DEPTH, B::LABELS;
  	using d_node = typename TraversalInterface< ts, coface_roots >::d_node;
  
  	// Constructors
  	coface_roots() : TraversalInterface< ts, coface_roots >() {};
  	coface_roots(const SimplexTree* st_, node_ptr start = nullptr) : TraversalInterface<  ts, coface_roots >(st_, start) {};
  	template< typename P1, typename P2 >
  	coface_roots(const SimplexTree* st_, node_ptr start, P1 pred1, P2 pred2) : TraversalInterface<  ts, coface_roots >(st_, start, pred1, pred2){};
  
  	// Iterator type
  	struct iterator : public TraversalInterface< ts, coface_roots >::iterator {
  		using Bit = typename B::iterator;
  		using Bit::current, Bit::labels, Bit::info, Bit::update_simplex, Bit::base, Bit::current_t_node, Bit::sentinel;
  
  		// coface-roots-specific structures 
  		simplex_t start_coface_s; 
  		// SimplexTree::level_map::value_type;
  		// vector< node_ptr >* current_cousins { nullptr }; 
  		// vector< node_ptr >::const_iterator current_cousin;
  		size_t c_level_key = 0; // the current level_map key
  		size_t c_level_idx = 0; 
  
  		// Iterator constructor
  		iterator(coface_roots& dd, node_ptr cn) : TraversalInterface< ts, coface_roots >::iterator(dd){
  			// if (cn == nullptr || cn == dd.st->root.get()){ throw std::invalid_argument("Invalid given coface."); };
  			const size_t c_depth = dd.st->depth(cn); 
  			current = std::make_tuple(cn, c_depth);
  			start_coface_s = dd.st->full_simplex(cn, c_depth);
  		}
  		
  		// Increment operator is all that is needed by default
  		auto& operator++(){
  			// If root was given, end the traversal immediately.
  			if (get< NP >(current) == base().st->root.get()){ 
  				current = sentinel();
  				return(*this);
  			}
  			// 1. Increment the depth until a cousin is found; else if exceed tree the depth, set current = { nullptr, 0 };
  			// 2. Assume current depth is <= tree depth. Iterate through cousins, checking whether they are coface roots.
  			while(get< DEPTH >(current) <= base().st->tree_max_depth && base().st->level_map.count(c_level_key) == 0) {
  				get< DEPTH >(current) = get< DEPTH >(current)+1;
  				
  				// Look for cousins at next depth going up 
  				auto tmp_key = base().st->encode_node(base().init->label, get< DEPTH >(current)); 
  				auto ni = base().st->level_map.find(tmp_key); 
  				if (ni == base().st->level_map.end()){
  					c_level_key = 0; 
  				} else {
  					c_level_key = tmp_key;
  				}
  				c_level_idx = 0; 
  			}
  
  			// If we've passed the tree's max depth, end the traversal
  			if (get< DEPTH >(current) > base().st->tree_max_depth || c_level_key == 0){
  				current = sentinel();
  			} else {
  				// Otherwise iterate through the cousins
  				auto cousins = base().st->level_map.at(c_level_key);
  				node_ptr c_cousin = cousins.at(c_level_idx);
  				if (base().st->is_face(start_coface_s, base().st->full_simplex(c_cousin, get< DEPTH >(current)))){ 
  						get< NP >(current) = c_cousin; // should be guarenteed to exist
  				}
  				if (cousins.size() >= (c_level_idx + 1)){
  					c_level_key = 0; 
  					c_level_idx = 0; 
  				} else {
  					c_level_idx++;
  				}
  			}
  			update_simplex();
  			return *this; 
  		}; // operator++
  
  		// Doesn't work with regular update method, so override
  		constexpr void update_simplex(){
  			if constexpr (ts){
  				labels = base().st->full_simplex(get< NP >(current));
  			}
  		};
  	}; // iterator 
  
  	// Only start at a non-root node, else return the end
  	auto begin(){ return iterator(*this, init); };
  	auto end(){ return iterator(*this, nullptr); };
  
  }; // end coface_roots 
  
  // ---- coface iterator ------
  template < bool ts = false > 
  struct cofaces : TraversalInterface< ts, cofaces > {
  	using B = TraversalInterface< ts, cofaces >;
  	using B::init, B::st, B::p1, B::p2, B::NP, B::DEPTH, B::LABELS;
  	using d_node = typename TraversalInterface< ts, cofaces >::d_node;
  
  	cofaces(const SimplexTree* st, node_ptr start) : TraversalInterface< ts, cofaces >(st, start){ }
  
  	struct iterator : public TraversalInterface< ts, cofaces >::iterator {
  		using Bit = typename B::iterator;
  		using Bit::current, Bit::labels, Bit::info, Bit::update_simplex, Bit::base, Bit::current_t_node, Bit::sentinel;
  
  		using preorder_it = decltype(std::declval< preorder< ts > >().begin());
  		using coface_root_it = decltype(std::declval< coface_roots< ts > >().begin());
  
  		// coface-roots-specific structures 
  		coface_roots< ts > roots; 
  		coface_root_it c_root;
  		preorder< ts > subtree; 
  		preorder_it c_node; 
  
  		// Iterator constructor
  		iterator(cofaces& dd, node_ptr cn) : TraversalInterface< ts, cofaces >::iterator(dd),
  		 roots(coface_roots< ts >(dd.st, cn)), c_root(roots, cn), subtree(preorder< ts >(dd.st, cn)), c_node(subtree.begin()) {
  			 current = std::make_tuple(cn, dd.st->depth(cn));
  		   update_simplex();
  		}
  		
  		// Increment operator is all that is needed by default
  		auto& operator++(){
  			// Need to increment iterator by one and return the result 
  			// Logically what needs to happen is: 
  			// - if next node is end of subtree 
  			//   - ... then if next root is end of coface_roots, we're finished return end 
  			//   - else next root is not end of coface_roots, advance root and reset subtree
  			// else while next node is not at end of subtree 
  			// - advance one in subtree, report it back
  			if (get< NP >(*c_root) == base().st->root.get()){ ++c_root; }
  			if (std::next(c_node) == subtree.end()){
  				if (c_root == roots.end()){
  					current = sentinel();
  				} else {
  					++c_root;
  					subtree.reset(get< NP >(*c_root));
  					c_node = subtree.begin();
  					current = { get< NP >(*c_node), get< DEPTH >(*c_node) };
  				}
  			} else {
  				++c_node;
  				current = { get< NP >(*c_node), get< DEPTH >(*c_node) };
  			}
  			update_simplex();
  			return(*this);
  		};
  		
  		// Doesn't work with regular update method, so override
    	constexpr void update_simplex(){
    		if constexpr (ts){
    			labels = base().st->full_simplex(get< NP >(current), get< DEPTH >(current));
    		}
    	};
  	}; // iterator
  	
  	// Only start at a non-root node, else return the end
  	auto begin(){ return iterator(*this, init); };
  	auto end(){ return iterator(*this, nullptr); };
  }; // cofaces
  
  // K-skeleton iterator
  template < bool ts = false > 
  struct k_skeleton : preorder< ts > {
  	using t_node = typename TraversalInterface< ts, preorder >::t_node;
  	static constexpr auto valid_eval = [](const size_t k) constexpr { 
  		return([k](t_node& cn) -> bool { return get< 1 >(cn) <= (k+1); });
  	};
  	static constexpr auto valid_children = [](const size_t k) constexpr { 
  		return([k](t_node& cn) -> bool{ return get< 1 >(cn) <= k; });
  	};
  	k_skeleton(const SimplexTree* st, node_ptr start, const size_t k) : preorder< ts >(st, start, valid_eval(k), valid_children(k)){}
  };
  
  template < bool ts = false > 
  struct max_skeleton : preorder< ts > {
  	using t_node = typename TraversalInterface< ts, preorder >::t_node;
  	//using d_node = tuple< node_ptr, idx_t >;
  	static constexpr auto valid_eval = [](const size_t k) constexpr { 
  		return([k](t_node& cn) -> bool { return get< 1 >(cn) == (k+1); });
  	};
  	static constexpr auto valid_children = [](const size_t k) constexpr { 
  		return([k](t_node& cn) -> bool{ return get< 1 >(cn) <= k; });
  	};
  	max_skeleton(const SimplexTree* st, node_ptr start, const size_t k) : preorder< ts >(st, start, valid_eval(k), valid_children(k)){}
  };
  
  
  template < bool ts = false > 
  struct maximal : preorder< ts > {
  	using t_node = typename TraversalInterface< ts, preorder >::t_node; 
  
  	// Check if a given node has no children
  	static constexpr bool has_no_children(const t_node& cn){ return(get< 0 >(cn)->children.empty()); }
  	
  	// Check cn is only coface 
  	static constexpr bool is_only_root(const SimplexTree* st, const t_node& cn){
  		auto cr = coface_roots< ts >(st, get< 0 >(cn));
  		return(std::next(cr.begin()) == cr.end());
  	}
  
  	// A valid member in a preorder-traversal is just 
  	static constexpr auto valid_eval = [](const SimplexTree* st) constexpr { 
  		return([st](t_node& cn) -> bool { 
  			return get<0>(cn) == nullptr ? false : has_no_children(cn) && is_only_root(st, cn);
  		});
  	};
  	static constexpr auto always_true = [](t_node& cn) -> bool { return true; };
  
  	// Constructor is all that is needed
  	maximal(const SimplexTree* st, node_ptr start) : 
  		preorder< ts >(st, start, valid_eval(st), always_true){}
  
  };
  
  // Checks for empty intersection
  inline bool empty_intersection(const vector< idx_t > x, const vector< idx_t > y){
    vector<idx_t>::const_iterator i = x.begin(), j = y.begin();
    while (i != x.end() && j != y.end()){
      if (*i<*j) ++i; else if(*j<*i) ++j; else return false;
    }
    return true;
  }
  
  template < bool ts = false > 
  struct link : preorder< ts > {
  	using t_node = typename TraversalInterface< ts, preorder >::t_node; 
  
  	static inline simplex_t get_simplex(const SimplexTree* st, t_node& cn){
  		if constexpr (ts){ return(get< 2 >(cn)); }
  		else{ return(st->full_simplex(get< 0 >(cn))); }
  	}
  
  	// A valid member in a preorder-traversal is just 
  	static inline auto valid_eval = [](const SimplexTree* st, node_ptr s_np) { 
  		const simplex_t s = st->full_simplex(s_np);
  		return([st, s](t_node& cn) -> bool { 
  			bool is_link = false; 
  			const simplex_t t = get_simplex(st, cn); 
  			bool is_disjoint = empty_intersection(t, s);
  			if (is_disjoint){
  				// potential_link.resize(0);
  				vector< idx_t > potential_link;
        	std::set_union(s.begin(), s.end(), t.begin(), t.end(), std::back_inserter(potential_link)); 
        	if (st->find_node(potential_link) != nullptr){ is_link = true;  }
  			}
  			return(is_link);
  		});
  	};
  	static constexpr auto always_true = [](t_node& cn) -> bool { return true; };
  
  	// Constructor is all that is needed
  	link(const SimplexTree* st, node_ptr start) : 
  		preorder< ts >(st, st->root.get(), valid_eval(st, start), always_true){}
  }; // link
  
  
  template < bool ts = false > 
  struct faces : level_order< ts > {
  	using t_node = typename TraversalInterface< ts, level_order >::t_node; 
  
  	// A valid member in a preorder-traversal is just 
  	static constexpr auto valid_eval = [](const SimplexTree* st, node_ptr start) { 
  		simplex_t sigma = st->full_simplex(start);
  		return([sigma](t_node& cn) -> bool { 
  		  return SimplexTree::is_face(get< 2 >(cn), sigma);
  		});
  	};
  	static constexpr auto valid_children = [](const SimplexTree* st, node_ptr start) constexpr { 
  	  idx_t k = st->depth(start);
  		return([k](t_node& cn) -> bool{ return get< 1 >(cn) <= k; });
  	};
  	// Constructor is all that is needed
  	faces(const SimplexTree* st, node_ptr start) : 
  		level_order< ts >(st, start, valid_eval(st, start), valid_children(st, start)){}
  };
  
  // Generic traversal function which unpacks the tuple and allows for early termination of the iterable
  template <class Iterable, class Lambda> 
  void traverse(Iterable traversal, Lambda f){
  	for (auto& cn: traversal){ 
  		bool should_continue = std::apply(f, cn);
  		if (!should_continue){ break; }
  	}
  }
  
  template <class T>
  constexpr auto get_node_ptr(T& cn){ return std::get< 0 >(cn); }
  
  template <class T>
  constexpr auto get_depth(T& cn){ return std::get< 1 >(cn); }
  
  template <class T>
  constexpr auto get_simplex(T& cn){ return std::get< 2 >(cn); }
  
  // template <class Iterable, class Lambda> 
  // auto traverse(Iterable traversal, Lambda f){
  // 	auto v = std::generate(traversal.begin(), traversal.end(), f);
  // 	return(v);
  // }
  // 
  
  
}; // end namespace st 


#endif 