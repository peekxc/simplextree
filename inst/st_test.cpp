
#include <simplextree.h>
#include <iostream> 
#include <array>

using node_ptr = SimplexTree::node_ptr;
using std::array;
using std::get;  

int main(){
	// Make simplex tree 
	const size_t n = 4; 
	auto st = SimplexTree();

	// Insert vertices + edges
	simplex_t sigma = simplex_t(1,0);
	for (size_t i = 0; i < n; ++i){
		sigma[0] = i;
		st.insert_simplex(sigma);
	}
	for (int i = 0; i < n; ++i){
		for (int j = i+1; j < n; ++j){
			simplex_t tau = {std::size_t(i), std::size_t(j)};
			st.insert_simplex(tau);
		}
	}
	st.expansion(4);
	st.print_tree();
	// st.print_cousins();

	// DFS 
	std::cout << "Preorder traversal: " << std::flush << std::endl;
	auto st_dfs = preorder< true >(&st);
	for (auto& cn: st_dfs){
		std::cout << "{ ";
		for (auto vi: get< 2 >(cn)){ std::cout << vi << ", "; };
		std::cout << " }" << std::endl;
		//st.print_simplex(get<0>(cn));
	}
	traverse(preorder< true >(&st), [](node_ptr cn, idx_t depth, simplex_t sigma){ 
		// node_ptr node = get<0>(cn);
		std::cout << depth << std::endl;
		return true; 
	});
	// traverse(&st, PRE_ORDER, [](auto& cn){
	// 	node_ptr node = get<0>(cn);
	// });

	std::cout << "Preorder traversal 2: " << std::flush << std::endl;
	auto st_dfs2 = st_dfs;
	std::cout << "Root before copy: " << st_dfs.st->root.get() << std::endl;
	std::cout << "Root after copy: " << st_dfs2.st->root.get() << std::endl;
	for (auto& cn: st_dfs2){
		st.print_simplex(get<0>(cn));
	}

	// 1-skeleton
	std::cout << "One skeleton: " << std::flush << std::endl;
	auto one_skeleton = k_skeleton< true >(&st, st.root.get(), 1);
	for (auto& cn: one_skeleton){
		st.print_simplex(get<0>(cn));
	}
	std::cout << std::flush;

	// coface roots
	st.clear();
	st.insert_simplex(simplex_t{ 1, 2, 3 });
	st.insert_simplex(simplex_t{ 2, 3, 4, 5 });
	st.print_tree();
	auto face = st.find_node(simplex_t{ 3, 4 });
	auto cr = coface_roots< false >(&st, face);
	std::cout << "Coface roots: " << std::flush << std::endl;
	for (auto& cn: cr){
		st.print_simplex(get<0>(cn));
	}

	std::cout << "Cofaces 1: " << std::flush << std::endl;
	auto cr1 = coface_roots< true >(&st, face);
	for (auto& cn: cr1){
		st.print_simplex(get<0>(cn));
		auto co1 = preorder< true >(&st, get< 0 >(cn));
		for (auto& cn2: co1){
			st.print_simplex(get<0>(cn2));
		}
	}

	std::cout << "Cofaces 2: " << std::flush << std::endl;
	auto co = cofaces< true >(&st, face);
	for (auto& cn: co){
		st.print_simplex(get<0>(cn));
	}

	std::cout << "Preorder traversal: " << std::flush << std::endl;
	st_dfs = preorder< true >(&st, st.root.get());
	for (auto& cn: st_dfs){
		st.print_simplex(get<0>(cn));
	}

	// auto dfs3 = preorder< false >(&st);
	// auto it = dfs3.begin();
	// for (; it != dfs3.end(); ++it){
	// 	auto it2 = it; 
	// 	st.print_simplex(get<0>(*it2));
	// }


	std::cout << "Maximal simplices: " << std::flush << std::endl;
	auto ma = maximal< false >(&st, st.root.get());
	for (auto it=ma.begin(); it != ma.end(); ++it){
		st.print_simplex(get<0>(*it));
		// std::cout << "here" << std::flush << std::endl;
	}
	// std::cout << int(ma.begin() == ma.end()) << std::endl;
	// for (auto& cn: ma){
	// 	st.print_simplex(get<0>(cn));
	// }


	std::cout << "BFS: " << std::flush << std::endl;
	auto bfs = level_order< false >(&st);
	for (auto& cn: bfs){
		st.print_simplex(get<0>(cn));
	}

	std::cout << "Link: " << std::flush << std::endl;
	auto one = st.find_node({ 1 });
	auto link_ = link(&st, one);
	for (auto& cn: link_){
		st.print_simplex(get<0>(cn));
	}

	// const auto postorder = [](){

	// };
	// visit_postorder(node->left);
  //           visit_postorder(node->right);
  //           dispatch_node(node);

	// DFS 
	// std::cout << "DFS: " << std::endl;
	// auto dfs_traversal = dfs< false >(&st, st.root.get());
	// for (auto& node: dfs_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }

	// BFS 
	// std::cout << "BFS: " << std::endl;
	// auto bfs_traversal = bfs< false >(&st);
	// for (auto& node: bfs_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }

	// simplex_t tau = { std::size_t(0), std::size_t(1), std::size_t(2) };
	// std::cout << "coface_roots: " << std::endl;
	// auto wut = st.find_node(tau); 
	// printf("%p \n", wut);
	// auto cr_traversal = coface_roots< false >(&st, wut);
	// for (auto& node: cr_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }

	// std::cout << "coface_roots of root: " << std::endl;
	// cr_traversal = coface_roots< false >(&st, st.root.get());
	// for (auto& node: cr_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }

	// simplex_t tau013 = { std::size_t(0), std::size_t(1), std::size_t(3) };
	// std::cout << "coface_roots of 0 1 3: " << std::endl;
	// cr_traversal = coface_roots< false >(&st, st.find_node(tau013));
	// for (auto& node: cr_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }
	
	// TODO: MUST SWITCH TO INNER ITERATOR CLASS TO PREVENT MULTIPEL CONSTRUCTORS BEING CALLED!


	// std::cout << "cofaces: " << std::endl;
	// auto cofaces_traversal = cofaces< false >(&st, st.find_node(tau));
	// for (auto& node: cofaces_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }

	// std::cout << "cofaces of root: " << std::endl;
	// cofaces_traversal = cofaces< false >(&st, st.root.get());
	// for (auto& node: cofaces_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }

	// std::cout << "maximal simplices: " << std::endl;
	// auto ms_traversal = maximal_simplices< false >(&st, st.root.get());
	// for (auto& node: ms_traversal){
	// 	auto [sigma, d] = node;
	// 	st.print_simplex(sigma);
	// }


	// 

	// st.expansion(4);
	// st.print_tree();
	// auto three_node = st.find_node({ std::size_t(3) });
	// std::cout << three_node->children.empty() << std::endl; 

	// simplex_t tau = { std::size_t(0), std::size_t(1), std::size_t(2) };
	// node_ptr cn = st.find_node(tau);
	// auto tau_labels = st.full_simplex_r< 3 >(cn);
	// tau_labels = st.full_simplex(cn, 3);
	
	// simplex_t tau2 = { std::size_t(0), std::size_t(1), std::size_t(2), std::size_t(3)  };
	// auto tau2_labels = st.full_simplex(st.find_node(tau2), 4);

	//auto cf = coface_roots< true >(&st, st.root.get());
	// auto cf = coface_roots< true >(&st, st.find_node(tau));
	// //auto cf = coface_roots< true >(&st, st.find_node(tau2));
	// for (auto& node: cf){
	// 	auto [sigma, d, s_labels] = node; 
	// 	expand_subtrees(sigma);
	// 	auto s_labels2 = st.full_simplex(sigma);
	//  	for (auto& label: s_labels2){ std::cout << label << ", "; }
	// 	std::cout << std::endl;
	// 	//for (auto& label: s_labels){ std::cout << label << ", "; }
	// 	// std::cout << std::endl;
	// }
	

	// auto dfs_traversal = dfs< false >(&st);
	// for (auto& node: dfs_traversal){
	// 	auto [sigma, d] = node; 
	// 	std::cout << sigma->label << ", ";
	// }

	// auto bfs_traversal = bfs< false >(&st);
	// for (auto& node: bfs_traversal){
	// 	auto [sigma, d] = node; 
	// 	std::cout << sigma->label << ", ";
	// }




	// idx_t (&c1)[3] = *static_cast< idx_t (*)[3]>(static_cast<void*>(tau_labels.data()));

	// do_something(tau_labels, [](auto label){
	// 	std::cout << label << ", "; 
	// });
	// do_something(c1, [](auto label){
	// 	std::cout << label << ", "; 
	// });
	// do_something(as_array_ref< 3, idx_t >(tau_labels), [](auto label){
	// 	std::cout << label << ", "; 
	// });

	// for (auto& el: tau_labels){
	// 	std::cout << el << ", "; 
	// }
	// std::cout << std::endl;

	// 
	// // v_t< 4 > tau2_labels = full_simplex< 4 >(st.find_node(tau2));
	// auto tau2_labels = full_simplex< 4 >(st.find_node(tau2));
	// std::cout << std::endl;
	// v_t< 20 > = 


	// idx_t (&c2)[4] = *static_cast< idx_t (*)[4]>(static_cast<void*>(tau2_labels.data()));

	// // std::cout << tau2_labels.size(); 
	// do_something(tau2_labels, [](auto label){
	// 	std::cout << label << ", "; 
	// });
	// do_something(c2, [](auto label){
	// 	std::cout << label << ", "; 
	// });

	// std::tuple_size_v< v_t< 4 > >
	// full_simplex3< 3 >(cn);
	


	// using d_node = std::pair< SimplexTree::node_ptr, idx_t >;
	// using d_node = std::tuple< SimplexTree::node_ptr, idx_t >;
	// using node_t = typename dfs< false >::t_node;
	// constexpr auto always_true = [](node_t& cn) -> bool { return true; };
	// const auto wut = [](node_t& cn) -> bool { return true; };
	// const auto is_edge = [](node_t& np) -> bool { return(std::get< 1 >(np) == 2); };
	// const auto depth_bound = [](node_t& np) -> bool { return(std::get< 1 >(np) <= 1); };

	// std::cout << "tree depth: " <<  st.tree_max_depth << std::endl; 
	

	// auto dfs_traversal = dfs< false >(&st, st.root.get(), is_edge, depth_bound);
	// for (auto& node: dfs_traversal){
	// 	// dfs_traversal.simplex;
	// 	// auto sigma = st.full_simplex(np.first);
	// 	auto [np, d] = node; 
	// 	//auto sigma = st.full_simplex(np);
	// 	for (auto& tau: sigma){
	// 		std::cout << tau << ", " << std::flush;
	// 	}
	// 	std::cout << std::endl; 
	// }

	// // st.traverse_skeleton(st.root.get(), [](node_ptr cn, size_t d){
	// // 	std::cout << cn->label << ","; 
	// // }, 1);

	// // TODO: abstract traversal on R side to accept as first argument list of options + function


	// st.remove_simplex( {std::size_t(1), std::size_t(2)} );
	// st.expansion(2);


	// auto ns = st.n_simplexes;
	// std::cout << std::endl;
	// for (auto s_dim: ns){
	// 	std::cout << s_dim << ", ";
	// }
	// std::cout << std::endl;
	return 0; 
}