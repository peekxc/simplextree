#include "simplextree.h"

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  using namespace Rcpp;
  Rcpp::class_<SimplexTree>("SimplexTree")
    .constructor()
    .field_readonly("n_simplexes", &SimplexTree::n_simplexes)
    .field_readonly("max_depth", &SimplexTree::tree_max_depth)
    .method( "as_XPtr", &SimplexTree::as_XPtr)
  // .method( "add_vertices", &SimplexTree::add_vertices)
  // .method( "remove_vertices", &SimplexTree::remove_vertices)
  // .method( "vertex_available", &SimplexTree::vertex_available)
     .method( "adjacent_vertices", &SimplexTree::adjacent_vertices)
     .method( "insert_simplex", &SimplexTree::insert_simplex)
     .method( "find_simplex", &SimplexTree::find_simplex)
  // .method( "remove_edge", &SimplexTree::remove_edge)
     .method( "remove_simplex", &SimplexTree::remove_simplex)
     .method( "print_tree", &SimplexTree::print_tree )
  // .method( "print_cofaces", &SimplexTree::print_cofaces )
  // .method( "print_cousins", &SimplexTree::print_cousins )
  // .method( "expansion", &SimplexTree::expansion )
     .method( "as_adjacency_matrix", &SimplexTree::as_adjacency_matrix )
     .method( "as_adjacency_list", &SimplexTree::as_adjacency_list)
     .method( "as_edge_list", &SimplexTree::as_edge_list)
     .method( "as_list", &SimplexTree::as_list)
  // .method( "get_vertices", &SimplexTree::get_vertices)
  // .method( "apply_dfs", &SimplexTree::apply_dfs)
     .method( "collapse", &SimplexTree::collapseR)
     .method( "contract", &SimplexTree::contract)
     .method( "is_face", &SimplexTree::is_face)
  // .method( "remove_subtree", &SimplexTree::remove_subtree2)
     .method( "apply", &SimplexTree::apply)
     .method( "serialize", &SimplexTree::serialize)
     .method( "unserialize", &SimplexTree::unserialize)
  // .method( "test", &SimplexTree::test)
  ;
}

/***R
## Inline tests
stree <- simplextree::simplex_tree()
stree$insert_simplex(c(1, 2, 3))
stree$insert_simplex(c(2, 3, 4, 5))

stree$apply(1L, print, "link")

stree$apply( c(1, 2), function(x){ print(x) }, "bfs")
stree$apply(c(1), function(x){ print(x) }, "cofaces")
stree$remove_simplex(c(2, 3))

stree$remove_simplex(c(1, 2, 3))
stree$print_cofaces(c(1, 2, 3))

## Testing collapses
stree <- simplextree::simplex_tree()
stree$insert_simplex(c(1, 2))
stree$test(1)
stree$apply(1L, print, "cofaces")

stree <- simplextree::simplex_tree()
stree$insert_simplex(c(1, 2, 3))
stree$collapse(c(1, 2), c(1, 2, 3))
# stree$test()

stree <- simplextree::simplex_tree()
stree$insert_simplex(c(1, 2, 3))
stree$print_tree()
stree$contract(c(1, 3))
stree$print_tree()

stree$insert_simplex(1)
stree$insert_simplex(2)
stree$insert_simplex(3)
stree$insert_simplex(c(1, 2))
stree$insert_simplex(c(2, 3))
stree$insert_simplex(c(1, 3))
stree$insert_simplex(c(1, 2, 3))


stree <- simplex_tree()
stree$insert_simplex(1:3)
stree$insert_simplex(2:5)
stree$insert_simplex(5:9)
stree$insert_simplex(7:8)
stree$insert_simplex(10)
stree$print_tree()
# 1 (h = 2): .( 2 3 )..( 3 )
# 2 (h = 3): .( 3 4 5 )..( 4 5 5 )...( 5 )
# 3 (h = 2): .( 4 5 )..( 5 )
# 4 (h = 1): .( 5 )
# 5 (h = 4): .( 6 7 8 9 )..( 7 8 9 8 9 9 )...( 8 9 9 9 )....( 9 )
# 6 (h = 3): .( 7 8 9 )..( 8 9 9 )...( 9 )
# 7 (h = 2): .( 8 9 )..( 9 )
# 8 (h = 1): .( 9 )
# 9 (h = 0): 
# 10 (h = 0): 

## Check cofaces
stree$apply(tau, print, "cofaces")
# [1] 2 3 4 5
# [1] 3 4 5

stree$collapse(tau, sigma)
# [1] TRUE

# 1 (h = 2): .( 2 3 )..( 3 )
# 2 (h = 2): .( 3 4 5 )..( 4 5 5 )
# 3 (h = 1): .( 4 5 )
# 4 (h = 1): .( 5 )
# 5 (h = 4): .( 6 7 8 9 )..( 7 8 9 8 9 9 )...( 8 9 9 9 )....( 9 )
# 6 (h = 3): .( 7 8 9 )..( 8 9 9 )...( 9 )
# 7 (h = 2): .( 8 9 )..( 9 )
# 8 (h = 1): .( 9 )
# 9 (h = 0): 
# 10 (h = 0): 
*/







