#include "simplextree.h"

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  using namespace Rcpp;
  Rcpp::class_<SimplexTree>("SimplexTree")
    .constructor()
    .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    .property("dimension", &SimplexTree::dimension)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &SimplexTree::get_vertices)
    .property("edges", &SimplexTree::get_edges)
    .property("triangles", &SimplexTree::get_triangles)
    .property("quads", &SimplexTree::get_quads)
    // .method("print_cousins", &SimplexTree::get_cousins)
    .method( "print_tree", &SimplexTree::print_tree )
    .method( "as_XPtr", &SimplexTree::as_XPtr)
    .method( "clear", &SimplexTree::clear)
    .method( "generate_ids", &SimplexTree::generate_ids)
    .method( "degree", &SimplexTree::degree)
    .method( "adjacent", &SimplexTree::adjacent_vertices)
    .method( "insert",  (void (SimplexTree::*)(SEXP))(&SimplexTree::insert))
    .method( "remove",  (void (SimplexTree::*)(SEXP))(&SimplexTree::remove))
    .method( "find",  (LogicalVector (SimplexTree::*)(SEXP))(&SimplexTree::find))
    // .method( "expansion", &SimplexTree::expansion )
    .method( "collapse", &SimplexTree::collapseR)
    .method( "collapse", &SimplexTree::vertex_collapseR)
    .method( "contract", &SimplexTree::contract)
    .method( "is_face", &SimplexTree::is_face)
    .method( "is_tree", &SimplexTree::is_tree)
    .method( "traverse", (void (SimplexTree::*)(Function, std::string))(&SimplexTree::traverse))
    .method( "traverse", (void (SimplexTree::*)(SEXP, Function, std::string))(&SimplexTree::traverse))
    .method( "traverse", (void (SimplexTree::*)(SEXP, Function, std::string, Rcpp::Nullable<List>))(&SimplexTree::traverse))
    .method( "ltraverse", (List (SimplexTree::*)(Function, std::string))(&SimplexTree::ltraverse))
    .method( "ltraverse", (List (SimplexTree::*)(SEXP, Function, std::string))(&SimplexTree::ltraverse))
    .method( "ltraverse", (List (SimplexTree::*)(SEXP, Function, std::string, Rcpp::Nullable<List>))(&SimplexTree::ltraverse))
    .method( "as_adjacency_matrix", &SimplexTree::as_adjacency_matrix)
    .method( "as_adjacency_list", &SimplexTree::as_adjacency_list)
    .method( "as_edge_list", &SimplexTree::as_edge_list)
    .method( "as_list", &SimplexTree::as_list)
    .method( "serialize", &SimplexTree::serialize)
    .method( "deserialize", &SimplexTree::deserialize)
    .method( "save", &SimplexTree::save)
    .method( "load", &SimplexTree::load)
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







