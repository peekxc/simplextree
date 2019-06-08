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
    .property("connected_components", &SimplexTree::connected_components)
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
    .method( "expand", &SimplexTree::expansion )
    .method( "collapse", &SimplexTree::collapseR)
    .method( "collapse", &SimplexTree::vertex_collapseR)
    .method( "contract", &SimplexTree::contract)
    .method( "reindex", (void (SimplexTree::*)(SEXP))(&SimplexTree::reindex))
    .method( "traverse", (void (SimplexTree::*)(Function, std::string))(&SimplexTree::traverse))
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





