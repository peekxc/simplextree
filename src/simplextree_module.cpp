#include "simplextree.h"
#include <Rcpp.h>
using namespace Rcpp; 

void faces(SimplexTree* st, vector< idx_t > sigma){
  st->traverse_faces_s(st->find_node(sigma), [](SimplexTree::simplex_t tau){
    IntegerVector tau_iv = wrap(tau);
    Rcout << tau_iv << std::endl; 
  });
}

// Creates a filtration from a neighborhood graph
void rips(SimplexTree* st, vector< double > edge_weights, const size_t k){
  st->rips(edge_weights, k);
}

// double filtration_index(SimplexTree* st){

// double filtration_index(SEXP st){
//   Rcpp::XPtr< SimplexTree > stree_ptr(st);
//   return filtration_index(static_cast< SimplexTree* >(stree_ptr));
// }

// For Debugging
void print_filtration(SimplexTree* st){
  Rcout << "Included vector: " << std::endl;
  LogicalVector is_included = wrap(st->included);
  Rcout << is_included << std::endl;
}

// Workaround to allow returning simplextree references when modify==TRUE 
// SEXP threshold(SimplexTree* st, std::string type, double eps_or_idx, bool modify){
//   if (type != "index" && type != "function"){
//     stop("Threshold must be given 'function' or 'index' as the threshold type.");
//   }
//   // Obtaining namespace of Matrix package
//   Environment pkg = Environment::namespace_env("simplextree");
//   Function f = pkg[".threshold_new"];
//   
//   // If mody is true, just modify the current tree based on the input
//   if (modify){
//     if (type == "index"){ st->threshold_index(static_cast< size_t >(eps_or_idx)); }
//     else if (type == "function"){ st->threshold_function(static_cast< size_t >(eps_or_idx)); }
//   } 
//   return f(st->as_XPtr(), type, eps_or_idx, modify);
// }

void insert_matrix(SimplexTree* st, NumericMatrix x){
  Rcout << "inserting matrix" << std::endl; 
  const size_t n = x.nrow();
  for (size_t i = 0; i < n; ++i){
    NumericVector cr = x(i,_);
    st->insert_simplex(as< vector< idx_t > >(cr));
  }
}
  
// Copies the contents of st1 to st2
// [[Rcpp::export]]
void copy_trees(SEXP st1, SEXP st2){
  Rcpp::XPtr<SimplexTree> st1_ptr(st1), st2_ptr(st2);
  *st2_ptr = static_cast< const SimplexTree& >(*st1_ptr);
}

// Traversal overload 1
void traverse_1(SimplexTree* st, Function f, std::string type){
  st->traverse_int(R_NilValue, f, type, R_NilValue, false);
}
// Traversal overload 2
void traverse_2(SimplexTree* st, SEXP simp, Function f, std::string type){
  st->traverse_int(simp, f, type, R_NilValue, false);
}
// Traversal overload 3
void traverse_3(SimplexTree* st, SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args) {
  st->traverse_int(simp, f, type, args, false);
}

// List-traversal overload 1
List ltraverse_1(SimplexTree* st, Function f, std::string type){
  return st->traverse_int(R_NilValue, f, type, R_NilValue, true);
}
// List-traversal overload 2
List ltraverse_2(SimplexTree* st, SEXP simp, Function f, std::string type){
  return st->traverse_int(simp, f, type, R_NilValue, true);
}
// List-traversal overload 3
List ltraverse_3(SimplexTree* st, SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args) {
  return st->traverse_int(simp, f, type, args, true);
}

// Vector-traversal overload 1
SEXP straverse_1(SimplexTree* st, Function f, std::string type) {
  Environment base("package:base"); 
  Function s2arr = base["simplify2array"];
  return s2arr(st->traverse_int(R_NilValue, f, type, R_NilValue, true));
}
// Vector-traversal overload 2
SEXP straverse_2(SimplexTree* st, SEXP simp, Function f, std::string type){
  Environment base("package:base"); 
  Function s2arr = base["simplify2array"];
  return s2arr(st->traverse_int(simp, f, type, R_NilValue, true));
}
// Vector-traversal overload 3
SEXP straverse_3(SimplexTree* st, SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args) {
  Environment base("package:base"); 
  Function s2arr = base["simplify2array"];
  return s2arr(st->traverse_int(simp, f, type, args, true));
}

// Workaround for internal const-member overloades
int degree_1(SimplexTree* st, int id){
  return st->degree(static_cast< idx_t >(id));
}

IntegerVector degree_2(SimplexTree* st, vector< idx_t > ids){
  return wrap(st->degree(ids));
}



// [[Rcpp::export]]
NumericVector profile(SEXP st){
  Rcpp::XPtr< SimplexTree > st_ptr(st);
  Timer timer;
  timer.step("start");
  st_ptr->expansion(2);
  timer.step("expansion");
  
  NumericVector res(timer);
  const size_t n = 1000;
  for (size_t i=0; i < res.size(); ++i) { res[i] = res[i] / n; }
  return res;
}

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  using namespace Rcpp;
  Rcpp::class_<SimplexTree>("SimplexTree")
    .constructor()
    .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    .property("dimension", &SimplexTree::dimension)
    .property("rips_index", &SimplexTree::rips_index)
    .property("rips_epsilon", &SimplexTree::rips_epsilon)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &SimplexTree::get_vertices)
    .property("edges", &SimplexTree::get_edges)
    .property("triangles", &SimplexTree::get_triangles)
    .property("quads", &SimplexTree::get_quads)
    .property("connected_components", &SimplexTree::connected_components)
    .property( "rips_simplices", &SimplexTree::rips_simplices)
    .property( "rips_weights", &SimplexTree::rips_weights)
    // .method("print_cousins", &SimplexTree::get_cousins)
    .const_method( "print_tree", &SimplexTree::print_tree )
    .method( "as_XPtr", &SimplexTree::as_XPtr)
    .method( "clear", &SimplexTree::clear)
    .method( "degree", &degree_1)
    .method( "degree", &degree_2)
    .method( "generate_ids", &SimplexTree::generate_ids)
    .method( "adjacent", &SimplexTree::adjacent_vertices)
    .method( "insert",  (void (SimplexTree::*)(SEXP))(&SimplexTree::insert))
    .method( "insert",  &insert_matrix)
    .method( "remove",  (void (SimplexTree::*)(SEXP))(&SimplexTree::remove))
    .method( "find",  (LogicalVector (SimplexTree::*)(SEXP))(&SimplexTree::find))
    .method( "expand", &SimplexTree::expansion )
    .method( "collapse", &SimplexTree::collapseR)
    .method( "collapse", &SimplexTree::vertex_collapseR)
    .method( "contract", &SimplexTree::contract)
    .method( "reindex", (void (SimplexTree::*)(SEXP))(&SimplexTree::reindex))
    .const_method( "is_face", &SimplexTree::is_face)
    .const_method( "is_tree", &SimplexTree::is_tree)
    .method( "traverse", &traverse_1)
    .method( "traverse", &traverse_2)
    .method( "traverse", &traverse_3)
    .method( "ltraverse", &ltraverse_1)
    .method( "ltraverse", &ltraverse_2)
    .method( "ltraverse", &ltraverse_3)
    .method( "straverse", &straverse_1)
    .method( "straverse", &straverse_2)
    .method( "straverse", &straverse_3)
    .const_method( "as_adjacency_matrix", &SimplexTree::as_adjacency_matrix)
    .const_method( "as_adjacency_list", &SimplexTree::as_adjacency_list)
    .const_method( "as_edge_list", &SimplexTree::as_edge_list)
    .const_method( "as_list", &SimplexTree::as_list)
    .const_method( "serialize", &SimplexTree::serialize)
    .method( "deserialize", &SimplexTree::deserialize)
    .const_method( "save", &SimplexTree::save)
    .method( "load", &SimplexTree::load)
    .method( "faces", &faces)
    .method( "rips", &rips)
    // .method( "threshold", &threshold)
    .method( "threshold_index", &SimplexTree::threshold_index)
    .method( "threshold_function", &SimplexTree::threshold_function)
    // .method( "print_filtration", &print_filtration )
    ;
}





