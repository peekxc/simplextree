#include <Rcpp.h>
using namespace Rcpp; 
// [[Rcpp::plugins(cpp17)]] 

#include "simplextree.h"

SEXP as_XPtr(SimplexTree* st){
  Rcpp::XPtr< SimplexTree > p(st, false);
  return(p);
}


// void faces(SimplexTree* st, vector< idx_t > sigma){
//   st->traverse_faces_s(st->find_node(sigma), [](SimplexTree::simplex_t tau){
//     IntegerVector tau_iv = wrap(tau);
//     Rcout << tau_iv << std::endl; 
//   });
// }

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
void copy_trees(SEXP st1, SEXP st2){
  Rcpp::XPtr<SimplexTree> st1_ptr(st1), st2_ptr(st2);
  *st2_ptr = static_cast< const SimplexTree& >(*st1_ptr);
}

// // Traversal overload 1
// void traverse_1(SimplexTree* st, Function f, std::string type){
//   st->traverse_int(R_NilValue, f, type, R_NilValue, false);
// }
// // Traversal overload 2
// void traverse_2(SimplexTree* st, SEXP simp, Function f, std::string type){
//   st->traverse_int(simp, f, type, R_NilValue, false);
// }
// // Traversal overload 3
// void traverse_3(SimplexTree* st, SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args) {
//   st->traverse_int(simp, f, type, args, false);
// }
// 
// // List-traversal overload 1
// List ltraverse_1(SimplexTree* st, Function f, std::string type){
//   return st->traverse_int(R_NilValue, f, type, R_NilValue, true);
// }
// // List-traversal overload 2
// List ltraverse_2(SimplexTree* st, SEXP simp, Function f, std::string type){
//   return st->traverse_int(simp, f, type, R_NilValue, true);
// }
// // List-traversal overload 3
// List ltraverse_3(SimplexTree* st, SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args) {
//   return st->traverse_int(simp, f, type, args, true);
// }
// 
// // Vector-traversal overload 1
// SEXP straverse_1(SimplexTree* st, Function f, std::string type) {
//   Environment base("package:base"); 
//   Function s2arr = base["simplify2array"];
//   return s2arr(st->traverse_int(R_NilValue, f, type, R_NilValue, true));
// }
// // Vector-traversal overload 2
// SEXP straverse_2(SimplexTree* st, SEXP simp, Function f, std::string type){
//   Environment base("package:base"); 
//   Function s2arr = base["simplify2array"];
//   return s2arr(st->traverse_int(simp, f, type, R_NilValue, true));
// }
// // Vector-traversal overload 3
// SEXP straverse_3(SimplexTree* st, SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args) {
//   Environment base("package:base"); 
//   Function s2arr = base["simplify2array"];
//   return s2arr(st->traverse_int(simp, f, type, args, true));
// }

// Workaround for internal const-member overloades
int degree_1(SimplexTree* st, int id){
  return st->degree(static_cast< idx_t >(id));
}

IntegerVector degree_2(SimplexTree* st, vector< idx_t > ids){
  return wrap(st->degree(ids));
}


#include <Rcpp/Benchmark/Timer.h>

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

// template <typename Lambda>
// inline void SimplexTree::trav_switch(node_ptr sigma, Lambda f, std::string type, List args) const{
//   if (type == "dfs") { traverse_dfs(sigma, f); }
//   else if (type == "bfs") { traverse_bfs(sigma, f); }
//   else if (type == "cofaces" || type == "star") { traverse_cofaces(sigma, f); } 
//   else if (type == "link"){ traverse_link(sigma, f); } 
//   else if (type == "skeleton" || type == "maximal-skeleton"){
//     if (args.size() == 0){ stop("Expecting dimension 'k' to be passed."); }
//     vector< std::string > args_str = as< vector< std::string > >(args.names());
//     bool contains_k = std::any_of(args_str.begin(), args_str.end(), [](const std::string arg){
//       return arg == "k";
//     });
//     if (!contains_k){ stop("Expecting dimension 'k' to be passed."); }
//     size_t k = args["k"];
//     if (type == "skeleton"){ traverse_skeleton(sigma, f, k); }
//     if (type == "maximal-skeleton"){ traverse_max_skeleton(sigma, f, k); }
//   } else if(type == "facets"){
//     traverse_facets(sigma, f);
//   } else { stop("Iteration 'type' is invalid. Please use one of: dfs, bfs, cofaces, star, link, skeleton, or maximal-skeleton"); }
// }

// Generic way to apply function to various types of simplices. 
// This acts the generic R-facing version.
// inline List SimplexTree::traverse_int(SEXP simp, Function f, std::string type, Rcpp::Nullable<List> args, bool save_res) const{
//   node_ptr sigma = bool(Rf_isNull(simp)) ? root : find_node( as< vector<idx_t> >(simp) );
//   if (sigma == nullptr){ return List::create(); }
//   
//   // Get the arguments
//   List actual_args = args.isNotNull() ? List(args) : List();
//   // vector< std::string > args_str;as< vector< std::string > >(actual_args.names());
//   
//   // Store results in a list
//   List res = List(); 
//   
//   // Default eval function: retrieve the full simplex, and call R function with that simplex
//   auto default_f = [this, &f, &res](const node_ptr cn, const size_t d = 0){
//     vector< idx_t > simplex = full_simplex(cn);
//     res.push_back(f(wrap(simplex)));
//   };
//   
//   // Do the traversal, possibly saving the results
//   trav_switch(sigma, default_f, type, actual_args);
//   return res; 
// }

// Exports the 1-skeleton as an adjacency matrix 
// inline IntegerMatrix SimplexTree::as_adjacency_matrix() const{
//   const size_t n = root->children.size();
//   IntegerMatrix res = IntegerMatrix(n, n);
//   
//   // Fill in the adjacency matrix
//   for (const node_ptr& vi: root->children){
//     const size_t i = vertex_index(vi->label);
//     for (const node_ptr& vj: vi->children){
//       const size_t j = vertex_index(vj->label);
//       res.at(i, j) = res.at(j, i) = 1; 
//     }
//   }
//   return(res);
// }

// // Exports the 1-skeleton as an adjacency list 
// inline List SimplexTree::as_adjacency_list() const{
//   const size_t n = root->children.size();
//   std::unordered_map< std::string, vector< idx_t > > res(n);
//   for (const node_ptr& vi: root->children){
//     for (const node_ptr& vj: vi->children){
//       res[std::to_string(vi->label)].push_back(vj->label);
//       res[std::to_string(vj->label)].push_back(vi->label);
//     }
//   }
//   return(wrap(res));
// }
// 
// // Exports the 1-skeleton as an edgelist 
// inline IntegerMatrix as_edge_list(SimplexTree* st) const{
//   return get_k_simplices(1);
// }
// 
// // Exports the k-skeleton as a list
// inline List as_list(SimplexTree* st) const{
//   List res = List();
//   vector<idx_t> all = vector<idx_t>(); 
//   idx_t d = 1;
//   std::for_each(begin_bfs(root), end_bfs(), [this, &d,&res, &all](const node_ptr sigma){
//     if ( sigma != nullptr && sigma != root){
//       vector<idx_t> si = full_simplex(sigma);
//       if (si.size() > d){ 
//         const size_t n = all.size() / d;
//         IntegerMatrix tmp = IntegerMatrix(n, d);
//         for (size_t i = 0; i < n; ++i){
//           IntegerVector row = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
//           tmp(i, _) = row;
//         }
//         res.push_back(tmp);
//         all.clear(); 
//         d = si.size(); 
//       }
//       all.insert(all.end(), si.begin(), si.end());
//     }
//   });
//   const size_t n = all.size() / d;
//   IntegerMatrix tmp = IntegerMatrix(n, d);
//   for (size_t i = 0; i < n; ++i){
//     IntegerVector row = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
//     tmp(i, _) = row;
//   }
//   res.push_back(tmp);
//   return res;
// }

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  using namespace Rcpp;
  Rcpp::class_<SimplexTree>("SimplexTree")
    .constructor()
    .method( "as_XPtr", &as_XPtr)
    .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    // .property("dimension", &SimplexTree::dimension)
    // .property("rips_index", &SimplexTree::rips_index)
    // .property("rips_epsilon", &SimplexTree::rips_epsilon)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &SimplexTree::get_vertices)
    // .property("edges", &SimplexTree::get_edges)
    // .property("triangles", &SimplexTree::get_triangles)
    // .property("quads", &SimplexTree::get_quads)
    .property("connected_components", &SimplexTree::connected_components)
    // .property( "rips_simplices", &SimplexTree::rips_simplices)
    // .property( "rips_weights", &SimplexTree::rips_weights)
    // .method("print_cousins", &SimplexTree::get_cousins)
    .const_method( "print_tree", &SimplexTree::print_tree )
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
    // .const_method( "is_face", &SimplexTree::is_face)
    // .const_method( "is_tree", &SimplexTree::is_tree)
    // .method( "traverse", &traverse_1)
    // .method( "traverse", &traverse_2)
    // .method( "traverse", &traverse_3)
    // .method( "ltraverse", &ltraverse_1)
    // .method( "ltraverse", &ltraverse_2)
    // .method( "ltraverse", &ltraverse_3)
    // .method( "straverse", &straverse_1)
    // .method( "straverse", &straverse_2)
    // .method( "straverse", &straverse_3)
    // .const_method( "as_adjacency_matrix", &SimplexTree::as_adjacency_matrix)
    // .const_method( "as_adjacency_list", &SimplexTree::as_adjacency_list)
    // .const_method( "as_edge_list", &SimplexTree::as_edge_list)
    // .const_method( "as_list", &SimplexTree::as_list)
    // .const_method( "serialize", &SimplexTree::serialize)
    // .method( "deserialize", &SimplexTree::deserialize)
    // .const_method( "save", &SimplexTree::save)
    // .method( "load", &SimplexTree::load)
    // .method( "faces", &faces)
    // .method( "rips", &rips)
    // .method( "threshold", &threshold)
    // .method( "threshold_index", &SimplexTree::threshold_index)
    // .method( "threshold_function", &SimplexTree::threshold_function)
    // .method( "print_filtration", &print_filtration )
    ;
}





