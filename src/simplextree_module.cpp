#include <Rcpp.h>
using namespace Rcpp; 

// [[Rcpp::plugins(cpp17)]] 
#include "simplextree.h"

using simplex_t = SimplexTree::simplex_t; 

SEXP as_XPtr(SimplexTree* st){
  Rcpp::XPtr< SimplexTree > p(st, false);
  return(p);
}

// Generic function to handle various vector types
template < bool check_rng = false, typename Lambda >
void vector_handler(SEXP sigma, Lambda&& f){
  const unsigned int s_type = TYPEOF(sigma);
  const auto check_valid = [](SEXP v) -> bool { 
    IntegerVector iv = v;
    return std::any_of(begin(iv), end(iv), [](int el) -> bool { 
      return(el < 0 || el > std::numeric_limits< idx_t >::max()); 
    });
  };
  if (!Rf_isNull(Rf_getAttrib(sigma, R_DimSymbol))){
    IntegerMatrix m = as< IntegerMatrix >(sigma);
    const size_t n = m.nrow();
    for (size_t i = 0; i < n; ++i){
      IntegerVector cr = m(i,_);
      if constexpr (check_rng) {
        if (check_valid(sigma)){ stop("Only unsigned integer simplices are supported."); }
      }
      f(as< simplex_t >(cr));
    }
  } else if (s_type == INTSXP || s_type == REALSXP){
    if constexpr (check_rng) {
      if (check_valid(sigma)){ stop("Only unsigned integer simplices are supported."); }
    }
    f(as< simplex_t >(sigma));
  } else if (s_type == LISTSXP || s_type == VECSXP){
    List simplices = List(sigma);
    const size_t n = simplices.size(); 
    for (size_t i = 0; i < n; ++i){
      if constexpr (check_rng) {
        if (check_valid(simplices.at(i))){ stop("Only unsigned integer simplices are supported."); }
      }
      f(as< simplex_t >(simplices[i]));
    }
  } else { stop("Unknown type passed, must be list type or vector type."); }
}

// R-facing Inserter 
void insert(SimplexTree* st, SEXP x, const bool check_overflow=false){
  if (check_overflow){
    vector_handler< true >(x, [st](simplex_t&& sigma){
      st->insert_simplex(sigma);
    });
  } else {
    vector_handler< false >(x, [st](simplex_t&& sigma){
      st->insert_simplex(sigma);
    });
  }
}

// R-facing remover 
void remove(SimplexTree* st, SEXP x, const bool check_overflow=false){
  if (check_overflow){
    vector_handler< true >(x, [st](simplex_t&& sigma){
      st->remove_simplex(sigma);
    });
  } else {
    vector_handler< false >(x, [st](simplex_t&& sigma){
      st->remove_simplex(sigma);
    });
  }
}

// Creates a filtration from a neighborhood graph
// void rips(SimplexTree* st, vector< double > edge_weights, const size_t k){
//   st->rips(edge_weights, k);
// }

// double filtration_index(SimplexTree* st){

// double filtration_index(SEXP st){
//   Rcpp::XPtr< SimplexTree > stree_ptr(st);
//   return filtration_index(static_cast< SimplexTree* >(stree_ptr));
// }

// For Debugging
// void print_filtration(SimplexTree* st){
//   Rcout << "Included vector: " << std::endl;
//   LogicalVector is_included = wrap(st->included);
//   Rcout << is_included << std::endl;
// }

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

// Copies the contents of st1 to st2
void copy_trees(SEXP st1, SEXP st2){
  Rcpp::XPtr<SimplexTree> st1_ptr(st1), st2_ptr(st2);
  *st2_ptr = static_cast< const SimplexTree& >(*st1_ptr);
}


// Workaround for internal const-member overloades
int degree_1(SimplexTree* st, int id){
  return st->degree(static_cast< idx_t >(id));
}

IntegerVector degree_2(SimplexTree* st, vector< idx_t > ids){
  return wrap(st->degree(ids));
}

// inline void SimplexTree::reindex(SEXP target_ids){
//   const unsigned int s_type = TYPEOF(target_ids);
//   if (s_type == INTSXP || s_type == REALSXP){
//     vector< idx_t > t_ids = as< vector< idx_t > >(target_ids);
//     reindex(t_ids);
//   } else if (s_type == LISTSXP || s_type == VECSXP){
//     List t_ids_lst = as< List >(target_ids); 
//     CharacterVector nm = t_ids_lst.names();
//     if (Rf_isNull(nm) || Rf_length(nm) == 0){ stop("target ids must be named if given as a list."); }
//     
//     // Do the mapping, then send to regular reindexing function
//     const vector< idx_t > base_vids = get_vertices();
//     vector< idx_t > new_vids = vector< idx_t >(begin(base_vids), end(base_vids));
//     for (size_t i = 0; i < nm.size(); ++i){
//       idx_t src_id = std::stoi(as< std::string >(nm.at(i)));
//       idx_t tgt_id = as< idx_t >(t_ids_lst.at(i));
//       const size_t idx = std::distance(begin(base_vids), std::lower_bound(begin(base_vids), end(base_vids), src_id));
//       new_vids.at(idx) = tgt_id;
//     }
//     reindex(new_vids);
//   }
// }
// 
inline void SimplexTree::reindex(vector< idx_t > target_ids){
  if (n_simplexes.at(0) != target_ids.size()){ stop("target id vector must match the size of the number of 0-simplices."); }
  vector< vector< idx_t > > minimal = serialize();
  vector< idx_t > vids = get_vertices();

  // Check the target ids are unique
  vector< idx_t > v_check(begin(target_ids), end(target_ids));
  std::sort(begin(v_check), end(v_check));
  auto it = std::unique(begin(v_check), end(v_check));
  if (std::distance(begin(v_check), it) != vids.size()){
    stop("target ids must all unique.");
  }

  clear(); // clear the tree now that it's been serialized
  for (simplex_t sigma: minimal){
    const size_t n = sigma.size();
    for (size_t i = 0; i < n; ++i){
      const size_t idx = std::distance(begin(vids), std::lower_bound(begin(vids), end(vids), sigma.at(i)));
      sigma.at(i) = target_ids.at(idx);
    }
    insert_simplex(sigma);
  }
}

LogicalVector find_R(SimplexTree* st, SEXP simplices){
  LogicalVector v; 
  vector_handler< false >(simplices, [st, &v](simplex_t&& sigma){
    v.push_back(st->find_simplex(sigma));
  });
  return(v);
}


IntegerMatrix get_k_simplices(SimplexTree* st, const size_t k) {
  if (st->n_simplexes.size() <= k){ return IntegerMatrix(0, k+1); }
  IntegerMatrix res = IntegerMatrix(st->n_simplexes.at(k), k+1);
  size_t i = 0; 
  auto tr = st::max_skeleton< true >(st, st->root.get(), k);
  traverse(tr, [&res, &i](node_ptr cn, idx_t depth, simplex_t sigma){
    res(i++, _) = IntegerVector(sigma.begin(), sigma.end());
    return true; 
  });
  return(res);
}

// Retrieve the vertices by their label
vector< idx_t > get_vertices(SimplexTree* st) {
  if (st->n_simplexes.size() == 0){ return vector< idx_t >(); } //IntegerVector(); }
  vector< idx_t > v;
  v.reserve(st->n_simplexes.at(0));
  for (auto& cn: st->root->children){ 
    v.push_back(cn->label); 
  }
  //return IntegerVector(v.begin(), v.end());
  return(v);
}
IntegerMatrix get_edges(SimplexTree* st) { return get_k_simplices(st, 1); }
IntegerMatrix get_triangles(SimplexTree* st) { return get_k_simplices(st, 2); }
IntegerMatrix get_quads(SimplexTree* st) { return get_k_simplices(st, 3); }


// Exports the 1-skeleton as an adjacency matrix 
IntegerMatrix as_adjacency_matrix(SimplexTree* st) {
  const size_t n = st->root->children.size();
  IntegerMatrix res = IntegerMatrix(n, n);

  // Fill in the adjacency matrix
  for (auto& vi: st->root->children){
    const size_t i = st->vertex_index(vi->label);
    for (auto& vj: vi->children){
      const size_t j = st->vertex_index(vj->label);
      res.at(i, j) = res.at(j, i) = 1;
    }
  }
  return(res);
}

// Exports the 1-skeleton as an adjacency list
List as_adjacency_list(SimplexTree* st) {
  const size_t n = st->root->children.size();
  std::unordered_map< std::string, vector< idx_t > > res(n);
  for (auto& vi: st->root->children){
    for (auto& vj: vi->children){
      res[std::to_string(vi->label)].push_back(vj->label);
      res[std::to_string(vj->label)].push_back(vi->label);
    }
  }
  return(wrap(res));
}

// Exports the 1-skeleton as an edgelist
IntegerMatrix as_edge_list(SimplexTree* st) {
  return get_k_simplices(st, 1);
}

// Exports the k-skeleton as a list
List as_list(SimplexTree* st){
  List res = List();
  vector< idx_t > all = vector< idx_t >();
  idx_t d = 1;
  auto bfs = st::level_order< true >(st);
  traverse(bfs, [&res, &d, &all](node_ptr cn, idx_t depth, simplex_t sigma){
    if (depth > d){
      const size_t n = all.size() / d;
      IntegerMatrix tmp = IntegerMatrix(n, d);
      for (size_t i = 0; i < n; ++i){
        IntegerVector row = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
        tmp(i, _) = row;
      }
      res.push_back(tmp);
      all.clear();
      d = sigma.size();
    }
    all.insert(all.end(), sigma.begin(), sigma.end());
    return true; 
  });
  const size_t n = all.size() / d;
  IntegerMatrix tmp = IntegerMatrix(n, d);
  for (size_t i = 0; i < n; ++i){
    IntegerVector row = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
    tmp(i, _) = row;
  }
  res.push_back(tmp);
  return res;
}

void load(SimplexTree* st, std::string filename){
  Function readRDS = Function("readRDS");
  List st_lst = readRDS(_["file"] = filename);
  const size_t n = st_lst.size();
  for (size_t i = 0; i < n; ++i){
    IntegerVector si = st_lst.at(i);
    simplex_t sigma(si.begin(), si.end());
    st->insert_simplex(sigma);
  }
}

void save(SimplexTree* st, std::string filename){
  using simplex_t = vector< idx_t >;
  Function saveRDS = Function("saveRDS");
  vector< simplex_t > minimal = st->serialize();
  List res = wrap(minimal);
  saveRDS(_["object"] = res, _["file"] = filename);
}
// 
// void print_cousins(SimplexTree* st){
//   st->print_cousins();
// }

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  using namespace Rcpp;
  Rcpp::class_<SimplexTree>("SimplexTree")
    .constructor()
    .method( "as_XPtr", &as_XPtr)
    .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    .property("dimension", &SimplexTree::dimension)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &get_vertices, "Returns the vertex labels as an integer vector.")
    .property("edges", &get_edges, "Returns the edges as an integer matrix.")
    .property("triangles", &get_triangles, "Returns the 2-simplices as an integer matrix.")
    .property("quads", &get_quads, "Returns the 3-simplices as an integer matrix.")
    .property("connected_components", &SimplexTree::connected_components)
    .const_method( "print_tree", &SimplexTree::print_tree )
    .const_method( "print_cousins", &SimplexTree::print_cousins )
    .method( "clear", &SimplexTree::clear)
    .method( "degree", &degree_1)
    .method( "degree", &degree_2)
    .method( "generate_ids", &SimplexTree::generate_ids)
    .method( "adjacent", &SimplexTree::adjacent_vertices)
    .method( "insert",  &insert)
    .method( "remove",  &remove)
    .method( "find", &find_R)
    .method( "expand", &SimplexTree::expansion )
    .method( "collapse", &SimplexTree::collapse)
    .method( "collapse", &SimplexTree::vertex_collapseR)
    .method( "contract", &SimplexTree::contract)
    .const_method( "is_tree", &SimplexTree::is_tree)
    .method( "as_adjacency_matrix", &as_adjacency_matrix)
    .method( "as_adjacency_list", &as_adjacency_list)
    .method( "as_edge_list", &as_edge_list)
    .method( "as_list", &as_list)
    .const_method( "serialize", &SimplexTree::serialize)
    .method( "deserialize", &SimplexTree::deserialize)
    .method( "save", &save)
    .method( "load", &load)
    ;
}

// typedef bool (*ValidConstructor)(SEXP*,int);

void init_filtration(Filtration* filt, SEXP st) {
  Rcpp::XPtr< SimplexTree > st_ptr(st);
  filt->initialize(*st_ptr);
}

List get_simplices(Filtration* st){
  auto si = st->simplices(); 
  List results = List();
  for (auto& sigma: si){ results.push_back(sigma); }
  return results;
}

void make_filtration(Filtration* st, const NumericVector& D){
  if (st->n_simplexes.size() <= 1){ return; }
  const size_t ne = st->n_simplexes.at(1);
  const auto v = st->get_vertices();
  const size_t N = BinomialCoefficient(v.size(), 2);
  if (ne == D.size()){
    vector< double > weights(D.begin(), D.end());
    st->flag(weights, false);
  } else if (D.size() == N){ // full distance vector passed in
    auto edge_iter = st::max_skeleton< true >(st, st->root.get(), 1);
    vector< double > weights;
    weights.reserve(ne);
    st::traverse(edge_iter, [&weights, &D, &v](node_ptr np, idx_t depth, simplex_t sigma){
      auto v1 = sigma[0], v2 = sigma[1];
      auto idx1 = std::distance(begin(v), std::lower_bound(begin(v), end(v), v1));
      auto idx2 = std::distance(begin(v), std::lower_bound(begin(v), end(v), v2));
      auto dist_idx = to_natural_2(idx1, idx2, v.size());
      weights.push_back(D[dist_idx]);
      Rcout << D[dist_idx] << ", ";
      return true; 
    });
    Rcout << std::endl;
    st->flag(weights, false);
  } else {
    throw std::invalid_argument("Flag filtrations require a vector of distances for each edge or a 'dist' object");
  }
}
// 
// void test_filtration(Filtration* st, const size_t i){
//   auto ind = st->simplex_idx(i);
//   for (auto idx: ind){ Rcout << idx << ", "; }
//   Rcout << std::endl;
//   auto sigma = st->expand_simplex(begin(ind), end(ind));
//   for (auto& label: sigma){ Rcout << label << ", "; }
//   Rcout << std::endl; 
// }



RCPP_MODULE(filtration_module) {
  using namespace Rcpp;
  Rcpp::class_< SimplexTree >("SimplexTree")
    .constructor()
    .method( "as_XPtr", &as_XPtr)
    .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    .property("dimension", &SimplexTree::dimension)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &get_vertices, "Returns the vertex labels as an integer vector.")
    .property("edges", &get_edges, "Returns the edges as an integer matrix.")
    .property("triangles", &get_triangles, "Returns the 2-simplices as an integer matrix.")
    .property("quads", &get_quads, "Returns the 3-simplices as an integer matrix.")
    .property("connected_components", &SimplexTree::connected_components)
    .const_method( "print_tree", &SimplexTree::print_tree )
    .const_method( "print_cousins", &SimplexTree::print_cousins )
    .method( "clear", &SimplexTree::clear)
    .method( "degree", &degree_1)
    .method( "degree", &degree_2)
    .method( "generate_ids", &SimplexTree::generate_ids)
    .method( "adjacent", &SimplexTree::adjacent_vertices)
    .method( "insert", &insert)
    .method( "remove", &remove)
    .method( "find", &find_R)
    .method( "expand", &SimplexTree::expansion )
    .method( "collapse", &SimplexTree::collapse)
    .method( "collapse", &SimplexTree::vertex_collapseR)
    .method( "contract", &SimplexTree::contract)
    .const_method( "is_tree", &SimplexTree::is_tree)
    .const_method( "serialize", &SimplexTree::serialize)
    .method( "deserialize", &SimplexTree::deserialize)
    .method( "as_adjacency_matrix", &as_adjacency_matrix)
    .method( "as_adjacency_list", &as_adjacency_list)
    .method( "as_edge_list", &as_edge_list)
    .method( "as_list", &as_list)
    .method( "save", &save)
    .method( "load", &load)
    ;
  Rcpp::class_< Filtration >("Filtration")
    .derives< SimplexTree >("SimplexTree")
    .constructor()
    .method("init_tree", &init_filtration)
    .field("included", &Filtration::included)
    .property("current_index", &Filtration::current_index)
    .property("current_value", &Filtration::current_value)
    .property("simplices", &get_simplices, "Returns the simplices in the filtration")
    .property("weights", &Filtration::weights, "Returns the weights in the filtration")
    .property("dimensions", &Filtration::dimensions, "Returns the dimensions of the simplices in the filtration")
    .method("flag", &make_filtration, "Constructs a flag filtration")
    .method("threshold_value", &Filtration::threshold_value)
    .method("threshold_index", &Filtration::threshold_index)
    ;
}

// --- Begin functional exports + helpers ---

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


bool contains_arg(vector< std::string > v, std::string arg_name){
  return(std::any_of(v.begin(), v.end(), [&arg_name](const std::string arg){ 
      return arg == arg_name;
  }));
};

// From: https://stackoverflow.com/questions/56465550/how-to-concatenate-lists-in-rcpp
List cLists(List x, List y) {
  int nsize = x.size(); 
  int msize = y.size(); 
  List out(nsize + msize);

  CharacterVector xnames = x.names();
  CharacterVector ynames = y.names();
  CharacterVector outnames(nsize + msize);
  out.attr("names") = outnames;
  for(int i = 0; i < nsize; i++) {
    out[i] = x[i];
    outnames[i] = xnames[i];
  }
  for(int i = 0; i < msize; i++) {
    out[nsize+i] = y[i];
    outnames[nsize+i] = ynames[i];
  }
  return(out);
}

// The types of traversal supported
enum TRAVERSAL_TYPE { 
  PREORDER, LEVEL_ORDER, FACES, COFACES, COFACE_ROOTS, K_SKELETON, MAX_SKELETON, MAXIMAL, LINK
}; 
const size_t N_TRAVERSALS = 9;

// Exports a list with the parameters for a preorder traversal
// [[Rcpp::export]]
List parameterize_R(SEXP st, IntegerVector sigma, std::string type, Rcpp::Nullable<List> args){
  List init_params = List::create(_[".ptr"] = st, _["sigma"] = sigma); 
  List param_res = args.isNotNull() ? cLists(List(args), init_params) : init_params;
  if (type == "preorder" || type == "dfs") { param_res["traversal_type"] = int(PREORDER); }
  else if (type == "level_order" || type == "bfs") { param_res["traversal_type"] = int(LEVEL_ORDER); }
  else if (type == "cofaces" || type == "star") { param_res["traversal_type"] = int(COFACES); }
  else if (type == "coface_roots") { param_res["traversal_type"] = int(COFACE_ROOTS); }
  else if (type == "link"){ param_res["traversal_type"] = int(LINK); }
  else if (type == "k_skeleton" || type == "skeleton"){ param_res["traversal_type"] = int(K_SKELETON); }
  else if (type == "k_simplices" || type == "maximal-skeleton"){ param_res["traversal_type"] = int(MAX_SKELETON); }
  else if (type == "maximal"){ param_res["traversal_type"] = int(MAXIMAL); }
  else if(type == "faces"){ param_res["traversal_type"] = int(FACES); }
  else { stop("Iteration 'type' is invalid. Please use one of: preorder, level_order, faces, cofaces, star, link, skeleton, or maximal-skeleton"); }
  param_res.attr("class") = "st_traversal";
  return(param_res);
}

using param_pack = typename std::tuple< SimplexTree*, node_ptr, TRAVERSAL_TYPE >;

template < class Lambda > 
void traverse_switch(param_pack&& pp, List args, Lambda&& f){
  auto args_str = as< vector< std::string > >(args.names());
  SimplexTree* st = get< 0 >(pp);
  node_ptr init = get< 1 >(pp);
  TRAVERSAL_TYPE tt = get< 2 >(pp);
  switch(tt){
    case PREORDER: {
      auto tr = st::preorder< true >(st, init);
      traverse(tr, f);
      break; 
    }
    case LEVEL_ORDER: {
      auto tr = st::level_order< true >(st, init);
      traverse(tr, f);
      break; 
    }
    case FACES: {
      auto tr = st::faces< true >(st, init);
      traverse(tr, f);
      break; 
    }
    case COFACES: {
      auto tr = st::cofaces< true >(st, init);
      traverse(tr, f);
      break; 
    }
    case COFACE_ROOTS: {
      auto tr = st::coface_roots< true >(st, init);
      traverse(tr, f);
      break; 
    }
    case K_SKELETON: {
      if (!contains_arg(args_str, "k")){ stop("Expecting dimension 'k' to be passed."); }
      idx_t k = args["k"];
      auto tr = st::k_skeleton< true >(st, init, k);
      traverse(tr, f);
      break; 
    }
    case MAX_SKELETON: {
      if (!contains_arg(args_str, "k")){ stop("Expecting dimension 'k' to be passed."); }
      idx_t k = args["k"];
      auto tr = st::max_skeleton< true >(st, init, k);
      traverse(tr, f);
      break; 
    }
    case MAXIMAL: {
      auto tr = st::maximal< true >(st, init);
      traverse(tr, f);
      break; 
    }
    case LINK: {
      auto tr = st::link< true >(st, init);
      traverse(tr, f);
      break; 
    }
  }
}

// To validate the traversal parameters
param_pack validate_params(List args){
  // Extract parameters 
  auto args_str = as< vector< std::string > >(args.names());
  
  // Extract tree
  if (!contains_arg(args_str, ".ptr")){ stop("Simplex tree pointer missing."); }
  SEXP xptr = args[".ptr"]; 
  if (TYPEOF(xptr) != EXTPTRSXP || R_ExternalPtrAddr(xptr) == NULL){
    stop("Invalid pointer to simplex tree.");
  }
  XPtr< SimplexTree > st(xptr); // Unwrap XPtr
  
  // Extract initial simplex
  node_ptr init = nullptr; 
  if (!contains_arg(args_str, "sigma")){ init = st->root.get(); }
  else {
    IntegerVector sigma = args["sigma"];
    init = st->find_node(as< simplex_t >(sigma)); 
    if (init == nullptr){ init = st->root.get(); }
  }
  if (init == nullptr){ stop("Invalid starting simplex"); }
  
  // Extract traversal type 
  int tt = args["traversal_type"];
  if (tt < 0 || tt >= N_TRAVERSALS){ stop("Unknown traversal type."); }
  
  return(std::make_tuple(static_cast< SimplexTree* >(st), init, static_cast< TRAVERSAL_TYPE >(tt)));
}

// [[Rcpp::export]]
void traverse_R(List args, Function f){
  const auto run_Rf = [&f](node_ptr cn, idx_t depth, simplex_t tau){ f(tau); return true; };
  traverse_switch(validate_params(args), args, run_Rf);
}

// [[Rcpp::export]]
List ltraverse_R(List args, Function f){
  List res = List(); 
  auto run_Rf = [&f, &res](node_ptr cn, idx_t d, simplex_t tau){
    res.push_back(f(wrap(tau)));
    return(true);
  };
  traverse_switch(validate_params(args), args, run_Rf);
  return(res);
}

// [[Rcpp::export]]
SEXP straverse_R(List args, Function f) {
  Environment base("package:base");
  Function s2arr = base["simplify2array"];
  return s2arr(ltraverse_R(args, f));
}



