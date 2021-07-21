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
template < typename Lambda >
void vector_handler(SEXP sigma, Lambda&& f){
  const unsigned int s_type = TYPEOF(sigma);
  if (!Rf_isNull(Rf_getAttrib(sigma, R_DimSymbol))){
    IntegerMatrix m = as< IntegerMatrix >(sigma);
    const size_t n = m.ncol();
    for (size_t i = 0; i < n; ++i){
      if (i % 1000 == 0){ Rcpp::checkUserInterrupt(); }
      IntegerMatrix::Column cc = m(_,i);
      f(simplex_t(cc.begin(), cc.end()));
    }
  } else if (s_type == INTSXP || s_type == REALSXP){
    f(as< simplex_t >(sigma));
  } else if (s_type == LISTSXP || s_type == VECSXP){
    List simplices = List(sigma);
    const size_t n = simplices.size();
    for (size_t i = 0; i < n; ++i){
      if (i % 1000 == 0){ Rcpp::checkUserInterrupt(); }
      f(as< simplex_t >(simplices[i]));
    }
  } else { stop("Unknown type passed, must be list type or vector type."); }
}

// R-facing Inserter 
void insert_R(SimplexTree* st, SEXP x){
  SimplexTree& st_ref = *st; 
  vector_handler(x, [&st_ref](simplex_t&& sigma){ st_ref.insert(sigma); });
}


// R-facing Inserter for lexicographically-sorted column matrix
void insert_lex(SimplexTree* st, const IntegerMatrix& simplices){
  SimplexTree& st_ref = *st;
  const size_t m = simplices.ncol();
  for (size_t i = 0; i < m; ++i){
    IntegerMatrix::ConstColumn cc = simplices(_,i);
    st_ref.insert_it< true >(cc.begin(), cc.end(), st_ref.root.get(), 0);
  }
}
// void insert_lex(SimplexTree* st, const IntegerMatrix& simplices){
//   SimplexTree& st_ref = *st;
//   const size_t m = simplices.ncol();
//   const size_t d = simplices.nrow();
//   
//   splex_alloc_t a; 
//   splex_t k_simplex{ a };
//   k_simplex.resize(d);
//   
//   // Start the search 
//   IntegerMatrix::ConstColumn cc = simplices(_,0);
//   // node_ptr np = st_ref.find_it(cc.begin(), cc.end(), st_ref.root.get());
//   
//   // Get the initial simplex
//   st_ref.full_simplex_out(st_ref.root.get(), d-1, begin(cc)));
//   // std::copy(k_simplex.begin(), k_simplex.end(), k_p1_simplex.begin());
//   
//   while(std::equal(begin(cc), end(cc), begin(cc)) && i < m){
//     if (np != nullptr){
//       auto new_it = np->children.emplace_hint(np->children.end(), std::make_unique< node >(label, np));
//       auto child_np = (*new_it).get();
//       if (d > 1){ // keep track of nodes which share ids at the same depth
//         if (depth_index(d+1) >= level_map.size()){ level_map.resize(depth_index(d+1) + 1); }
// 	      auto& label_map = level_map[depth_index(d+1)][child_np->label];
//   	    label_map.push_back(child_np);
//       }
//       record_new_simplexes(d, 1);
//     }
//     cc = simplices(_,++i);
//     
//   }
//   
//       
//   for (size_t i = 1; i < m; ++i){
//     
//     
//     if (np != nullptr){
//       auto new_it = np->children.emplace_hint(np->children.end(), std::make_unique< node >(label, np));
//       auto child_np = (*new_it).get();
//       if (d > 1){ // keep track of nodes which share ids at the same depth
//         if (depth_index(d+1) >= level_map.size()){ level_map.resize(depth_index(d+1) + 1); }
// 	      auto& label_map = level_map[depth_index(d+1)][child_np->label];
//   	    label_map.push_back(child_np);
//       }
//       record_new_simplexes(d, 1);
//     }
//   }
// }

// R-facing remover 
void remove_R(SimplexTree* st, SEXP x){
  vector_handler(x, [&st](simplex_t&& sigma){ 
    st->remove(st->find(sigma)); 
  });
}

LogicalVector find_R(SimplexTree* st, SEXP simplices){
  LogicalVector v; 
  vector_handler(simplices, [&st, &v](simplex_t&& sigma){
    node_ptr np = st->find(sigma);
    v.push_back(np != st->root.get() && np != nullptr);
  });
  return(v);
}

bool collapse_R(SimplexTree* st, IntegerVector tau, IntegerVector sigma){
  return st->collapse(st->find(tau), st->find(sigma));
}

// Copies the contents of st1 to st2
void copy_trees(SEXP st1, SEXP st2){
  Rcpp::XPtr<SimplexTree> st1_ptr(st1), st2_ptr(st2);
  *st2_ptr = static_cast< const SimplexTree& >(*st1_ptr);
}

IntegerMatrix get_k_simplices(SimplexTree* st, const size_t k) {
  if (st->n_simplexes.size() <= k){ return IntegerMatrix(0, k+1); }
  IntegerMatrix res = IntegerMatrix(st->n_simplexes.at(k), k+1);
  size_t i = 0; 
  auto tr = st::k_simplices< true >(st, st->root.get(), k);
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
IntegerMatrix as_adjacency_matrix(SimplexTree* stp) {
  const SimplexTree& st = *stp; 
  const auto& vertices = st.root->children; 
  const size_t n = vertices.size();
  IntegerMatrix res = IntegerMatrix(n, n);

  // Fill in the adjacency matrix
  size_t i = 0; 
  for (auto& vi: vertices){
    for (auto& vj: vi->children){
      auto it = std::lower_bound(begin(vertices), end(vertices), vj->label, [](const node_uptr& cn, const idx_t label){
        return cn->label < label; 
      });
      const size_t j = std::distance(begin(vertices), it); 
      res.at(i, j) = res.at(j, i) = 1;
    }
    ++i;
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
// List as_list(SimplexTree* st){
//   List res = List();
//   vector< idx_t > all = vector< idx_t >();
//   idx_t d = 1;
//   auto bfs = st::level_order< true >(st);
//   traverse(bfs, [&res, &d, &all](node_ptr cn, idx_t depth, simplex_t sigma){
//     if (depth > d){
//       const size_t n = all.size() / d;
//       IntegerMatrix tmp = IntegerMatrix(d,n);
//       for (size_t i = 0; i < n; ++i){
//         IntegerVector col = IntegerVector(all.begin() + i*d, all.begin() + (i+1)*d);
//         tmp(_, i) = col;
//       }
//       res.push_back(tmp);
//       all.clear();
//       d = sigma.size();
//     }
//     all.insert(all.end(), sigma.begin(), sigma.end());
//     return true; 
//   });
//   return res;
// }



IntegerVector degree_R(SimplexTree* st, vector< idx_t > ids){
  // if (ids.size() == 0){ ids = wrap(get_vertices(st); }
  // Nullable< IntegerVector > ids_ = R_NilValue
  // IntegerVector ids = ids_.isUsable() ? IntegerVector(ids_) : (IntegerVector) wrap(get_vertices(st)); 
  // if (){ ; } else { ids = wrap(get_vertices(st)); }
  IntegerVector res(ids.size());
  std::transform(begin(ids), end(ids), begin(res), [&st](int id){
    return st->degree(static_cast< idx_t >(id));
  });
  return res;
}
IntegerVector degree_R_default(SimplexTree* st){
  return(degree_R(st, get_vertices(st)));
}

List adjacent_R(SimplexTree* st, vector< idx_t > ids = vector< idx_t >()){
  if (ids.size() == 0){ ids = get_vertices(st); }
  List res(ids.size());
  std::transform(begin(ids), end(ids), begin(res), [&st](int id){
    return st->adjacent_vertices(static_cast< idx_t >(id));
  });
  return res;
}


void print_tree(SimplexTree* st){ st->print_tree(Rcout); }
void print_cousins(SimplexTree* st){ st->print_cousins(Rcout); }

IntegerVector simplex_counts(SimplexTree* st){
  auto zero_it = std::find(st->n_simplexes.begin(), st->n_simplexes.end(), 0); 
  auto ne = std::vector< size_t >(st->n_simplexes.begin(), zero_it);
  return(wrap(ne));
}

// Exposed Rcpp Module 
RCPP_MODULE(simplex_tree_module) {
  using namespace Rcpp;
  Rcpp::class_<SimplexTree>("SimplexTree")
    .constructor()
    .method( "as_XPtr", &as_XPtr)
    .property("n_simplices", &simplex_counts, "Gets simplex counts")
    // .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    .property("dimension", &SimplexTree::dimension)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &get_vertices, "Returns the vertex labels as an integer vector.")
    .property("edges", &get_edges, "Returns the edges as an integer matrix.")
    .property("triangles", &get_triangles, "Returns the 2-simplices as an integer matrix.")
    .property("quads", &get_quads, "Returns the 3-simplices as an integer matrix.")
    .property("connected_components", &SimplexTree::connected_components)
    .method( "print_tree", &print_tree )
    .method( "print_cousins", &print_cousins )
    .method( "clear", &SimplexTree::clear)
    // .method( "degree", (IntegerVector (SimplexTree::*)())(&degree_R_default))
    .method( "degree", &degree_R)
    .method( "insert",  &insert_R)
    .method( "insert_lex", &insert_lex)
    .method( "remove",  &remove_R)
    .method( "find", &find_R)
    .method( "generate_ids", &SimplexTree::generate_ids)
    .method( "reindex", &SimplexTree::reindex)
    .method( "adjacent", &adjacent_R)
    .method( "expand", &SimplexTree::expansion )
    .method( "collapse", &collapse_R)
    .method( "vertex_collapse", (bool (SimplexTree::*)(idx_t, idx_t, idx_t))(&SimplexTree::vertex_collapse))
    .method( "contract", &SimplexTree::contract)
    .const_method( "is_tree", &SimplexTree::is_tree)
    .method( "as_adjacency_matrix", &as_adjacency_matrix)
    .method( "as_adjacency_list", &as_adjacency_list)
    .method( "as_edge_list", &as_edge_list)
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

void make_flag_filtration(Filtration* st, const NumericVector& D){
  if (st->n_simplexes.size() <= 1){ return; }
  const size_t ne = st->n_simplexes.at(1);
  const auto v = st->get_vertices();
  const size_t N = BinomialCoefficient(v.size(), 2);
  if (size_t(ne) == size_t(D.size())){
    vector< double > weights(D.begin(), D.end());
    st->flag_filtration(weights, false);
  } else if (size_t(D.size()) == size_t(N)){ // full distance vector passed in
    auto edge_iter = st::k_simplices< true >(st, st->root.get(), 1);
    vector< double > weights;
    weights.reserve(ne);
    st::traverse(edge_iter, [&weights, &D, &v](node_ptr np, idx_t depth, simplex_t sigma){
      auto v1 = sigma[0], v2 = sigma[1];
      auto idx1 = std::distance(begin(v), std::lower_bound(begin(v), end(v), v1));
      auto idx2 = std::distance(begin(v), std::lower_bound(begin(v), end(v), v2));
      auto dist_idx = to_natural_2(idx1, idx2, v.size());
      weights.push_back(D[dist_idx]);
      return true; 
    });
    st->flag_filtration(weights, false);
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
    .property("n_simplices", &simplex_counts, "Gets simplex counts")
    // .field_readonly("n_simplices", &SimplexTree::n_simplexes)
    .property("dimension", &SimplexTree::dimension)
    .property("id_policy", &SimplexTree::get_id_policy, &SimplexTree::set_id_policy)
    .property("vertices", &get_vertices, "Returns the vertex labels as an integer vector.")
    .property("edges", &get_edges, "Returns the edges as an integer matrix.")
    .property("triangles", &get_triangles, "Returns the 2-simplices as an integer matrix.")
    .property("quads", &get_quads, "Returns the 3-simplices as an integer matrix.")
    .property("connected_components", &SimplexTree::connected_components)
    .method( "print_tree", &print_tree )
    .method( "print_cousins", &print_cousins )
    .method( "clear", &SimplexTree::clear)
    .method( "generate_ids", &SimplexTree::generate_ids)
    .method( "reindex", &SimplexTree::reindex)
    .method( "adjacent", &adjacent_R)
    .method( "degree", &degree_R)
    .method( "insert",  &insert_R)
    .method( "insert_lex", &insert_lex)
    .method( "remove",  &remove_R)
    .method( "find", &find_R)
    .method( "expand", &SimplexTree::expansion )
    .method( "collapse", &collapse_R)
    .method( "vertex_collapse", (bool (SimplexTree::*)(idx_t, idx_t, idx_t))(&SimplexTree::vertex_collapse))
    .method( "contract", &SimplexTree::contract)
    .const_method( "is_tree", &SimplexTree::is_tree)
    .method( "as_adjacency_matrix", &as_adjacency_matrix)
    .method( "as_adjacency_list", &as_adjacency_list)
    .method( "as_edge_list", &as_edge_list)
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
    .method("flag_filtration", &make_flag_filtration, "Constructs a flag filtration")
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
  for (size_t i=0; i < size_t(res.size()); ++i) { res[i] = res[i] / n; }
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
  PREORDER = 0, LEVEL_ORDER = 1, FACES = 2, COFACES = 3, COFACE_ROOTS = 4, K_SKELETON = 5, 
  K_SIMPLICES = 6, MAXIMAL = 7, LINK = 8
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
  else if (type == "k_simplices" || type == "maximal-skeleton"){ param_res["traversal_type"] = int(K_SIMPLICES); }
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
    case K_SIMPLICES: {
      if (!contains_arg(args_str, "k")){ stop("Expecting dimension 'k' to be passed."); }
      idx_t k = args["k"];
      auto tr = st::k_simplices< true >(st, init, k);
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
    init = st->find(sigma); 
    if (init == nullptr){ init = st->root.get(); }
  }
  if (init == nullptr){ stop("Invalid starting simplex"); }
  
  // Extract traversal type 
  size_t tt = (size_t) args["traversal_type"];
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



