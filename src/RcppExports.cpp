// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/simplextree.h"
#include <Rcpp.h>

using namespace Rcpp;

// bench_preorder
IntegerVector bench_preorder(SEXP stree);
RcppExport SEXP _simplextree_bench_preorder(SEXP streeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stree(streeSEXP);
    rcpp_result_gen = Rcpp::wrap(bench_preorder(stree));
    return rcpp_result_gen;
END_RCPP
}
// bench_levelorder
IntegerVector bench_levelorder(SEXP stree);
RcppExport SEXP _simplextree_bench_levelorder(SEXP streeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stree(streeSEXP);
    rcpp_result_gen = Rcpp::wrap(bench_levelorder(stree));
    return rcpp_result_gen;
END_RCPP
}
// cech_complex
void cech_complex(SEXP stx, const size_t k, Function seb_f);
RcppExport SEXP _simplextree_cech_complex(SEXP stxSEXP, SEXP kSEXP, SEXP seb_fSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stx(stxSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    Rcpp::traits::input_parameter< Function >::type seb_f(seb_fSEXP);
    cech_complex(stx, k, seb_f);
    return R_NilValue;
END_RCPP
}
// n_choose_k
size_t n_choose_k(const size_t n, const size_t k);
RcppExport SEXP _simplextree_n_choose_k(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(n_choose_k(n, k));
    return rcpp_result_gen;
END_RCPP
}
// inv_choose_2_R
size_t inv_choose_2_R(const size_t x);
RcppExport SEXP _simplextree_inv_choose_2_R(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const size_t >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_choose_2_R(x));
    return rcpp_result_gen;
END_RCPP
}
// to_subscript_R
IntegerMatrix to_subscript_R(IntegerVector numbers, const size_t n, const size_t k);
RcppExport SEXP _simplextree_to_subscript_R(SEXP numbersSEXP, SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type numbers(numbersSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(to_subscript_R(numbers, n, k));
    return rcpp_result_gen;
END_RCPP
}
// to_natural_R
IntegerVector to_natural_R(IntegerMatrix m, const size_t n);
RcppExport SEXP _simplextree_to_natural_R(SEXP mSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(to_natural_R(m, n));
    return rcpp_result_gen;
END_RCPP
}
// nfold_intersection
bool nfold_intersection(vector< vector< int > > x, const size_t n);
RcppExport SEXP _simplextree_nfold_intersection(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vector< vector< int > > >::type x(xSEXP);
    Rcpp::traits::input_parameter< const size_t >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(nfold_intersection(x, n));
    return rcpp_result_gen;
END_RCPP
}
// nerve_expand
void nerve_expand(SEXP stx, vector< size_t > ids, vector< vector< int > > cover, const size_t k, const size_t threshold);
RcppExport SEXP _simplextree_nerve_expand(SEXP stxSEXP, SEXP idsSEXP, SEXP coverSEXP, SEXP kSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stx(stxSEXP);
    Rcpp::traits::input_parameter< vector< size_t > >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< vector< vector< int > > >::type cover(coverSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    Rcpp::traits::input_parameter< const size_t >::type threshold(thresholdSEXP);
    nerve_expand(stx, ids, cover, k, threshold);
    return R_NilValue;
END_RCPP
}
// nerve_expand_f
void nerve_expand_f(SEXP stx, vector< size_t > ids, Function include_f, const size_t k);
RcppExport SEXP _simplextree_nerve_expand_f(SEXP stxSEXP, SEXP idsSEXP, SEXP include_fSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stx(stxSEXP);
    Rcpp::traits::input_parameter< vector< size_t > >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< Function >::type include_f(include_fSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    nerve_expand_f(stx, ids, include_f, k);
    return R_NilValue;
END_RCPP
}
// expand_f_bernoulli
void expand_f_bernoulli(SEXP stx, const size_t k, const double p);
RcppExport SEXP _simplextree_expand_f_bernoulli(SEXP stxSEXP, SEXP kSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type stx(stxSEXP);
    Rcpp::traits::input_parameter< const size_t >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    expand_f_bernoulli(stx, k, p);
    return R_NilValue;
END_RCPP
}
// profile
NumericVector profile(SEXP st);
RcppExport SEXP _simplextree_profile(SEXP stSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type st(stSEXP);
    rcpp_result_gen = Rcpp::wrap(profile(st));
    return rcpp_result_gen;
END_RCPP
}
// parameterize_R
List parameterize_R(SEXP st, IntegerVector sigma, std::string type, Rcpp::Nullable<List> args);
RcppExport SEXP _simplextree_parameterize_R(SEXP stSEXP, SEXP sigmaSEXP, SEXP typeSEXP, SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type st(stSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<List> >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(parameterize_R(st, sigma, type, args));
    return rcpp_result_gen;
END_RCPP
}
// traverse_R
void traverse_R(List args, Function f);
RcppExport SEXP _simplextree_traverse_R(SEXP argsSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    traverse_R(args, f);
    return R_NilValue;
END_RCPP
}
// ltraverse_R
List ltraverse_R(List args, Function f);
RcppExport SEXP _simplextree_ltraverse_R(SEXP argsSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(ltraverse_R(args, f));
    return rcpp_result_gen;
END_RCPP
}
// straverse_R
SEXP straverse_R(List args, Function f);
RcppExport SEXP _simplextree_straverse_R(SEXP argsSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type args(argsSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(straverse_R(args, f));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_union_find_module();
RcppExport SEXP _rcpp_module_boot_simplex_tree_module();
RcppExport SEXP _rcpp_module_boot_filtration_module();

static const R_CallMethodDef CallEntries[] = {
    {"_simplextree_bench_preorder", (DL_FUNC) &_simplextree_bench_preorder, 1},
    {"_simplextree_bench_levelorder", (DL_FUNC) &_simplextree_bench_levelorder, 1},
    {"_simplextree_cech_complex", (DL_FUNC) &_simplextree_cech_complex, 3},
    {"_simplextree_n_choose_k", (DL_FUNC) &_simplextree_n_choose_k, 2},
    {"_simplextree_inv_choose_2_R", (DL_FUNC) &_simplextree_inv_choose_2_R, 1},
    {"_simplextree_to_subscript_R", (DL_FUNC) &_simplextree_to_subscript_R, 3},
    {"_simplextree_to_natural_R", (DL_FUNC) &_simplextree_to_natural_R, 2},
    {"_simplextree_nfold_intersection", (DL_FUNC) &_simplextree_nfold_intersection, 2},
    {"_simplextree_nerve_expand", (DL_FUNC) &_simplextree_nerve_expand, 5},
    {"_simplextree_nerve_expand_f", (DL_FUNC) &_simplextree_nerve_expand_f, 4},
    {"_simplextree_expand_f_bernoulli", (DL_FUNC) &_simplextree_expand_f_bernoulli, 3},
    {"_simplextree_profile", (DL_FUNC) &_simplextree_profile, 1},
    {"_simplextree_parameterize_R", (DL_FUNC) &_simplextree_parameterize_R, 4},
    {"_simplextree_traverse_R", (DL_FUNC) &_simplextree_traverse_R, 2},
    {"_simplextree_ltraverse_R", (DL_FUNC) &_simplextree_ltraverse_R, 2},
    {"_simplextree_straverse_R", (DL_FUNC) &_simplextree_straverse_R, 2},
    {"_rcpp_module_boot_union_find_module", (DL_FUNC) &_rcpp_module_boot_union_find_module, 0},
    {"_rcpp_module_boot_simplex_tree_module", (DL_FUNC) &_rcpp_module_boot_simplex_tree_module, 0},
    {"_rcpp_module_boot_filtration_module", (DL_FUNC) &_rcpp_module_boot_filtration_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_simplextree(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
