# simplextree
`simplextree` is an [R](https://www.r-project.org/) package aimed at simplifying computation for general [simplicial complexes](https://en.wikipedia.org/wiki/Simplicial_complex). This package facilitates this aim by providing an R-bindings to a _Simplex Tree_ data structure, implemented using [C++11](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3690.pdf) and exported as a [Rcpp module](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-modules.pdf). 

The Simplex Tree was originally introduced in the following paper: 

> Boissonnat, Jean-Daniel, and Clément Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

A _Simplex Tree_ is an ordered, [trie](https://en.wikipedia.org/wiki/Trie)-like structure. Here's a picture of a simplicial complex (left) and its corresponding _Simplex Tree_ (right):

![simplex tree picture](./man/figures/simplextree.png)
 
## Installation 
The current development version can be installed with the [devtools](https://github.com/r-lib/devtools) package: 

```R
require("devtools")
devtools::install_gitub("peekxc/simplextree")
```

A stable CRAN release is planned for the future. 

## Quickstart

```R
library(simplextree)
st <- simplex_tree() ## instantiation wrapper
st$insert(list(1:3, 4:5, 6)) ## Inserts { 1, 2, 3 }, { 4, 5 }, and { 6 }

## Summary of complex
print(st) 
# Simplex Tree with (6, 4, 1) (0, 1, 2)-simplices

## More detailed look at structure
st$print_tree()
# 1 (h = 2): .( 2 3 )..( 3 )
# 2 (h = 1): .( 3 )
# 3 (h = 0): 
# 4 (h = 1): .( 5 )
# 5 (h = 0): 
# 6 (h = 0): 

## Print the set of simplices making up the star of the simplex '2'
st$traverse(2, function(simplex){ print(simplex) }, "star")
# [1] 1 2
# [1] 1 2 3
# [1] 2
# [1] 2 3

## Retrieves list of all simplices in DFS order, starting with the empty face 
dfs_list <- st$ltraverse(empty_face, identity, "dfs")

## Get connected components 
print(st$connected_components)
# [1] 1 1 1 4 4 5

## Serialization/persistent storage options available
print(st$serialize())
# [[1]]
# [1] 1 2 3
# [[2]]
# [1] 4 5
# [[3]]
# [1] 6

## As are various export options
list_of_simplices <- st$as_list()
adj_matrix <- st$as_adjacency_matrix()
# ... see also as_adjacency_list(), as_edge_list
```

## API Reference 

Below is the API reference for the R-bindings of the _SimplexTree_ class. A _SimplexTree_ object can also be passed, manipulated, and modified via C++ in another R package as well. See [Usage with Rcpp](#Usage with Rcpp). 

### SimplexTree

[#](#simplex_tree) __simplex\_tree__()

Trivial wrapper for constructing a _SimplexTree_ instance. See below. 

[#](#SimplexTree) __SimplexTree__ (_C++ Class_)

Exposed binding for an internal _SimplexTree_ C++ class. The binding is exposed as an [Rcpp Module](). Instantiation returns an object of type _Rcpp\_SimplexTree_. Instantiate with either `simplex_tree()` or `new(SimplexTree)`.

[#](#n_simplexes) _SimplexTree_ $ **n\_simplices**

An integer vector, where each index *k* denotes the number of (*k-1*)-simplices. Read-only. 

[#](#max_depth) _SimplexTree_ $ **dimension**

The dimension of a simplicial complex *K* is the highest dimension of any of *K*'s simplices. Equivalently, this is the maximum depth of any subtree in the _SimplexTree_, where the _depth_ of a *k*-simplex is *k+1*. The root node has depth 0.  

#### Modifying the tree 

[#](#insert_simplex) _SimplexTree_ $ **insert_simplex**(\[*simplex*\])

Inserts the specified _simplex_ into the simplex tree, if it doesn't already exist. The _simplex_ is ordered prior to insertion. If the _simplex_ already exists, the tree is not modified. 

Any proper faces of _simplex_ not in the tree already are also inserted.  

Note that the _SimplexTree_ does not track orientation, e.g. the simplices _(1, 2, 3)_ and _(2, 1, 3)_ are considered identical. 

[#](#insert_simplex) _SimplexTree_ $ **insert_simplices**(\[*simplices*\]))

[Inserts](#insert_simplex) every simplex in the list _simplices_ into the simplex tree.

[#](#remove_simplex) _SimplexTree_ $ **remove_simplex**(\[*simplex*\])

Removes the specified _simplex_ from the simplex tree, if it exists. If the _simplex_ doesn't exist, the tree is not modified. 

Removing a simplex will also remove all its cofaces as well. 

[#](#contract) _SimplexTree_ $ **contract**(\[*a, b*\])

Performs and *edge contraction*, contracting vertex *b* to vertex *a*. This is equivalent to removing vertex *b* from the simplex tree and augmenting the link of vertex *a* with the link of vertex *b*. If the edge does not exist in the tree, the tree is not modified.

[#](#collapse) _SimplexTree_ $ **collapse**(...)

1. (\[_tau_\], \[_sigma_\])
2. (_u_, _v_, _w_)

Performs an _elementary collapse_. There are multiple simplifying operations that have been referred to as elementary collapses; this method provides two such operations.

(1) elementary collapse ( from [1](#simplex-tree-paper) ) 

Collapses _tau_ through its coface _sigma_ if _sigma_ is the only coface of _tau_. 

(2) vertex collapse ( from [2](#simplex-tree-paper))

Collapses a free pair (_u_, _v_) -> (_w_), mapping a pair of vertices to a single vertex. 

Technically, if you want to do an _elementary_ collapse, it's required that either _u_ = _w_, such that (_u_, _v_) -> (_u_), or _v_ = _w_, such that (_u_, _v_) -> (_v_). However, if (_u_, _v_) -> (_w_) is specified, where _u_ != _w_ and _v_ != _w_ , the collapse is decomposed into two elementary collapses: (_u_, _w_) -> (_w_) and (_v_, _w_) -> (_w_). 

[#](#as_XPtr) _SimplexTree_ $ **as_XPtr**()

Converts the simplex tree into an [XPtr](https://github.com/RcppCore/Rcpp/blob/master/inst/include/Rcpp/XPtr.h), sometimes called an _external pointer_. _XPtr_'s can be passed as an [SEXP](https://cran.r-project.org/doc/manuals/r-release/R-ints.html#SEXPs) to other C/C++ functions via R's C API or Rcpp.

This method does _not_ register a delete finalizer. 

#### Querying the tree 

[#](#print_tree) _SimplexTree_ $ **print_tree**()

Prints the simplicial complex to _standard out_. By default, this is set to R's buffered output, which is shown in the R console. The printed format is: 

> \[*vertex*\] (h = \[*subtree height*\]): \[*subtree depth*\](\[*subtree*\])

Where each line corresponds to a *vertex* and its corresponding subtree. The *subtree depth* represents the set of _sibling_ *k*-simplices at that level in tree. represented by a sequence of dots ('**.**').  

[#](#find_simplex) _SimplexTree_ $ **find_simplex**(\[*simplex*\])

Traverses the simplex tree downard starting at the root by successively using each of the ordered labels in _simplex_ as keys. Returns a logical indicating whether _simplex_ was found. _simplex_ is sorted prior to traversing. 

[#](#degree) _SimplexTree_ $ **degree**(_[vertices]_)

Returns the degree of a given vertex or set of vertices, i.e. the number of 1-simplices each vertex is a face of. Higher order simplices are not counted, for that see #cofaces.

[#](#adjacent_vertices) _SimplexTree_ $ **adjacent_vertices**()

[#](#is_face) _SimplexTree_ $ **is_face**(_[tau]_, _[sigma]_)

Returns a logical indicating whether tau is a face of sigma. 

[#](#is_face) _SimplexTree_ $ **is_tree**()

Returns a logical indicating whether the simplex tree is fully connected and acyclic. 

Note that although the simplex tree is an ordered structure, in this implementation, simplex orientation (or direction, in the graph sense) is not maintained. 

#### Traversals

The _SimplexTree_ data structure supports various types of _traversals_. A _traversal_ is a (possibly optimized) path that allows iteration through a subset of the _SimplexTree_. The traversal _type_  determines the subset and path to iterate through. _f_ is called with a simplex as its first argument for each simplex in the traversal. 

[#](#traverse) _SimplexTree_ $ **traverse** <br /> 

1. (_f_, _type_)
2. (\[_simplex_\], _f_, _type_)
3. (\[_simplex_\], _f_, _type_, _params_)

The **traverse** method has three overloads, based on the traversal _type_ and intended usage of _f_.

(1) applies _f_ to each simplex in the traversal path _type_, starting at the root of the tree. 

(2) applies _f_ to each simplex in the traversal path _type_ , starting at the specified _simplex_ in the tree. The root simplex (empty face) may be specified using the _NULL_ keyword. 

(3)  applies _f_ to each simplex in the traversal path _type_, starting at the specified _simplex_ in the tree, and supplying the necessary parameters to the traversal _type_ as a list _params_. The root simplex (empty face) may be specified using the _NULL_ keyword. 

For example, to traverse the simplicial complex in a depth-first manner:

```R
st <- simplex_tree()
st$insert_simplex(1:3)
st$traverse(message, "dfs") # equivalent to 'st$traverse(NULL, message, "dfs")'
# 1
# 12
# 123
# 13
# 2
# 23
# 3
```
Or to print the cofaces of the vertex with label '1': 

```R
st$traverse(1, message, "cofaces")
# 1
# 12
# 123
# 13
```
[#](#ltraverse) _SimplexTree_ $ **ltraverse**(...)

Performs a _traversal_, returning a list of the same length as the traversal path, with each element containing the result of _f_. The parameters *...* are the same as in [traverse](#traverse). **ltraverse** is meant to used in a similar way as lapply.

##### Traversal types

The _type_ parameter passed to the traverse family of algorithms determines the subset and corresponding path that is enumerated in simplex tree. A traversal _type_ is specified by a string, and its corresponding _params_ are specified in a list. The currently supported traversal types are as follows: 

[#](#dfs) _type_ = "**dfs**" 

Performs a depth-first traversal of the _SimplexTree_ starting at _simplex_. If _simplex_ is not supplied, the traversal starts at the root node.

[#](#bfs) _type_ = "**bfs**" 

Performs a breadth-first traversal of the _SimplexTree_ starting at _simplex_. If _simplex_ is not supplied, the traversal starts at the root node.

[#](#bfs) _type_ = "**cofaces**" or _type_ = "**star**"

Traverse all of the cofaces (the star) of _simplex_ in the _SimplexTree_. If _simplex_ is not supplied, the traversal starts at the root node.

[#](#bfs) _type_ = "**link**"

Traverse all of the link of _simplex_ in the _SimplexTree_. If _simplex_ is not supplied, the traversal starts at the root node.

[#](#bfs) _type_ = "**skeleton**"

Traverses all of simplices in the the _k-skeleton_ of the _SimplexTree_, where the dimension _k_ must be supplied via _params_. If _simplex_ is not supplied, the traversal starts at the root node.

[#](#bfs) _type_ = "**maximal-skeleton**"

Traverses all of simplices in the the _maximal k-skeleton_ of the _SimplexTree_, where the dimension _k_ must be supplied via _params_. If _simplex_ is not supplied, the traversal starts at the root node.

#### Import / Export options 

[#](#serialize) _SimplexTree_ $ **serialize**()

[#](#unserialize) _SimplexTree_ $ **deserialize**()

[#](#save) _SimplexTree_ $ **save**()

[#](#load) _SimplexTree_ $ **load**()

#### Conversions

The full simplicial complex can always be converted to a list of matrices. The 1-skeleton can also be converted to any of the standard graph representations. 

[#](#as_list) _SimplexTree_ $ **as\_list**()

[#](#as_edge_list) _SimplexTree_ $ **as\_edge_list**()

[#](#as_adjacency_matrix) _SimplexTree_ $ **as\_adjacency_list**()

[#](#as_adjacency_matrix) _SimplexTree_ $ **as\_adjacency_matrix**()

## Usage with Rcpp

There are two ways to use a _SimplexTree_ object with Rcpp. 

#### 1. Pure Rcpp

If you're developing purely in Rcpp, you can just use the _SimplexTree_ class directly. The _SimplexTree_ is header-only, and can be imported via the **Rcpp::depends** attribute (see [Rcpp Attributes](http://dirk.eddelbuettel.com/code/rcpp/Rcpp-attributes.pdf))

```C
#include "Rcpp.h"

// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"

void my_function(){
  SimplexTree st = SimplexTree();
  ...
}

```
If you're developing using a package, make sure to add `simplextree` to the `LinkingTo` list to ensure the header is properly included (e.g. `-I"<...>/simplextree/include"` should appear in the build steps).

#### 2. Moving between R and Rcpp

A _SimplexTree_ Rcpp module can be passed directly from R to any Rcpp method. To do so, export the object as an external pointer (XPtr), pass as an SEXP, then reinterpret cast the SEXP on the C++ side. 

For example, on the C++ side, one might do:
 
```C
// my_source.cpp
#include "Rcpp.h"

// [[Rcpp::depends(simplextree)]]
#include "simplextree.h"

[[Rcpp::export]]
void modify_tree(SEXP stree){
  Rcpp::XPtr<SimplexTree> stree_ptr(stree);
  stree_ptr->print_tree();
  ....
}
```
Then on the R-side, use [as\_XPtr](#as\_XPtr) method to get an [XPtr](https://cran.r-project.org/doc/manuals/R-exts.html#External-pointers-and-weak-references). 

```R
# my_source.R
stree <- simplextree::simplex_tree()
modify_tree(stree$as_XPtr())
```  
Note that the C++ class contains a superset of the functionality exported to R, however they do not necessarily have the same bindings. See the header file for a complete list. 


## References 

> <a name="simplex-tree-paper">1.</a> Boissonnat, Jean-Daniel, and Clément Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427. 

> <a name="simplicial-map-paper">2.</a> Dey, Tamal K., Fengtao Fan, and Yusu Wang. "Computing topological persistence for simplicial maps." Proceedings of the thirtieth annual symposium on Computational geometry. ACM, 2014.

## More 
```R
## see ?simplextree
```
...TODO 
