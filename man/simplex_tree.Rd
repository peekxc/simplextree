% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplex-tree.R
\docType{class}
\name{simplex_tree}
\alias{simplex_tree}
\alias{SimplexTree}
\alias{Rcpp_SimplexTree}
\title{Simplex Tree}
\usage{
simplex_tree(simplices = NULL)
}
\arguments{
\item{simplices}{optional simplices to initialize the simplex tree with. See \code{\link[=insert]{insert()}}.}
}
\value{
A queryable simplex tree, as an object (Rcpp module) of class \code{"Rcpp_SimplexTree"}.
}
\description{
Simplex tree class exposed as an Rcpp Module.
}
\details{
A simplex tree is an ordered trie-like structure specialized for storing and doing general computation
simplicial complexes. Here is figure of a simplex tree, taken from the original paper (see Boissonnat and Maria, 2014):

\if{html}{\figure{simplextree.png}{options: width="80\%" alt="Figure: simplextree.png"}}
\if{latex}{\figure{simplextree.pdf}{options: width=12cm}}

The current implementation provides a subset of the functionality described in the paper.
}
\section{Fields}{

\describe{
\item{\code{n_simplices}}{A vector, where each index \eqn{k} denotes the number \eqn{(k-1)}-simplices.}

\item{\code{dimension}}{The dimension of the simplicial complex.}
}}

\section{Properties}{

Properties are actively bound shortcuts to various methods of the simplex tree that may be thought of as fields.
Unlike fields, however, properties are not explicitly stored: they are generated on access.
\describe{
\item{\verb{$id_policy}}{ The policy used to generate new vertex ids. May be assigned \code{"compressed"} or \code{"unique"}. See \code{\link[=generate_ids]{generate_ids()}}. }
\item{\verb{$vertices}}{ The 0-simplices of the simplicial complex, as a matrix. }
\item{\verb{$edges}}{ The 1-simplices of the simplicial complex, as a matrix. }
\item{\verb{$triangles}}{ The 2-simplices of the simplicial complex, as a matrix. }
\item{\verb{$quads}}{ The 3-simplices of the simplicial complex, as a matrix. }
\item{\verb{$connected_components}}{ The connected components of the simplicial complex. }
}
}

\section{Methods}{

\describe{
\item{\verb{$as_XPtr()}}{ Creates an external pointer. }
\item{\code{$}\code{\link[=generate_ids]{generate_ids()}}}{ Generates new vertex ids according to the set policy. }
\item{\code{$}\code{\link[=insert]{insert()}}}{ Inserts a simplex into the trie. }
\item{\code{$}\code{\link[=remove]{remove()}}}{ Removes a simplex from the trie. }
\item{\code{$}\code{\link[=find]{find()}}}{ Returns whether a simplex exists in the trie. }
\item{\code{$}\code{\link[=degree]{degree()}}}{ Returns the degree of each given vertex. }
\item{\code{$}\code{\link[=adjacent]{adjacent()}}}{ Returns vertices adjacent to a given vertex. }
\item{\code{$}\code{\link[=clear]{clear()}}}{ Clears the simplex tree. }
\item{\code{$}\code{\link[=expand]{expand()}}}{ Performs an \eqn{k}-expansion. }
\item{\code{$}\code{\link[=collapse]{collapse()}}}{ Performs an elementary collapse. }
\item{\code{$}\code{\link[=contract]{contract()}}}{ Performs an edge contraction. }
\item{\code{$}\code{\link[=traverse]{traverse()}}}{ Traverses a subset of the simplex tree, applying a function to each simplex. }
\item{\code{$}\code{\link[=ltraverse]{ltraverse()}}}{ Traverses a subset of the simplex tree, applying a function to each simplex and returning the result as a list. }
\item{\code{$}\code{\link[=is_face]{is_face()}}}{ Checks for faces. }
\item{\code{$}\code{\link[=is_tree]{is_tree()}}}{ Checks if the simplicial complex is a tree. }
\item{\verb{$as_list()}}{ Converts the simplicial complex to a list. }
\item{\verb{$as_adjacency_matrix()}}{ Converts the 1-skeleton to an adjacency matrix. }
\item{\verb{$as_adjacency_list()}}{ Converts the 1-skeleton to an adjacency list. }
\item{\verb{$as_edgelist()}}{ Converts the 1-skeleton to an edgelist. }
}
}

\examples{
## Recreating simplex tree from figure. 
st <- simplex_tree()
st \%>\% insert(list(1:3, 2:5, c(6, 7, 9), 7:8, 10))
plot(st)

## Example insertion
st <- simplex_tree(list(1:3, 4:5, 6)) ## Inserts one 2-simplex, one 1-simplex, and one 0-simplex
print(st) 
# Simplex Tree with (6, 4, 1) (0, 1, 2)-simplices

## More detailed look at structure
print_simplices(st, "tree")
# 1 (h = 2): .( 2 3 )..( 3 )
# 2 (h = 1): .( 3 )
# 3 (h = 0): 
# 4 (h = 1): .( 5 )
# 5 (h = 0): 
# 6 (h = 0): 
## Print the set of simplices making up the star of the simplex '2'
print_simplices(st \%>\% cofaces(2))
# 2, 2 3, 1 2, 1 2 3

## Retrieves list of all simplices in DFS order, starting with the empty face 
dfs_list <- ltraverse(st \%>\% preorder(empty_face), identity)

## Get connected components 
print(st$connected_components)
# [1] 1 1 1 4 4 5

## Use clone() to make copies of the complex (don't use the assignment `<-`)
new_st <- st \%>\% clone()

## Other more internal methods available via `$` 
list_of_simplices <- st$as_list()
adj_matrix <- st$as_adjacency_matrix()
# ... see also as_adjacency_list(), as_edge_list(), etc 
}
\references{
Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: An
efficient data structure for general simplicial complexes." Algorithmica
70.3 (2014): 406-427.
}
\seealso{
Other simplicial complex constructors: 
\code{\link{flag}()},
\code{\link{nerve}()},
\code{\link{rips}()}
}
\author{
Matt Piekenbrock
}
\concept{simplicial complex constructors}
