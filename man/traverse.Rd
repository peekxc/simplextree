% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplex_tree.R
\name{traverse}
\alias{traverse}
\alias{straverse}
\alias{ltraverse}
\title{traverse}
\usage{
traverse(traversal, f, ...)

straverse(traversal, f, ...)

ltraverse(traversal, f, ...)
}
\arguments{
\item{traversal}{the type of traversal.}

\item{f}{the function to apply to each simplex.}

\item{...}{unused.}
}
\value{
NULL; for list or vector-valued returns, use \code{ltraverse} and \code{straverse} respectively.
}
\description{
Traverses specific subsets of a simplicial complex.
}
\details{
\code{\link{traverse}} allows for traversing ordered subsets of the simplex tree. 
The specific subset and order are determined by the choice of \emph{traversal}: examples include 
the \code{\link{preorder}} traversal, the \code{\link{cofaces}} traversal, etc. See the links below. 
Each simplex in the traversal is passed as the first and only argument to \code{f}, one per simplex in the traversal.
\code{\link{traverse}} does nothing with the result; if you want to collect the results of applying \code{f} to each simplex 
into a list, use \code{\link{ltraverse}} (or \code{\link{straverse}}), which are meant to be used like \code{\link{lapply}} 
and \code{\link{sapply}}, respectively.
}
\examples{
## Starter example complex 
st <- simplex_tree()
st \%>\% insert(list(1:3, 2:5))

## Print out complex using depth-first traversal. 
st \%>\% preorder() \%>\% traverse(print)

## Collect the last labels of each simplex in the tree. 
last_labels <- st \%>\% preorder() \%>\% straverse(function(simplex){ tail(simplex, 1) })
}
\concept{traversals}
