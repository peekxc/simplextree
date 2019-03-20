# simplextree
This package provides an interface for simplifying basic computational operations on general simplicial complexes. This is accomplished by providing a _Simplex Tree_ data structure, implemented as a [Rcpp module](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-modules.pdf). The Simplex Tree was introduced in the following paper: 

> Boissonnat, Jean-Daniel, and Clément Maria. "The simplex tree: An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

## Installation 
The current development version can be installed with the [devtools](https://github.com/r-lib/devtools) package: 
```R
require("devtools")
devtools::install_gitub("peekxc/simplextree")
```

A stable CRAN release is planned for the future. 

## Getting Started

```R
library(simplextree)
st <- simplex_tree()
st$insert_simplex(c(1, 2, 3))
st$print_tree()
# 1 (h = 2): .( 2 3 )..( 3 )
# 2 (h = 1): .( 3 )
# 3 (h = 0): 
```

## API 
```R
## see ?simplextree
```
...TODO 
