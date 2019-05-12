## ---- echo=FALSE, fig.width=2--------------------------------------------
fig_path <- file.path(system.file(package = "simplextree"), "help", "figures", "simplextree.png")
knitr::include_graphics(path = fig_path, auto_pdf = TRUE, dpi = 330)

## ------------------------------------------------------------------------
library("simplextree")
st <- simplex_tree()

## ------------------------------------------------------------------------
## Inserts the 2-simplex { 1, 2, 3 }. Also inserts its faces.
st$insert(c(1, 2, 3))
print(st)

## ------------------------------------------------------------------------
## Inserts two 0-simplices, { 4 } and { 5 }, and one 1-simplex, { 5, 6 }
st$insert(list(4, 5, c(5, 6)))
print(st)

## ------------------------------------------------------------------------
st$find(c(2, 3)) ## TRUE
st$remove(list(4, 5))
print(st)

## ------------------------------------------------------------------------
print(st) # also see $dimension and $n_simplices
st$print_tree()

## ------------------------------------------------------------------------
st$traverse(print, "dfs") # equivalent to st$traverse(empty_face, print, "dfs")

## ------------------------------------------------------------------------
## Prints the cofaces of { 2 }
st$traverse(2L, print, "cofaces")

## ------------------------------------------------------------------------
st$traverse(empty_face, print, "maximal-skeleton", list(k = 1L))

## ------------------------------------------------------------------------
## Get the simplices in the link of the 1-simplex { 2, 3 }
st$ltraverse(3, identity, "link")

## ------------------------------------------------------------------------
## Prints the cofaces of the 0-simplices
st$traverse(empty_face, function(simplex){ 
  cofaces <- st$ltraverse(simplex, identity, "cofaces")
  cat(sprintf("Cofaces of %s: %s\n", simplex, paste0(cofaces, collapse = ", ")))
}, "maximal-skeleton", list(k = 0))

## ------------------------------------------------------------------------
## Collapse the free pair ({ 3 }, { 6 }) -> { 6 }
st$collapse(3, 6, 6)
st$print_tree()

## ------------------------------------------------------------------------
st1 <- simplex_tree()
st2 <- simplex_tree()

## Inserts 11 simplices
st1$insert(list(1, 2, 3, c(1, 2), c(2, 3), c(1, 3), c(1, 2, 3), 4, 5, c(4, 5), 6))

## serialize exports the maximal faces needed to restore the complex 
st_serialized <- st1$serialize()
writeLines(paste0(st_serialized))

## deserializing inserts the faces
st2$deserialize(st_serialized)
all.equal(st1, st2)

