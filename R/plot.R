# ---- plot_Rcpp_SimplexTree ----
#' @name plot_simplextree
#' @title Plots the simplex tree
#' @param x a simplex tree.
#' @param coords Optional (n x 2) matrix of coordinates, where n is the number of 0-simplices. 
#' @param vertex_opt Optional parameters to modify default vertex plotting options. Passed to \code{\link[graphics]{points}}.
#' @param text_opt Optional parameters to modify default vertex text plotting options. Passed to \code{\link[graphics]{text}}.
#' @param edge_opt Optional parameters to modify default edge plotting options. Passed to \code{\link[graphics]{segments}}.
#' @param polygon_opt Optional parameters to modify default k-simplex plotting options for k > 1. Passed to \code{\link[graphics]{polygon}}.
#' @param color_pal Optional vector of colors. See details.
#' @param maximal Whether to draw only the maximal faces of the complex. Defaults to true. 
#' @param by_dim Whether to apply (and recycle or truncate) the color palette to the dimensions rather than to the individual simplices. Defaults to true.
#' @param add Whether to add to the plot or redraw. Defaults to false. See details.
## @param clip_polygons Whether to clip the polygons. Useful when visualizing large complexes. See details. 
#' @param ... unused
#' @details This function allows generic plotting of simplicial complexes using base \code{\link[graphics:graphics-package]{graphics}}.\cr
#' If not (x,y) coordinates are supplied via \code{coords}, a default layout is generated via phyllotaxis arrangement. This layout is 
#' not in general does not optimize the embedding towards any usual visualization criteria e.g. it doesn't try to separate connected components, 
#' minimize the number of crossings, etc. For those, the user is recommended to look in existing code graph drawing libraries, e.g. igraphs 'layout.auto' function, etc. 
#' The primary benefit of the default phyllotaxis arrangement is that it is deterministic and fast to generate. 
#' \cr
#' All parameters passed via list to \code{vertex_opt}, \code{text_opt}, \code{edge_opt}, \code{polygon_opt} 
#' override default parameters and are passed to \code{\link[graphics]{points}}, \code{\link[graphics]{text}}, \code{\link[graphics]{segments}}, 
#' and \code{\link[graphics]{polygon}}, respectively.\cr
#' \cr
#' If \code{add} is true, the plot is not redrawn. \cr
#' \cr
#' If \code{maximal} is true, only the maximal simplices are drawn. \cr
#' \cr
#' The \code{color_pal} argument controls how the simplicial complex is colored. It can be specified in multiple ways.
#' \enumerate{
#'   \item A vector of colors of length \emph{dim+1}, where \emph{dim}=\code{x$dimension}
#'   \item A vector of colors of length \emph{n}, where \emph{n}=\code{sum(x$n_simplices)}
#'   \item A named list of colors
#' }
#' Option (1) assigns every simplex a color based on its dimension. \cr
#' \cr
#' Option (2) assigns each individual simplex a color. The vector must be specified in level-order 
#' (see \code{\link{ltraverse}} or examples below). \cr
#' \cr
#' Option (3) allows specifying individual simplices to draw. It expects a named list, where the names
#' must correspond to simplices in \code{x} as comma-separated strings and whose values are colors. If 
#' option (3) is specified, this method will \emph{only} draw the simplices given in \code{color_pal}.\cr
#' \cr
#' If \code{length(color_pal)} does not match the dimension or the number of simplices in the complex, 
#' the color palette is recyled and simplices are as such. 
#' @importFrom utils modifyList
#' @examples 
#'
#' ## Simple 3-simplex 
#' st <- simplex_tree() %>% insert(list(1:4))
#' 
#' ## Default is categorical colors w/ diminishing opacity
#' plot(st)
#' 
#' ## If supplied colors have alpha defined, use that 
#' vpal <- rainbow(st$dimension + 1)
#' plot(st, color_pal = vpal)
#' 
#' ## If alpha not supplied, decreasing opacity applied
#' plot(st, color_pal = substring(vpal, first=1, last=7))
#' 
#' ## Bigger example; observe only maximal faces (+vertices and edges) are drawn
#' st <- simplex_tree(list(1:3, 2:5, 5:9, 7:8, 10))
#' plot(st, color_pal = rainbow(st$dimension + 1))
#' 
#' ## If maximal == FALSE, every simplex is drawn (even on top of each other)
#' vpal <- rainbow(st$dimension + 1)[c(1,2,5,4,3)]
#' pal_alpha <- c(1, 1, 0.2, 0.35, 0.35)
#' vpal <- sapply(seq_along(vpal), function(i) adjustcolor(vpal[i], alpha.f = pal_alpha[i]))
#' plot(st, color_pal = vpal, maximal = FALSE)
#' 
#' ## You can also color each simplex individually by supplying a vector 
#' ## of the same length as the number of simplices. 
#' plot(st, color_pal = sample(rainbow(sum(st$n_simplices))))
#' 
#' ## The order is assumed to follow the level order traversal (first 0-simplices, 1-, etc.)
#' ## This example colors simplices on a rainbow gradient based on the sum of their labels
#' si_sum <- straverse(st %>% level_order, sum) 
#' rbw_pal <- rev(rainbow(50, start=0,end=4/6))
#' plot(st, color_pal=rbw_pal[cut(si_sum, breaks=50, labels = FALSE)])
#' 
#' ## This also makes highlighting simplicial operations fairly trivial 
#' four_cofaces <- as.list(cofaces(st, 4))
#' coface_pal <- straverse(level_order(st), function(simplex){ 
#'     ifelse(list(simplex) %in% four_cofaces, "orange", "blue") 
#' })
#' plot(st, color_pal=unlist(coface_pal))
#' 
#' ## You can also give a named list to draw individual simplices. 
#' ## **Only the maximal simplices in the list are drawn** 
#' blue_vertices <- structure(as.list(rep("blue", 5)), names=as.character(seq(5, 9)))
#' plot(st, color_pal=append(blue_vertices, list("5,6,7,8,9"="red")))
#' @export
plot.Rcpp_SimplexTree <- function(x, coords = NULL, vertex_opt=NULL, text_opt=NULL, edge_opt=NULL, polygon_opt=NULL, color_pal=NULL, maximal=TRUE, by_dim=TRUE, add=FALSE,...) { # nocov start
  stopifnot(class(x) %in% .st_classes)
  if (sum(x$n_simplices) == 0){ graphics::plot.new(); return() } 

  
  ## Default color palette; categorical diverging if (# colors) <= 9, o/w rainbow
  if (missing(color_pal) || is.null(color_pal)){  
    n_colors <- if (by_dim) x$dimension+1 else sum(x$n_simplices)
    if (n_colors <= 9){
      color_pal <- .default_st_colors[seq(n_colors)]
    } else {
      color_pal <- substr(rev(grDevices::rainbow(n_colors, start=0, end=4/6)), start=1,stop=7)
    }
  }
  
  ## Regardless of type of palette given, the result is parsed into a vector of hexadecimal colors 
  ## or length (# simplices) in breadth-first order, not including the empty face. 
  simplex_colors <- NULL # placeholder
  draw_simplex <- rep(TRUE, sum(x$n_simplices))
  dim_idx <- straverse(level_order(x), length)-1L
  is_char_vec <- all(is.character(color_pal))
  is_in <- function(lst){ 
    return(function(element) { any(sapply(lst, function(x) (all.equal(x, element) == TRUE))) })
  }
  
  ## If the maximal faces are requested, set non-maximal `draw_simplex` indices to FALSE 
  if (maximal){
    all_simplices <- as.list(level_order(x))
    max_idx <- match(as.list(maximal(x)), all_simplices)
    draw_simplex <- vector("logical", length=sum(x$n_simplices))
    draw_simplex[max_idx] <- TRUE
    draw_simplex[dim_idx %in% c(0L, 1L)] <- TRUE ## always draw points and edges
  }
  
  ## Converts non-hex colors to hex. Additionally, any 7-length hex has 
  ## simplex dimension opacity scaling applied to it
  col_to_hex <- function(cp){
    is_hex <- (substring(cp,first=1,last=1) == "#")
    is_rgb <- (nchar(cp) == 7)
    is_col <- (!is_hex | (is_hex & is_rgb))
    cp[is_col] <- apply(grDevices::col2rgb(cp[is_col]), 2, function(col){ do.call(grDevices::rgb, as.list(col/255)) })
    cp[is_col] <- alpha4sc(cp)[is_col]
    return(cp)
  }
  
  ## Case 1: color_pal is a named list where each name is a comma-separated simplex 
  if (is.list(color_pal)){
    stopifnot(is.character(names(color_pal)))
    
    ## Extract simplices in names. Check named labels are ordered + simplices exist.
    simplices <- lapply(lapply(strsplit(names(color_pal), ","), as.integer), sort)
    names(color_pal) <- sapply(simplices, function(simplex){ paste0(simplex, collapse=",") })
    stopifnot(all(x$find(simplices)))
    
    ## Color named simplex w/ color if given, otherwise use default
    si_in <- is_in(simplices)
    si_color <- function(simplex){ ifelse(!is.null(simplex) && si_in(simplex), color_pal[[paste0(simplex, collapse=",")]], NA) }
    draw_simplex <- straverse(level_order(x), si_in)
    simplex_colors <- straverse(level_order(x), si_color)
  } else if (is_char_vec && (length(color_pal) == sum(x$n_simplices))){
    ## Case 2: color_pal is character vector w/ length == # simplices
    simplex_colors <- col_to_hex(color_pal)
  } else if (is_char_vec && ((length(color_pal) == x$dimension+1L) || by_dim)){
    ## Case 3: color_pal is character vector, recycled to dimensions
    color_pal <- rep(color_pal, length.out = x$dimension+1L)
    simplex_colors <- col_to_hex(color_pal)[dim_idx+1L]
  } else if (is_char_vec){
    simplex_colors <- rep(color_pal, length.out=sum(x$n_simplices))
  } else {
    stop("Invalid color palette given. Must be either a character vector or named list. See `?plot.simplextree`.")
  }
  
  ## Get coordinates of vertices
  if (!missing(coords)){ stopifnot(is.matrix(coords) && all(dim(coords) == c(x$n_simplices[1], 2))) }
  else {
    ## If no coordinates specified, use phyllotaxis arrangement as default 
    t <- seq(x$n_simplices[1])*pi * (3 - sqrt(5))
    coords <- cbind(sin(t)*t, cos(t)*t)
  }
  
  ## Create a new plot by default unless specified otherwise 
  if (!add){
    params <- list(...)
    rel_params <- intersect(c("xlim", "ylim", "log", "asp", "xaxs", "yaxs", "lab"), names(params))
    graphics::plot.new()
    if (length(rel_params) > 0){
      default_p <- list(xlim=range(coords[,1]), ylim=range(coords[,2]))
      do.call(graphics::plot.window, modifyList(default_p, params[rel_params]))
    } else {
      graphics::plot.window(xlim=range(coords[,1]), ylim=range(coords[,2])) 
    }
  }
  
  # plot polygons for simplices of dimension 2+; omits edges and vertices
  # this just plots the triangles
  v <- x$vertices # cache vertices
  if (x$dimension >= 2L){
    for (d in seq(x$dimension, 2)){
      if (any(draw_simplex[dim_idx == d])){
        
        ## TODO: this code union's the polygons in the plane using the polyclip package,
        ## dramatically reducing the amount of rendering that needs to occur. However, the processing 
        ## time to perform the union is sometimes larger on my machine than the naive method of rendering 
        ## many polygons, making this code ineffective. Nonetheless, the logic is retained here in case this can 
        ## be improved upon later. 
        # safe_to_clip <- is_char_vec && ((length(color_pal) == x$dimension+1L) || by_dim)
        # if (clip_polygons && safe_to_clip){
        #   polys <- ltraverse(k_simplices(x, k=d), function(simplex){ 
        #     poly <- coords[match(simplex, v),] 
        #     poly <- poly[grDevices::chull(poly),,drop=FALSE]
        #     list(x=poly[,1], y=poly[,2])
        #   })
        #   subset <- (draw_simplex & (dim_idx==d))
        #   polys <- polys[subset[dim_idx==d]]
        #   clipped_polys <- Reduce(f = function(A, B){ polyclip::polyclip(A, B, op = "union") }, 
        #                           x = polys[-1], init = polys[[1]])
        #   clipped_polys <- do.call(rbind, lapply(clipped_polys, function(p){ rbind(cbind(p$x, p$y), c(NA,NA))}))
        #   params <- list(x=clipped_polys, border=NA, col=simplex_colors[head(which(dim_idx==d),1)])
        #   do.call(graphics::polygon, modifyList(params, as.list(polygon_opt)))
        # } else {
          polys <- ltraverse(k_simplices(x, k=d), function(simplex){ 
            poly <- coords[match(simplex, v),] 
            rbind(poly[grDevices::chull(poly),,drop=FALSE], c(NA, NA))
          })
          subset <- (draw_simplex & (dim_idx==d))
          polys_to_draw <- polys[subset[dim_idx==d]]
          params <- list(x=do.call(rbind, polys_to_draw), border=NA, col=simplex_colors[dim_idx==d])
          do.call(graphics::polygon, modifyList(params, as.list(polygon_opt)))
        # }
      }
    }
  }
  # plot segments for edges
  if (length(x$n_simplices) >= 2 && any(draw_simplex[dim_idx == 1L])){
    lc <- apply(x$edges, 1, function(e){ t(coords[match(e, x$vertices),,drop=FALSE]) })
    subset <- (draw_simplex & (dim_idx==1L))
    e_subset <- subset[dim_idx==1L]
    params <- list(x0=lc[1,e_subset], y0=lc[2,e_subset], x1=lc[3,e_subset], y1=lc[4,e_subset], lwd=2, col=simplex_colors[subset])
    do.call(graphics::segments, modifyList(params, as.list(edge_opt)))
  }
  # plot vertices
  if (length(x$n_simplices) >= 1 && any(draw_simplex[dim_idx == 0L])){
    subset <- (draw_simplex & (dim_idx==0L))
    v_subset <- subset[dim_idx==0L]
    do.call(graphics::points, modifyList(list(x=coords[v_subset,,drop=FALSE], pch=21, bg=simplex_colors[subset], col=simplex_colors[subset], cex=2), as.list(vertex_opt)))
    do.call(graphics::text, modifyList(list(x=coords[v_subset,,drop=FALSE], labels=as.character(x$vertices)[v_subset], col="white", cex=0.75), as.list(text_opt))) 
  }
} # nocov end


#' @describeIn plot_simplextree family of plotting methods. 
#' @export
plot.Rcpp_Filtration <- function(...){
  plot.Rcpp_SimplexTree(...)
}

.default_st_colors <- c("#e41a1c", "#377eb8", "#ffff33", "#984ea3", "#ff7f00", "#4daf4a", "#a65628", "#f781bf", "#999999")

# Adjusts the simplex alphas for each dimension; expects hexadecimal
alpha4sc <- function(col_pal) {
  nc <- length(col_pal)
  if (nc == 0) { return(col_pal) }
  ext <- if (nc > 2){ seq(0.80, 0.45, length.out = nc-2) } else { NULL }
  si_alpha <- c(1, 1, ext) 
  sapply(seq_along(col_pal), function(i){ grDevices::adjustcolor(col_pal[i], alpha.f = si_alpha[i]) })
}
