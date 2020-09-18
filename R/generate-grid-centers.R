#' Compute center coordinates of cells in a regular grid
#'
#' Given a rectangular extent and a resolution in both the x- and y-directions, this function computes the coordinates of the centers of the grid cells.
#'
#' @param x.lower,x.upper,y.lower,y.upper Specifies the boundaries of the rectangular
#'     region that the GRF is defined on.
#' @param resolution.x,resolution.y The number of grid cells used for the regular grid,
#'     in the x- and y-direction, respectively.
#' @param params The parameters above can instead be specified in the list \code{params}.
#' @param separate If true, the function returns the x- and y-coordinates of the centers in separate vectors.
#'
#' @return A \code{(resolution.x*resolution.y) x 2} matrix specifying the centers of the grid cells. If \code{separate = FALSE}, a list of two vectors of length resolution.x and resolution.y.
#'
#' @examples
#' # Create a 15 x 15 regular grid on [0, 1] x [0, 1]
#' reg.grid = generate.grid.centers(0, 1, 0, 1, 15, 15)
#' plot(reg.grid, asp = 1, pch = 16)
#' abline(h = seq(0, 1, length.out = 16),
#'        v = seq(0, 1, length.out = 16))
#'
#' @export
#' @author Mathias Isaksen
generate.grid.centers = function(x.lower, x.upper, y.lower, y.upper, resolution.x, resolution.y, params, separate = FALSE) {
  if (!missing(params)) {
    # Unpack parameters to function scope
    list2env(params, environment())
  }
  h.x = (x.upper-x.lower)/resolution.x
  h.y = (y.upper-y.lower)/resolution.y
  x.grid = seq(x.lower+h.x/2, x.upper-h.x/2, length.out = resolution.x)
  y.grid = seq(y.lower+h.y/2, y.upper-h.y/2, length.out = resolution.y)
  if (separate) {
    return(list(x = x.grid, y = y.grid))
  } else {
    s = expand.grid(x.grid, y.grid)
    colnames(s) = c("x", "y")
    return(s)
  }
}
