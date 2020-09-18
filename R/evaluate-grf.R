#' Evaluate the value of a GRF realization in specified locations
#'
#' @param locations An N x 2 matrix containing the coordinates of the N locations of interest.
#' @param grf.object The GRF object containing the realization.
#' @param realization.number The index of the realization of interest.
#' @param periodic If true, locations outside the rectangular domain of \code{grf.object} are mapped to the corresponding point inside the domain. Otherwise, the locations are mapped to the nearest point on the boundary.
#' @param method Either \code{"nearest"} (which maps the location to the nearest grid cell), \code{"linear"} (which performs bilinear interpolation of the nearest four grid cells), or \code{"smooth"} (which is a smoothed version of bilinear interpolation, using the smoothstep as a weight function).
#' @param smooth.func If \code{method == "smooth"}, then \code{smooth.func} is an optional parameter specifying the function used for weighting the values in the nearest four grid cells. Should satisfy \code{smooth.func(0) = 0} and \code{smooth.func(1) = 1}.
#' @param rescale If true, a min-max scaling is applied to the interpolated values. This ensures that the values are between 0 and 1.
#'
#' @return A \code{numeric} of length N containing the interpolated values of the GRF in the specified locations.
#'
#' @note In its current state, \code{method = "smooth"} seems to have some artifacts visible at high resolutions. Therefore, it should be used with caution.
#'
#' @examples
#' # Prepare GRF object
#' grf.object = generate.grf.object(0, 1, 0, 1, 25, 25,
#'                                  strength.parameter = 2,
#'                                  direction.parameter = pi/4,
#'                                  initial.seed = 1000)
#' # Add realization to object
#' grf.object = add.grf.realization(grf.object)
#' # Extract data frame containing the 50 x 50 grid coordinates and
#' # the value of the GRF in each grid cell.
#' original.df = get.realization.df(grf.object)
#' original.df$facet.var = "Original"
#' # Interpolate GRF on finer grid for comparison
#' interp.locations = generate.grid.centers(0, 1, 0, 1, 500, 500)
#' interp.values = evaluate.grf(interp.locations, grf.object, method = "linear")
#' interp.df = data.frame(x = interp.locations[, 1],
#'                        y = interp.locations[, 2],
#'                        z = interp.values)
#' interp.df$facet.var = "Interpolation"
#'
#' library(ggplot2)
#' ggplot()+
#'   geom_raster(data = rbind(original.df, interp.df),
#'               aes(x = x, y = y, fill = z))+
#'   coord_fixed()+
#'   facet_grid(cols = vars(facet.var))
#'
#'
#' @export
#' @author Mathias Isaksen
evaluate.grf = function(locations, grf.object, realization.number = 1, periodic = TRUE, method = "linear", smooth.func = NULL, rescale = FALSE) {
  if(is.null(grf.object$realizations)) {
    stop("No realizations have been generated.")
  }
  # Unpack parameters to function scope
  list2env(grf.object$params, environment())
  norm.locs = normalize.locations(locations, grf.object)
  norm.grid = normalize.locations(grf.object$grid, grf.object)
  if (periodic) {
    norm.locs = map.locations.periodic(norm.locs)
  } else {
    norm.locs = map.locations.boundary(norm.locs)
  }
  if (method == "smooth") {
    if (is.null(smooth.func)){
      weight.func = function(x) ifelse(x <= 1 & x >= 0, 3*x^2-2*x^3, ifelse(x > 1, 1, 0))
    } else {
      weight.func = smooth.func
    }
  } else if (method == "linear") {
    weight.func = function(x) ifelse(x <= 1 & x >= 0, x, ifelse(x > 1, 1, 0))
  } else if (method == "nearest") {
    weight.func = function(x) ifelse(x >= 0.5, 1, 0)
  } else {
    stop("method must be either smooth, linear or nearest.")
  }
  # Get values of desired realization
  values = grf.object$realizations[[realization.number]]
  # Convert onto matrix form for simpler syntax
  value.matrix = matrix(values, nrow = resolution.y, ncol = resolution.x, byrow = TRUE)
  # Pad matrix for easier interpolation of locations close to left/bottom boundaries
  value.matrix = cbind(value.matrix, value.matrix[, 1])
  value.matrix = rbind(value.matrix, value.matrix[1, c(1:resolution.x, 1)])

  # Compute separate x- and y-components of grid
  grid.vectors = generate.grid.centers(0, 1, 0, 1, resolution.x, resolution.y, separate = TRUE)
  x.g = grid.vectors$x
  y.g = grid.vectors$y
  h.x = x.g[2] - x.g[1]
  h.y = y.g[2] - y.g[1]
  # Add another grid cell for dealing with periodicity
  x.g = c(x.g, x.g[resolution.x] + h.x)
  y.g = c(y.g, y.g[resolution.y] + h.y)

  # Move locations that are located to the left or underneath to padded area
  norm.locs[, 1] = ifelse(norm.locs[, 1] < min(x.g), max(x.g) - (min(x.g)-norm.locs[, 1]), norm.locs[, 1])
  norm.locs[, 2] = ifelse(norm.locs[, 2] < min(y.g), max(y.g) - (min(y.g)-norm.locs[, 2]), norm.locs[, 2])

  # Finds the indices of the columns that are to the left and right of each location
  left = round(norm.locs[, 1]/h.x)
  right = left + 1

  # Finds the indices of the rows that are to the left and right of each location
  bottom = round(norm.locs[, 2]/h.y)
  top = bottom + 1

  # Distance between grid cells, in this case constant
  width.x = x.g[right] - x.g[left]
  width.y = y.g[top] - y.g[bottom]

  # The input to the weight function is a number between 0 and 1, and the same holds for the output
  weight.x = weight.func((norm.locs[, 1]-x.g[left])/width.x)
  weight.y = weight.func((norm.locs[, 2]-y.g[bottom])/width.y)

  # For each location, compute a weighted mean of the values in the the four closest grid cells
  interpolated.values = (1-weight.x)*(1-weight.y)*value.matrix[cbind(bottom, left)] +
    weight.x*(1-weight.y)*value.matrix[cbind(bottom, right)] +
    (1-weight.x)*weight.y*value.matrix[cbind(top, left)] +
    weight.x*weight.y*value.matrix[cbind(top, right)]
  if (rescale) {
    interpolated.values = (interpolated.values - min(values))/(max(values) - min(values))
  }
  return(interpolated.values)
}
