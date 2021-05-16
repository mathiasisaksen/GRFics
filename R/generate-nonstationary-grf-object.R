#' Create object used for approximate computations with a non-stationary GRF
#'
#' Creates a list containing the elements necessary for doing essential computations with the non-stationary GRF approximation.
#'
#' @note For an explanation of how `range.function`, `scale.function` and `vector.field` affect the GRF, see the corresponding parameters `range.parameter`, `scale.parameter`, `strength.parameter` and `direction.parameter` in \code{\link{generate.grf.object}}.
#' While the parameter interpretations for the stationary GRF only hold as an approximation in the non-stationary case, it still provides useful intuition.
#'
#' @param x.lower,x.upper,y.lower,y.upper Specifies the boundaries of the rectangular
#'     region that the GRF is defined on.
#' @param resolution.x,resolution.y The number of grid cells used for the regular grid,
#'     in the x- and y-direction, respectively.
#' @param range.function A strictly positive function that specifies the range of the GRF in each location. The function must be vectorized: It takes in an n x 2 matrix of locations and returns a vector of length n, containing the value of the range in each location.
#' @param scale.function A strictly positive function that is used for controlling the scale (standard deviation) of the GRF in each location. The function must be vectorized: It takes in an n x 2 matrix of locations and returns a vector of length n, containing the value of the scale in each location.
#' @param vector.field A vector field that controls nature of the non-stationarity in each location.
#' For a given location, the length and the direction of the corresponding vector plays the same role as `strength.parameter` and `direction.parameter` in \code{\link{generate.grf.object}}.
#' More precisely: Let X be a location. If `r = range.function(X)`, and `s` and `theta` are the length and direction of `vector.field(X)`, then the
#' longest range is going to be approximately `r*(1+s)`, along the direction of `theta`. The shortest range is approximately `r`, along the direction perpendicular to `theta`.
#' @param params Optional, a list containing all of the above parameters.
#' @param initial.seed The initial seed value used when generating functions.
#' @param num.functions The number of functions to initialize the GRF object with.
#' @param function.names The names of the functions to initialize the GRF object with.
#'
#' @details If the GRF object is to be initialized with functions, either \code{num.functions} or \code{function.names} must be specified.
#'
#' @return An object (technically, a list) containing:
#' \item{model.components}{A list containing the components defining the GRF: The precision matrix \code{Q}, and the Cholesky decomposition of \code{Q}, \code{chol.object}. The latter is computed and added when a function is added to the GRF object.}
#' \item{grid}{A \code{(resolution.x*resolution.y) x 2} matrix containing the center coordinates of the regular grid that the approximation is defined on.}
#' \item{initial.seed}{The initial seed value used when generating functions.}
#' \item{params}{A list of the remaining input parameters listed above.}
#'
#' @examples
#' # Example of non-stationary GRF on the square [-1, 1]^2
#' library(GRFics)
#'
#' # A range function that increases linearly with distance.
#' # The range is 0.01 in the origin (0, 0), and 0.5 in the corners (±1, ±1).
#' range.function = function(X) {
#'   r = sqrt(X[, 1]^2+X[, 2]^2)
#'   return(0.01 + (0.5-0.01)*r/sqrt(2))
#' }
#'
#' # A spiral-like vector field, where the length increases with distance from origin.
#' vector.field = function(X) {
#'   r = sqrt(X[, 1]^2+X[, 2]^2)
#'   theta = atan2(X[, 2], X[, 1]) + pi/4
#'   return(data.frame(vx = 3*r*cos(theta), vy = 3*r*sin(theta)))
#' }
#'
#' library(ggplot2)
#' # Plots of range.function and vector.field
#' plot.grid = generate.grid.centers(-1, 1, -1, 1, 20, 20)
#'
#' range.df = cbind(plot.grid, range = range.function((plot.grid)))
#' # Length of vector.field is scaled by 0.02 for plotting
#' vector.df = cbind(plot.grid, 0.02*vector.field(plot.grid))
#' vector.length = sqrt(vector.df$vx^2+vector.df$vy^2)
#' # range.function is shown using color, while vector.field is shown using arrows
#' ggplot()+
#'   geom_raster(data = range.df, aes(x = x, y = y, fill = range))+
#'   geom_segment(data = vector.df, aes(x = x-vx/2, xend = x+vx/2,
#'                                      y = y-vy/2, yend = y+vy/2),
#'                arrow = arrow(length = unit(0.2*vector.length, "npc")))+
#'   scale_fill_gradient(low = "blue", high = "white")+
#'   coord_fixed()
#'
#' # Prepare GRF object
#' grf.object = generate.nonstationary.grf.object(-1, 1, -1, 1, 100, 100,
#'                                  range.function = range.function,
#'                                  vector.field = vector.field,
#'                                  function.names = "Original",
#'                                  initial.seed = 1000)
#'
#' # Extract data frame containing the 100 x 100 grid coordinates and
#' # the value of the GRF in each grid cell.
#' original.df = get.function.df(grf.object)
#' # Interpolate GRF on finer grid for comparison
#' interp.locations = generate.grid.centers(-1, 1, -1, 1, 500, 500)
#' interp.values = evaluate.grf(interp.locations, grf.object, rescale.method = "none")
#' interp.df = data.frame(x = interp.locations[, 1],
#'                        y = interp.locations[, 2],
#'                        z = interp.values,
#'                        name = "Interpolation")
#'
#' # Plot of original and interpolated function
#' ggplot()+
#'   geom_raster(data = rbind(original.df, interp.df),
#'               aes(x = x, y = y, fill = z))+
#'   coord_fixed()+
#'   facet_grid(cols = vars(name))+
#'   scale_fill_gradient(low = "black", high = "white")
#'
#' # Create new vector field, where the direction is given by the non-stationary GRF
#' vector.grid = grf.object$grid # Use same 100 x 100 grid as the GRF is defined on
#' # Use uniform.transform = T to ensure that values are between 0 and 1 and all directions are equally respresented,
#' # and multiply by 2pi to get directions between 0 and 2pi (in radians)
#' vector.direction = 2*pi*evaluate.grf(vector.grid, grf.object, rescale.method = "uniform")
#' nonstat.vector.df = data.frame(
#'   x = vector.grid$x,
#'   y = vector.grid$y,
#'   vx = 0.02*cos(vector.direction),
#'   vy = 0.02*sin(vector.direction))
#' # Plot of random vector field (should be viewed in a large window)
#' ggplot()+
#'   geom_segment(data = nonstat.vector.df, aes(x = x-vx/2, xend = x+vx/2,
#'                                      y = y-vy/2, yend = y+vy/2),
#'                arrow = arrow(length = unit(0.003, "npc")))+
#'   coord_fixed()
#'
#' @export
#' @importFrom Rdpack reprompt
#' @author Mathias Isaksen \email{mathiasleanderi@@gmail.com}
generate.nonstationary.grf.object = function(x.lower = -1, x.upper = 1, y.lower = -1, y.upper = 1,
                               resolution.x = 100, resolution.y = 100,
                               range.function = 1, scale.function = 1,
                               vector.field = c(0, 0),
                               initial.seed = 0,
                               params = NULL,
                               num.functions = NULL,
                               function.names = NULL) {
  required.parameters = c("x.lower", "x.upper", "y.lower", "y.upper",
                          "resolution.x", "resolution.y", "range.function", "scale.function",
                          "vector.field")
  if(!missing(params)) {
    cat("Using parameters given in params.")
    missing.parameters = setdiff(required.parameters, names(params))
    if (length(missing.parameters) > 0) {
      stop(sprintf("The following values are missing in params: %s", paste(missing.parameters, collapse = ", ")))
    }
  } else {
    passed = names(as.list(match.call())[-1])
    missing.parameters = setdiff(required.parameters, passed)
    if (length(missing.parameters) > 0) {
      cat(sprintf("Note: Default values are used for the following parameters: %s", paste(missing.parameters, collapse = ", ")))
    }
    params = list(x.lower = x.lower, x.upper = x.upper, y.lower = y.lower, y.upper = y.upper,
                  resolution.x = resolution.x, resolution.y = resolution.y,
                  range.function = range.function, scale.function = scale.function,
                  vector.field = vector.field)
  }

  Q = construct.Q.nonstationary(params = params)
  grid = generate.grid.centers(params = params)
  grf.object = list(model.components = list(Q = Q), grid = grid, initial.seed = initial.seed, params = params)
  if (!(is.null(num.functions) & is.null(function.names))) {
    grf.object = add.multiple.grf.functions(grf.object, num.functions = num.functions, function.names = function.names)
  }
  return(grf.object)
}
