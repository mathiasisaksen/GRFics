#' Create object used for approximate computations with a stationary GRF
#'
#' Creates a list containing the elements necessary for doing essential computations with the stationary GRF approximation.
#'
#' @param x.lower,x.upper,y.lower,y.upper Specifies the boundaries of the rectangular
#'     region that the GRF is defined on.
#' @param resolution.x,resolution.y The number of grid cells used for the regular grid,
#'     in the x- and y-direction, respectively.
#' @param range.parameter The range is a positive number. A small range gives functions
#'     that tend to have more variation over shorter distances, while a long range leads to functions where the value changes more gradually.
#'     The reciprocal 1/\code{range.parameter} has a similar interpretation as the frequency parameter in Perlin/simplex noise.
#'     In technical terms: The range is the distance at which the correlation between two locations is approximately 0.14.
#' @param scale.parameter The scale is a positive number that controls the scale of the
#'     functions. Setting \code{scale.parameter = a} is equivalent to setting
#'     \code{scale.parameter = 1} and then multiplying the functions by \code{a}.
#'     If you plan on using \code{\link{evaluate.grf}} with \code{rescale = TRUE} or \code{uniform.transform = TRUE}, this parameter can be ignored.
#'     In more technical terms, this is the marginal standard deviation of the GRF.
#' @param strength.parameter,direction.parameter These parameters are used to specify a range that depends on direction. If \code{range.parameter = r}, \code{strength.parameter = p} and \code{direction.parameter = theta}, then the range is \code{r*(1+p)} along the direction given by \code{theta}, and \code{r} along the direction perpendicular to \code{theta}. \code{direction.parameter} is the angle formed with the x-axis measured in radians.
#' @param params Optional, a list containing all of the above parameters.
#' @param initial.seed The initial seed value used when generating functions.
#' @param num.functions The number of functions to initialize the object with.
#' @param function.names The names of the functions to initialize the object with.
#'
#' @details If the object is to be initialized with functions, either \code{num.functions} or \code{function.names} must be specified.
#'
#' @return An object (technically, a list) containing:
#' \item{model.components}{A list containing the components defining the GRF: The precision matrix \code{Q}, and the Cholesky decomposition of \code{Q}, \code{chol.object}. The latter is computed and added when a function is added to the GRF object.}
#' \item{grid}{A \code{(resolution.x*resolution.y) x 2} matrix containing the center coordinates of the regular grid that the approximation is defined on.}
#' \item{initial.seed}{The initial seed value used when generating functions.}
#' \item{params}{A list of the remaining input parameters listed above (except \code{num.functions} and \code{function.names}).}
#'
#' @examples
#' library(GRFics)
#' # Prepare GRF object on the square [-1, 1]^2 on a 25 x 25 grid
#' grf.object = generate.grf.object(-1, 1, -1, 1, 25, 25,
#'                                  strength.parameter = 1,
#'                                  direction.parameter = pi/4,
#'                                  function.names = "Original",
#'                                  initial.seed = 1000)
#' # Extract data frame containing the 25 x 25 grid coordinates and
#' # the value of the GRF in each grid cell.
#' original.df = get.function.df(grf.object)
#' # Interpolate GRF on finer grid for comparison
#' interp.locations = generate.grid.centers(-1, 1, -1, 1, 500, 500)
#' interp.df = evaluate.grf(interp.locations, grf.object, rescale.method = "none", return.df = TRUE)
#' interp.df$name = "Interpolation"
#'
#' library(ggplot2)
#' ggplot()+
#'   geom_raster(data = rbind(original.df, interp.df),
#'               aes(x = x, y = y, fill = z))+
#'   coord_fixed()+
#'   facet_grid(cols = vars(name))
#'
#' # Create and plot vector field, using the GRF to specify direction
#' vector.df = generate.grid.centers(-1, 1, -1, 1, 25, 25)
#' # Getting an angle that is between 0 and 2*pi. Using rescale.method = "uniform" ensures that every direction is equally represented.
#' vector.df$theta = 2*pi*evaluate.grf(vector.df, grf.object, rescale.method = "uniform")
#' vector.df$vx = 0.06*cos(vector.df$theta)
#' vector.df$vy = 0.06*sin(vector.df$theta)
#' # Plot of random vector field (should be viewed in a large window)
#' ggplot()+
#'   geom_segment(data = vector.df, aes(x = x-vx/2, xend = x+vx/2,
#'                                      y = y-vy/2, yend = y+vy/2),
#'                arrow = arrow(length = unit(0.012, "npc")))+
#'   coord_fixed()
#'
#' @export
#' @importFrom Rdpack reprompt
#' @author Mathias Isaksen \email{mathiasleanderi@@gmail.com}
generate.grf.object = function(x.lower = -1, x.upper = 1, y.lower = -1, y.upper = 1,
                               resolution.x = 100, resolution.y = 100,
                               range.parameter = 1, scale.parameter = 1,
                               strength.parameter = 0, direction.parameter = 0,
                               initial.seed = 0,
                               params = NULL,
                               num.functions = NULL,
                               function.names = NULL) {
  required.parameters = c("x.lower", "x.upper", "y.lower", "y.upper",
                          "resolution.x", "resolution.y", "range.parameter", "scale.parameter",
                          "strength.parameter", "direction.parameter")
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
                  range.parameter = range.parameter, scale.parameter = scale.parameter,
                  strength.parameter = strength.parameter, direction.parameter = direction.parameter)
  }

  Q = construct.Q.stationary(params = params)
  grid = generate.grid.centers(params = params)
  grf.object = list(model.components = list(Q = Q), grid = grid, initial.seed = initial.seed, params = params)
  if (!(is.null(num.functions) & is.null(function.names))) {
    grf.object = add.multiple.grf.functions(grf.object, num.functions = num.functions, function.names = function.names)
  }
  return(grf.object)
}

