#' Create object used for generating realizations from an approximate GRF
#'
#' Creates a list containing the elements necessary for doing essential computations with the GRF approximation.
#'
#' @param x.lower,x.upper,y.lower,y.upper Specifies the boundaries of the rectangular
#'     region that the GRF is defined on.
#' @param resolution.x,resolution.y The number of grid cells used for the regular grid,
#'     in the x- and y-direction, respectively.
#' @param range.parameter The range is a positive number. A small range gives realizations
#'     that tend to have more variation over shorter distances, while a long range leads
#'     to realizations that have less variation and appear smoother.
#' @param scale.parameter The scale is a positive number that controls the scale of the
#'     realizations. Setting \code{scale.parameter = a} is equivalent to setting
#'     \code{scale.parameter = 1} and then multiplying the realizations by \code{a}. In more technical terms, this is the marginal standard deviation of the GRF.
#' @param strength.parameter,direction.parameter These parameters are used to specify a range that depends on direction. If \code{strength.parameter = p} and \code{direction.parameter = theta}, then the range is \code{p} times longer along the direction given by \code{theta}, when compared to the range along the direction perpendicular to \code{theta}. \code{direction.parameter} is the angle formed with the x-axis measured in radians.
#' @param params Optional, a list containing all of the above parameters.
#' @param initial.seed The initial seed value used when generating realizations.
#'
#' @return A list containing:
#' \item{Q}{The precision matrix of the GRF approximation.}
#' \item{grid}{A \code{(resolution.x*resolution.y) x 2} matrix containing the center coordinates of the regular grid that the approximation is defined on.}
#' \item{initial.seed}{The initial seed value used when generating realizations.}
#' \item{params}{A list of the remaining input parameters listed above.}
#'
#' @examples
#' # Prepare GRF object
#' grf.object = generate.grf.object(0, 1, 0, 1, 25, 25,
#'                                  strength.parameter = 2,
#'                                  direction.parameter = pi/4,
#'                                  initial.seed = 1000)
#' # Add realization to object
#' grf.object = add.grf.realization(grf.object)
#' # Extract data frame containing the 25 x 25 grid coordinates and
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
#' # Create and plot vector field, using the GRF to specify direction
#' vector.df = original.df
#' # Getting an angle that is between 0 and 2*pi
#' vector.df$theta = with(vector.df, 2*pi*(z - min(z))/(max(z) - min(z)))
#' vector.df$vx = 0.04*cos(vector.df$theta)
#' vector.df$vy = 0.04*sin(vector.df$theta)
#' ggplot()+
#'   geom_segment(data = vector.df, aes(x = x-vx/2, xend = x+vx/2,
#'                                      y = y-vy/2, yend = y+vy/2),
#'                arrow = arrow(length = unit(4, "pt")))+
#'   coord_fixed()
#'
#' @export
#' @importFrom Rdpack reprompt
#' @author Mathias Isaksen
generate.grf.object = function(x.lower = -1, x.upper = 1, y.lower = -1, y.upper = 1,
                               resolution.x = 100, resolution.y = 100,
                               range.parameter = 1, scale.parameter = 1,
                               strength.parameter = 1, direction.parameter = 0,
                               initial.seed = 0,
                               params = NULL) {
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

  Q = construct.Q(params = params)
  grid = generate.grid.centers(params = params)
  return(list(Q = Q, grid = grid, initial.seed = initial.seed, params = params))
}

