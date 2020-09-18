#' Get a data frame containing grid cell centers and corresponding values of the GRF
#'
#' This function returns a \code{(resolution.x*resolution.y) x 3} data frame, containing the coordinates of the centers of the grid cells in the regular grid and the value of the GRF in each location.
#'
#' @param grf.object The GRF object of interest.
#' @param realization.number The index of the realization of interest.
#'
#' @return A \code{(resolution.x*resolution.y) x 3} \code{data.frame}.
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
#' # Interpolate GRF on finer grid for comparison
#' interp.locations = generate.grid.centers(0, 1, 0, 1, 500, 500)
#' interp.values = evaluate.grf(interp.locations, grf.object, method = "linear")
#' interp.df = data.frame(x = interp.locations[, 1],
#'                        y = interp.locations[, 2],
#'                        z = interp.values)
#'
#' library(ggplot2)
#' ggplot()+
#'     geom_raster(data = interp.df, aes(x = x, y = y, fill = z))+
#'     coord_fixed()
#'     # Create and plot vector field, using the GRF to specify direction
#' vector.df = original.df
#'
#' @export
get.realization.df = function(grf.object, realization.number = 1) {
  return(data.frame(x = grf.object$grid[, 1],
                    y = grf.object$grid[, 2],
                    z = grf.object$realizations[[realization.number]]))
}
