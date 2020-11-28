#' Get a data frame containing grid cell centers and corresponding values of multiple realizations.
#'
#' This function returns a data frame containing 4 columns: The x- and y-coordinates of the centers of the grid cells in the regular grid, the value of the GRF in each location, and the name of the realization.
#'
#' @param grf.object The GRF object of interest.
#' @param realization.numbers The indices of the realizations of interest.
#' @param realization.names The names of the realizations of interest.
#' @param all.realizations If true, the data frame will include all of the realizations.
#'
#' @details If more than one of \code{realization.names}, \code{all.realizations} and \code{realization.numbers} is specified, the prioritization is \code{realization.names} > \code{all.realizations} > \code{realization.numbers}.

#' @return A \code{(num.realizations*resolution.x*resolution.y) x 4} \code{data.frame}.
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

get.realization.df = function(grf.object, realization.numbers = 1, realization.names = NULL, all.realizations = FALSE) {
  if (!is.null(realization.names)) {
    chosen.realizations = realization.names
    chosen.names = realization.names
  } else if (all.realizations) {
    chosen.realizations = 1:length(grf.object$realizations)
    chosen.names = names(grf.object$realizations)
  } else {
    chosen.realizations = realization.numbers
    chosen.names = names(grf.object$realizations)[realization.numbers]
  }

  z.total = NULL
  name.total = NULL
  grid.size = nrow(grf.object$grid)

  for (i in 1:length(chosen.realizations)) {
    z.total = c(z.total, grf.object$realizations[[chosen.realizations[i]]])
    name.total = c(name.total, rep(chosen.names[i], grid.size))
  }

  return(data.frame(x = grf.object$grid[, 1],
                    y = grf.object$grid[, 2],
                    z = z.total,
                    name = name.total))
}
