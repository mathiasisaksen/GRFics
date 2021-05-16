#' Get a data frame containing grid cell centers and corresponding values of multiple functions.
#'
#' This function returns a data frame containing 4 columns: The x- and y-coordinates of the centers of the grid cells in the regular grid, the value of the GRF in each location, and the name of the function.
#'
#' @param grf.object The GRF object of interest.
#' @param function.numbers The indices of the functions of interest.
#' @param function.names The names of the functions of interest.
#' @param all.functions If true, the data frame will include all of the functions.
#'
#' @details If more than one of \code{function.names}, \code{all.functions} and \code{function.numbers} is specified, the prioritization is \code{function.names} > \code{all.functions} > \code{function.numbers}.
#'
#' @return A \code{(num.functions*resolution.x*resolution.y) x 4} \code{data.frame}.
#'
#' @examples
#' library(GRFics)
#' # Prepare GRF object
#' grf.object = generate.grf.object(0, 1, 0, 1, 25, 25,
#'                                  strength.parameter = 2,
#'                                  direction.parameter = pi/4,
#'                                  initial.seed = 1000,
#'                                  num.functions = 1)
#' # Extract data frame containing the 50 x 50 grid coordinates and
#' # the value of the GRF in each grid cell.
#' original.df = get.function.df(grf.object)
#' # Print first portion of original.df
#' head(original.df)
#'
#' @export
#' @author Mathias Isaksen \email{mathiasleanderi@@gmail.com}

get.function.df = function(grf.object, function.numbers = 1, function.names = NULL, all.functions = FALSE) {
  if (!is.null(function.names)) {
    chosen.functions = function.names
    chosen.names = function.names
  } else if (all.functions) {
    chosen.functions = 1:length(grf.object$functions)
    chosen.names = names(grf.object$functions)
  } else {
    chosen.functions = function.numbers
    chosen.names = names(grf.object$functions)[function.numbers]
  }

  z.total = NULL
  name.total = NULL
  grid.size = nrow(grf.object$grid)

  for (i in 1:length(chosen.functions)) {
    z.total = c(z.total, grf.object$functions[[chosen.functions[i]]])
    name.total = c(name.total, rep(chosen.names[i], grid.size))
  }

  return(data.frame(x = grf.object$grid[, 1],
                    y = grf.object$grid[, 2],
                    z = z.total,
                    name = name.total))
}
