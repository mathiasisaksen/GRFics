# Note: These functions are considered internal, and are therefore not documented.

normalize.locations = function(locations, grf.object) {
  # Unpack parameters to function scope
  list2env(grf.object$params, environment())
  if ("x" %in% colnames(locations) && "y" %in% colnames((locations))) {
    x = locations[, "x"]
    y = locations[, "y"]
  } else {
    x = locations[, 1]
    y = locations[, 2]
  }
  normalized.locations = data.frame(x = (x - x.lower)/(x.upper - x.lower),
                                    y = (y - y.lower)/(y.upper - y.lower))
  return(normalized.locations)
}

map.locations.periodic = function(normalized.locations) {
  result = normalized.locations
  result[, "x"] = result[, "x"] - floor(result[, "x"])
  result[, "y"] = result[, "y"] - floor(result[, "y"])
  return(result)
}

map.locations.boundary = function(normalized.locations) {
  result = normalized.locations
  # Should be made more efficient, loops through dataframe 4 times
  result[result[, "x"] > 1, "x"] = 1
  result[result[, "x"] < 0, "x"] = 0
  result[result[, "y"] > 1, "y"] = 1
  result[result[, "y"] < 0, "y"] = 0
  return(result)
}

# Given a grid and locations, returns the index of the grid element closest to each location
match.location.to.grid = function(locations, grf.object = NULL, params = NULL, periodic = TRUE) {
  if (!is.null(params)) {
    # Unpack parameters to function scope
    list2env(params, environment())
  } else if (!is.null(grf.object)) {
    list2env(grf.object$params, environment())
  } else {
    stop("Either grf.object or params must be specified.")
  }
  normalized.locations = normalize.locations(locations, grf.object)
  if (periodic) {
    normalized.locations = map.locations.periodic(normalized.locations)
  } else {
    normalized.locations = map.locations.boundary(normalized.locations)
  }
  col.row.indices = round(cbind((normalized.locations[, 1]-1/(2*resolution.x))*resolution.x + 1,
                                (normalized.locations[, 2]-1/(2*resolution.y))*resolution.y + 1))
  col.row.indices[col.row.indices[, 1] < 1, 1] = 1
  col.row.indices[col.row.indices[, 1] > resolution.x, 1] = resolution.x
  col.row.indices[col.row.indices[, 2] < 1, 2] = 1
  col.row.indices[col.row.indices[, 2] > resolution.y, 2] = resolution.y
  grid.indices = as.integer((col.row.indices[, 2]-1) * resolution.x + col.row.indices[, 1])
  return(grid.indices)
}

match.location.to.grid.manual = function(grid, locations) {
  x.l = min(grid[, 1])
  x.u = max(grid[, 1])
  y.l = min(grid[, 2])
  y.u = max(grid[, 2])
  n.x = length(unique(grid[, 1]))
  n.y = length(unique(grid[, 2]))

  locations.normalized = cbind((locations[, 1] - x.l)/(x.u - x.l), (locations[, 2] - y.l)/(y.u - y.l))
  locations.normalized[locations.normalized[, 1] > 1, 1] = 1
  locations.normalized[locations.normalized[, 1] < 0, 1] = 0
  locations.normalized[locations.normalized[, 2] > 1, 2] = 1
  locations.normalized[locations.normalized[, 2] < 0, 2] = 0
  row.col.indices = t(apply(locations.normalized, 1, function(x) round(c(x[1] * (n.x - 1), x[2] * (n.y - 1))) + 1))
  grid.indices = as.integer((row.col.indices[, 2] - 1) * n.x + row.col.indices[, 1])
  return(grid.indices)
}

generate.grid.centers.HS = function(x.lower, x.upper, y.lower, y.upper, resolution.x, resolution.y, params, separate = FALSE) {
  if (!missing(params)) {
    # Unpack parameters to function scope
    list2env(params, environment())
  }
  x.grid.HS = seq(x.lower, x.upper, length.out = 2*resolution.x + 1)
  y.grid.HS = seq(y.lower, y.upper, length.out = 2*resolution.y + 1)
  if (separate) {
    return(list(x = x.grid.HS, y = y.grid.HS))
  } else {
    s.HS = expand.grid(x.grid.HS, y.grid.HS)
    colnames(s.HS) = c("x", "y")
    return(s.HS)
  }
}

construct.Q.stationary = function(params) {
  required.parameters = c("x.lower", "x.upper", "y.lower", "y.upper",
                          "resolution.x", "resolution.y", "range.parameter", "scale.parameter",
                          "strength.parameter", "direction.parameter")
  missing.parameters = setdiff(required.parameters, names(params))
  if (length(missing.parameters) > 0) {
    stop(sprintf("The following values are missing in params: %s", paste(missing.parameters, collapse = ", ")))
  }
  # Unpack params into function scope (a bit obfuscating, should be fixed)
  list2env(params, environment())
  direction.vector = c(cos(direction.parameter), sin(direction.parameter))
  H = range.parameter^2/8*diag(2) + range.parameter^2/8*strength.parameter*(strength.parameter+2)*direction.vector %*% t(direction.vector)

  stopifnot((x.upper > x.lower) && (y.upper > y.lower))

  n = resolution.x*resolution.y # Total number of cells in grid
  h.x = (x.upper-x.lower)/resolution.x # Step length in x direction
  h.y = (x.upper-x.lower)/resolution.x # Step length in y direction
  V = h.x*h.y # Area of a single grid cell
  # Centers of grid cells
  s = generate.grid.centers(params = params)

  # The suffix HS indicates half-spaced grid
  resolution.x.HS = 2*resolution.x+1
  resolution.y.HS = 2*resolution.y+1
  n.HS = resolution.x.HS*resolution.y.HS
  s.HS = generate.grid.centers.HS(params = params)

  H = matrix(H, nrow=4, ncol=n.HS, byrow = FALSE)
  rownames(H) = c("11", "12", "21", "22")

  dims = c(n, n) # Final dimensions of Q

  # Computing matrix A_H:

  # Start by computing 0-indexed row and column indices for non-zero elements
  # See p. 28 - 29 in Fuglstad et al., 2014 for details
  i = rep(0:(resolution.x-1), times = resolution.y)
  j = rep(0:(resolution.y-1), each = resolution.x)
  ind = j*resolution.x + i

  i.left = (i - 1) %% resolution.x
  ind.left = j*resolution.x + i.left

  i.right = (i + 1) %% resolution.x
  ind.right = j*resolution.x + i.right

  j.up = (j + 1) %% resolution.y
  ind.up = j.up*resolution.x + i

  j.down = (j - 1) %% resolution.y
  ind.down = j.down*resolution.x + i

  A.i = rep(ind, times = 9)
  A.j = c(ind, ind.left, ind.right, ind.up, ind.down,
          j.down*resolution.x + i.left, j.down*resolution.x + i.right,
          j.up*resolution.x + i.left, j.up*resolution.x + i.right)

  # Compute indices on half-spaced grid for picking out correct elements from H
  i.HS = rep(2*(0:(resolution.x-1))+1, times = resolution.y)
  j.HS = rep(2*(0:(resolution.y-1))+1, each = resolution.x)
  ind.HS = j.HS*resolution.x.HS + i.HS + 1

  i.left.HS = (i.HS - 1) %% resolution.x.HS
  ind.left.HS = j.HS*resolution.x.HS + i.left.HS + 1

  i.right.HS = (i.HS + 1) %% resolution.x.HS
  ind.right.HS = j.HS*resolution.x.HS + i.right.HS + 1

  j.up.HS = (j.HS + 1) %% resolution.y.HS
  ind.up.HS = j.up.HS*resolution.x.HS + i.HS + 1

  j.down.HS = (j.HS - 1) %% resolution.y.HS
  ind.down.HS = j.down.HS*resolution.x.HS + i.HS + 1

  # Values of elements, corresponding to order in Fuglstad et. al 2014:
  A.1 = -h.y/h.x*(H["11", ind.right.HS] + H["11", ind.left.HS]) - h.x/h.y*(H["22", ind.up.HS] + H["22", ind.down.HS])

  A.2 = h.y/h.x*H["11", ind.left.HS] - 1/4*(H["12", ind.up.HS] - H["12", ind.down.HS])
  A.3 = h.y/h.x*H["11", ind.right.HS] + 1/4*(H["12", ind.up.HS] - H["12", ind.down.HS])
  A.4 = h.x/h.y*H["22", ind.up.HS] + 1/4*(H["21", ind.right.HS] - H["21", ind.left.HS])
  A.5 = h.x/h.y*H["22", ind.down.HS] - 1/4*(H["21", ind.right.HS] - H["12", ind.left.HS])

  A.6 = 1/4 * (H["12", ind.down.HS] + H["21", ind.left.HS])
  A.7 = -1/4 * (H["12", ind.down.HS] + H["21", ind.right.HS])
  A.8 = -1/4 * (H["12", ind.up.HS] + H["21", ind.left.HS])
  A.9 = 1/4 * (H["12", ind.up.HS] + H["21", ind.right.HS])

  # Make sure to use same order as in A.i and A.j
  A.r = c(A.1, A.2, A.3, A.4, A.5, A.6, A.7, A.8, A.9)

  A.H = Matrix::sparseMatrix(i = A.i, j = A.j, x = A.r, dims = dims, giveCsparse = FALSE, index1 = FALSE)
  A = Matrix::Diagonal(x = rep(V, n)) - A.H
  Q = Matrix::crossprod(A)/V

  # Correct marginal variance
  original.variance = 2/(pi*range.parameter^2*(strength.parameter+1))
  desired.variance = scale.parameter^2
  Q = Q * original.variance/desired.variance
  return(Q)
}

construct.Q.nonstationary = function(params) {
  required.parameters = c("x.lower", "x.upper", "y.lower", "y.upper",
                          "resolution.x", "resolution.y",
                          "range.function", "scale.function",
                          "vector.field")
  missing.parameters = setdiff(required.parameters, names(params))
  if (length(missing.parameters) > 0) {
    stop(sprintf("The following values are missing in params: %s", paste(missing.parameters, collapse = ", ")))
  }
  list2env(params, environment())
  if (!is.function(range.function)) {
    if (is.numeric(range.function) & length(range.function) == 1) {
      range.value = range.function
      range.function = function(x) return(rep(range.value, nrow(x)))
    } else {
      stop("range.function must be either a function or a number.")
    }
  }
  if (!is.function(scale.function)) {
    if (is.numeric(scale.function) & length(scale.function) == 1) {
      scale.value = scale.function
      scale.function = function(x) return(rep(scale.value, nrow(x)))
    } else {
      stop("scale.function must be either a function or a number.")
    }
  }
  if (!is.function(vector.field)) {
    if (is.numeric(vector.field) & length(vector.field) == 2) {
      vector.value = matrix(vector.field, nrow = 1)
      vector.field = function(x) return(vector.value[rep(1, nrow(x)), ])
    } else {
      stop("vector.field must be either a function or a 2D vector.")
    }
  }
  # Unpack params into function scope (a bit obfuscating, should be fixed)

  stopifnot((x.upper > x.lower) && (y.upper > y.lower))

  n = resolution.x*resolution.y # Total number of cells in grid
  h.x = (x.upper-x.lower)/resolution.x # Step length in x direction
  h.y = (x.upper-x.lower)/resolution.x # Step length in y direction
  V = h.x*h.y # Area of a single grid cell
  # Centers of grid cells
  s = generate.grid.centers(params = params)

  # The suffix HS indicates half-spaced grid
  resolution.x.HS = 2*resolution.x+1
  resolution.y.HS = 2*resolution.y+1
  n.HS = resolution.x.HS*resolution.y.HS
  s.HS = generate.grid.centers.HS(params = params)

  range.grid = range.function(s.HS)
  vector.grid = vector.field(s.HS)

  direction.grid = atan2(vector.grid[, 2], vector.grid[, 1])
  strength.grid = sqrt(vector.grid[, 1]^2 + vector.grid[, 2]^2)

  H = 1/8*matrix(c(range.grid^2*(1+strength.grid*(strength.grid+2)*cos(direction.grid)^2),
               range.grid^2*strength.grid*(strength.grid+2)*sin(direction.grid)*cos(direction.grid),
               range.grid^2*strength.grid*(strength.grid+2)*sin(direction.grid)*cos(direction.grid),
               range.grid^2*(1+strength.grid*(strength.grid+2)*sin(direction.grid)^2)),
             byrow = TRUE, nrow = 4)

  rownames(H) = c("11", "12", "21", "22")


  dims = c(n, n) # Final dimensions of Q

  # Computing matrix A_H:

  # Start by computing 0-indexed row and column indices for non-zero elements
  # See p. 28 - 29 in Fuglstad et al., 2014 for details
  i = rep(0:(resolution.x-1), times = resolution.y)
  j = rep(0:(resolution.y-1), each = resolution.x)
  ind = j*resolution.x + i

  i.left = (i - 1) %% resolution.x
  ind.left = j*resolution.x + i.left

  i.right = (i + 1) %% resolution.x
  ind.right = j*resolution.x + i.right

  j.up = (j + 1) %% resolution.y
  ind.up = j.up*resolution.x + i

  j.down = (j - 1) %% resolution.y
  ind.down = j.down*resolution.x + i

  A.i = rep(ind, times = 9)
  A.j = c(ind, ind.left, ind.right, ind.up, ind.down,
          j.down*resolution.x + i.left, j.down*resolution.x + i.right,
          j.up*resolution.x + i.left, j.up*resolution.x + i.right)

  # Compute indices on half-spaced grid for picking out correct elements from H
  i.HS = rep(2*(0:(resolution.x-1))+1, times = resolution.y)
  j.HS = rep(2*(0:(resolution.y-1))+1, each = resolution.x)
  ind.HS = j.HS*resolution.x.HS + i.HS + 1

  i.left.HS = (i.HS - 1) %% resolution.x.HS
  ind.left.HS = j.HS*resolution.x.HS + i.left.HS + 1

  i.right.HS = (i.HS + 1) %% resolution.x.HS
  ind.right.HS = j.HS*resolution.x.HS + i.right.HS + 1

  j.up.HS = (j.HS + 1) %% resolution.y.HS
  ind.up.HS = j.up.HS*resolution.x.HS + i.HS + 1

  j.down.HS = (j.HS - 1) %% resolution.y.HS
  ind.down.HS = j.down.HS*resolution.x.HS + i.HS + 1

  # Values of elements, corresponding to order in Fuglstad et. al 2014:
  A.1 = -h.y/h.x*(H["11", ind.right.HS] + H["11", ind.left.HS]) - h.x/h.y*(H["22", ind.up.HS] + H["22", ind.down.HS])

  A.2 = h.y/h.x*H["11", ind.left.HS] - 1/4*(H["12", ind.up.HS] - H["12", ind.down.HS])
  A.3 = h.y/h.x*H["11", ind.right.HS] + 1/4*(H["12", ind.up.HS] - H["12", ind.down.HS])
  A.4 = h.x/h.y*H["22", ind.up.HS] + 1/4*(H["21", ind.right.HS] - H["21", ind.left.HS])
  A.5 = h.x/h.y*H["22", ind.down.HS] - 1/4*(H["21", ind.right.HS] - H["12", ind.left.HS])

  A.6 = 1/4 * (H["12", ind.down.HS] + H["21", ind.left.HS])
  A.7 = -1/4 * (H["12", ind.down.HS] + H["21", ind.right.HS])
  A.8 = -1/4 * (H["12", ind.up.HS] + H["21", ind.left.HS])
  A.9 = 1/4 * (H["12", ind.up.HS] + H["21", ind.right.HS])

  # Make sure to use same order as in A.i and A.j
  A.r = c(A.1, A.2, A.3, A.4, A.5, A.6, A.7, A.8, A.9)

  A.H = Matrix::sparseMatrix(i = A.i, j = A.j, x = A.r, dims = dims, giveCsparse = FALSE, index1 = FALSE)
  A = Matrix::Diagonal(x = rep(V, n)) - A.H
  Q = Matrix::crossprod(A)/V

  # Correct marginal variance
  half.spaced.indices = get.half.spaced.indices(resolution.x, resolution.y)
  approx.original.stddev.HS = sqrt(2/(pi*range.grid^2*(strength.grid+1)))
  approx.original.stddev = approx.original.stddev.HS[half.spaced.indices]
  desired.stddev = scale.function(s)
  D = Matrix::Diagonal(x = approx.original.stddev / desired.stddev)
  Q = D %*% Q %*% D
  return(Q)
}

get.half.spaced.indices = function(resolution.x, resolution.y) {
  half.spaced.indices = (2*resolution.x+1)*rep(2*(1:resolution.y)-1, each=resolution.x) + 2*rep(1:resolution.x, times=resolution.y)
  return(half.spaced.indices)
}

permuted.cholesky.decomp = function(Q) {
  cholesky.decomp = Matrix::Cholesky(Q, perm=TRUE, LDL=FALSE)
  L.perm = as(cholesky.decomp, "Matrix")
  P = as(cholesky.decomp, "pMatrix")
  return(list(L.perm=L.perm, P=P))
}

generate.gmrf.realization = function(chol.object, Q, n.sim = 1, z = NULL, seed = NULL) {
  if (!missing(Q)) {
    chol.object = permuted.cholesky.decomp(Q)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  n = nrow(chol.object$L.perm)
  u = matrix(NA, nrow=n, ncol=n.sim)
  if (is.null(z)) {
    z = matrix(rnorm(n*n.sim), nrow=n, ncol=n.sim)
  } else {
    z = matrix(z, nrow=n, ncol=n.sim, byrow=FALSE)
  }
  for (i in 1:n.sim) {
    y = Matrix::solve(Matrix::t(chol.object$L.perm), z[, i])
    u[, i] = Matrix::solve(chol.object$P, y)
  }
  if (n.sim == 1 | ncol(z) == 1) {
    return(as.numeric(u))
  } else {
    return(u)
  }
}

standardize = function(values, lower = NULL, upper = NULL) {
  lower = ifelse(is.null(lower), min(values), lower)
  upper = ifelse(is.null(upper), max(values), upper)
  return((values - lower)/(upper - lower))
}
