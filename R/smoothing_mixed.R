#' @export
smooth.FEM.mixed <- function(locations = NULL, observations, FEMbasis,
                             covariates, random_effect = NULL, PDE_parameters = NULL, BC = NULL,
                             incidence_matrix = NULL, areal.data.avg = TRUE,
                             FLAG_ITERATIVE = FALSE, threshold = 10^(-4), threshold_residual = 1e-8, max.steps = 50,
                             search = "tree", bary.locations = NULL, lambda = NULL,
                             lambda.selection.criterion = "grid", DOF.evaluation = NULL, lambda.selection.lossfunction = NULL,
                             DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0,
                             nrealizations = 100, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05,
                             verbose = FALSE, anderson_memory = 10L) {
  if (!is.logical(verbose)) {
    stop("'verbose' is not logical")
  }
  if (class(FEMbasis$mesh) == "mesh.2D") {
    ndim <- 2
    mydim <- 2
  } else if (class(FEMbasis$mesh) == "mesh.2.5D") {
    ndim <- 3
    mydim <- 2
  } else if (class(FEMbasis$mesh) == "mesh.3D") {
    ndim <- 3
    mydim <- 3
  } else {
    stop("Unknown mesh class")
  }
  ##################### Checking parameters, sizes and conversion ################################

  # Preliminary consistency of optimization parameters
  if (lambda.selection.criterion == "grid") {
    optim <- 0
  } else if (lambda.selection.criterion == "newton") {
    optim <- 1
  } else if (lambda.selection.criterion == "newton_fd") {
    optim <- 2
  } else {
    stop("'lambda.selection.criterion' must belong to the following list: 'none', 'grid', 'newton', 'newton_fd'.")
  }

  if (is.null(DOF.evaluation)) {
    optim <- c(optim, 0)
  } else if (DOF.evaluation == "stochastic") {
    optim <- c(optim, 1)
  } else if (DOF.evaluation == "exact") {
    optim <- c(optim, 2)
  } else {
    stop("'DOF.evaluation' must be NULL, 'stochastic' or 'exact'.")
  }

  if (is.null(lambda.selection.lossfunction)) {
    optim <- c(optim, 0)
  } else if (lambda.selection.lossfunction == "GCV") {
    optim <- c(optim, 1)
  } else {
    stop("'lambda.selection.lossfunction' has to be 'GCV'.")
  }

  # --> General consistency rules
  if (optim[2] != 0 & optim[3] != 1) {
    warning("Dof are computed, setting 'lambda.selection.lossfunction' to 'GCV'")
    optim[3] <- 1
  }
  if (optim[1] == 1 & optim[2] != 2) {
    warning("This method needs evaluate DOF in an 'exact' way, selecting 'DOF.evaluation'='exact'")
    optim[2] <- 2
  }
  if (!is.null(BC) & optim[1] == 1) {
    warning("'newton' 'lambda.selection.criterion' can't be performed with non-NULL boundary conditions, using 'newton_fd' instead")
    optim[1] <- 2
  }
  if ((optim[1] == 2 & optim[2] == 0) || (optim[1] == 0 & optim[2] == 0 & optim[3] == 1 & is.null(DOF.matrix))) {
    warning("This method needs evaluate DOF, selecting 'DOF.evaluation'='stochastic'")
    optim[2] <- 1
  }
  if (optim[1] != 0 & optim[3] == 0) {
    warning("An optimized method needs a loss function to perform the evaluation, selecting 'lambda.selection.lossfunction' as 'GCV'")
    optim[3] <- 1
  }

  if (is.null(lambda) & optim[1] == 0) {
    warning("the lambda passed is NULL, passing to default optimized methods")
    optim <- c(2, 1, 1)
  }


  if (any(lambda <= 0)) {
    stop("'lambda' can not be less than or equal to 0")
  }


  # Search algorithm
  if (search == "naive") {
    search <- 1
  } else if (search == "tree") {
    search <- 2
  } else if (search == "walking" & class(FEMbasis$mesh) == "mesh.2.5D") {
    stop("walking search is not available for mesh class mesh.2.5D.")
  } else if (search == "walking" & class(FEMbasis$mesh) != "mesh.2.5D") {
    search <- 3
  } else {
    stop("'search' must must belong to the following list: 'naive', 'tree' or 'walking'.")
  }

  # If locations is null but bary.locations is not null, use the locations in bary.locations
  if (is.null(locations) & !is.null(bary.locations)) {
    locations <- bary.locations$locations
    locations <- as.matrix(locations)
  }

  ## Converting to format for internal usage
  if (!is.null(locations)) {
    locations <- as.matrix(locations)
  }

  observations <- as.matrix(observations)
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
  }
  if (!is.null(incidence_matrix)) {
    incidence_matrix <- as.matrix(incidence_matrix)
  }
  if (!is.null(BC)) {
    BC$BC_indices <- as.matrix(BC$BC_indices)
    BC$BC_values <- as.matrix(BC$BC_values)
  }
  if (!is.null(lambda)) {
    lambda <- as.matrix(lambda)
  }
  if (!is.null(DOF.matrix)) {
    DOF.matrix <- as.matrix(DOF.matrix)
  }

  # Checking random effect
  if (is.null(random_effect)) warning(" you called a mixed solver with no random effects, it is therefore a fixed-effects model")
  random_effect <- sort(as.vector(random_effect))
  if (!is.null(random_effect)) {
    if (ncol(covariates) < length(random_effect)) {
      stop("'random_effect' has out of range index")
    }

    # check whether the index in random_effect vectors are appropriate
    cov_ind <- seq_len(dim(covariates)[2])
    for (i in random_effect) {
      if (!(i %in% cov_ind)) {
        stop("'random_effect' has out of range index")
      }
    }
  }

  space_varying <- checkSmoothingParameters_mixed(locations = locations, observations = observations, FEMbasis = FEMbasis, covariates = covariates, PDE_parameters = PDE_parameters, BC = BC, incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg, search = search, bary.locations = bary.locations, optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed, DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance, FLAG_ITERATIVE = FLAG_ITERATIVE, threshold = threshold, max.steps = max.steps, threshold_residual = threshold_residual)

  # if I have PDE non-sv case I need (constant) matrices as parameters
  if (!is.null(PDE_parameters) & space_varying == FALSE) {
    PDE_parameters$K <- as.matrix(PDE_parameters$K)
    PDE_parameters$b <- as.matrix(PDE_parameters$b)
    PDE_parameters$c <- as.matrix(PDE_parameters$c)
  }

  checkSmoothingParametersSize(
    locations = locations, observations = observations, FEMbasis = FEMbasis,
    covariates = covariates, PDE_parameters = PDE_parameters, incidence_matrix = incidence_matrix,
    BC = BC, space_varying = space_varying, ndim = ndim, mydim = mydim,
    lambda = lambda, DOF.matrix = DOF.matrix
  )

  # Check covariates sizes (it is better to keep it out of the above function for generalization of the code)
  if (nrow(covariates) != nrow(observations) * ncol(observations)) {
    stop("'covariates' and 'observations' have incompatible size;")
  }


  # Check whether the locations coincide with the mesh nodes (should be put after all the validations)
  if (!is.null(locations)) {
    if (dim(locations)[1] == dim(FEMbasis$mesh$nodes)[1] & dim(locations)[2] == dim(FEMbasis$mesh$nodes)[2]) {
      sum1 <- 0
      sum2 <- 0
      for (i in 1:nrow(locations)) {
        sum1 <- sum1 + abs(locations[i, 1] - FEMbasis$mesh$nodes[i, 1])
        sum2 <- sum2 + abs(locations[i, 2] - FEMbasis$mesh$nodes[i, 2])
      }
      if (sum1 == 0 & sum2 == 0) {
        message("No search algorithm is used because the locations coincide with the nodes.")
        locations <- NULL # In principle, R uses pass-by-value semantics in its function calls. So put ouside of checkSmoothingParameters function.
      }
    }
  }

  num_units <- dim(observations)[2] # number of statistical units, m

  # transform matrix data to vector data
  observations <- as.vector(observations)

  if ((class(FEMbasis$mesh) == "mesh.2D" | class(FEMbasis$mesh) == "mesh.2.5D") & is.null(PDE_parameters)) {
    bigsol <- NULL
    print("C++ Code Execution")

    bigsol <- CPP_smooth.FEM.mixed(
      locations = locations, observations = observations, FEMbasis = FEMbasis,
      covariates = covariates, ndim = ndim, mydim = mydim, BC = BC, num_units = num_units, random_effect = random_effect,
      incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
      search = search, bary.locations = bary.locations,
      optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
      DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance, FLAG_ITERATIVE = FLAG_ITERATIVE, threshold = threshold, max.steps = max.steps, threshold_residual = threshold_residual, verbose = verbose,
      anderson_memory = anderson_memory
    )
  } else if (class(FEMbasis$mesh) == "mesh.2D" & !is.null(PDE_parameters) & space_varying == FALSE) {
    bigsol <- NULL
    stop("still to be reimplemented")
  } else if (class(FEMbasis$mesh) == "mesh.2D" & !is.null(PDE_parameters) & space_varying == TRUE) {
    bigsol <- NULL
    stop("still to be reimplemented")
  }
  else if (class(FEMbasis$mesh) == "mesh.3D" & is.null(PDE_parameters)) {
    bigsol <- NULL
    print("C++ Code Execution")
    
    bigsol <- CPP_smooth.volume.FEM.mixed(
      locations = locations, observations = observations, FEMbasis = FEMbasis,
      covariates = covariates, ndim = ndim, mydim = mydim, BC = BC, num_units = num_units, random_effect = random_effect,
      incidence_matrix = incidence_matrix, areal.data.avg = areal.data.avg,
      search = search, bary.locations = bary.locations,
      optim = optim, lambda = lambda, DOF.stochastic.realizations = DOF.stochastic.realizations, DOF.stochastic.seed = DOF.stochastic.seed,
      DOF.matrix = DOF.matrix, GCV.inflation.factor = GCV.inflation.factor, lambda.optimization.tolerance = lambda.optimization.tolerance, FLAG_ITERATIVE = FLAG_ITERATIVE, threshold = threshold, max.steps = max.steps, threshold_residual = threshold_residual, verbose = verbose,
      anderson_memory = anderson_memory
    )
  }

  f <- bigsol[[1]][1:(num_units * nrow(FEMbasis$mesh$nodes)), ]
  g <- bigsol[[1]][(num_units * nrow(FEMbasis$mesh$nodes) + 1):(2 * num_units * nrow(FEMbasis$mesh$nodes)), ]
  dof <- bigsol[[2]]
  GCV_ <- bigsol[[3]]
  bestlambda <- bigsol[[4]] + 1

  # Start coefficient conversion
  p <- length(random_effect) # num of random-effect coeff
  q <- dim(covariates)[2] # num of common-effect coeff
  matrixcoeff <- matrix(data = bigsol[[5]], nrow = q - p + num_units * p, ncol = ifelse(is.null(lambda), 1, length(lambda))) # implementative ver. (length: (q-p) + m*p)

  # convert into official coeff (length: q + m*p)
  if (p < q) {
    # random-effect as subset
    betaPrime <- matrixcoeff[1:(q - p), , drop = FALSE]
    b_iPrime <- matrixcoeff[-(1:(q - p)), , drop = FALSE] # split matrixcoeff into 2 matrices
  }
  else if (p == q) {
    # random-effect as full set
    b_iPrime <- matrixcoeff
  }

  # convert fixed-effect and random effect
  beta <- matrix(0, nrow = q, ncol = ifelse(is.null(lambda), 1, length(lambda)))
  indBeta <- 1
  indBi <- 1

  for (i in 1:q){
    if (!is.element(i, random_effect)) { # beta as it is
      beta[i, ] <- betaPrime[indBeta, , drop = FALSE]
      indBeta <- indBeta + 1
    }else { # convert beta prime to original prime
      temp <- numeric(ifelse(is.null(lambda), 1, length(lambda))) # it initializes all to 0
      for (j in 1:num_units){
        temp <- temp + b_iPrime[indBi + (j - 1) * p, ]
      }
      beta[i, ] <- temp / num_units
      indBi <- indBi + 1
    }
  }
  b_i <- matrix(0, nrow = num_units * p, ncol = ifelse(is.null(lambda), 1, length(lambda)))
  if (p != 0) {
    indRanEff <- 1 # this index will be cycled according to random_effect elements
    for (i in 1:(num_units * p)){
      b_i[i, ] <- b_iPrime[i, ] - beta[random_effect[ifelse(indRanEff != 0, indRanEff, p)], ]
      indRanEff <- (indRanEff + 1) %% p
    }

    # change the name of the rows
    rname <- c()
    for (i in 1:num_units) {
      temp <- paste("b_", as.character(i), sep = "")
      for (j in 1:p) {
        temp2 <- paste(temp, as.character(j), sep = "")
        rname <- c(rname, temp2)
      }
    }
    rownames(b_i) <- rname
  }

  # End of conversion


  # Save information of Tree Mesh
  tree_mesh <- list(
    treelev = bigsol[[6]][1],
    header_orig = bigsol[[7]],
    header_scale = bigsol[[8]],
    node_id = bigsol[[9]][, 1],
    node_left_child = bigsol[[9]][, 2],
    node_right_child = bigsol[[9]][, 3],
    node_box = bigsol[[10]]
  )

  # Reconstruct FEMbasis with treeed<-f mesh
  mesh.class <- class(FEMbasis$mesh)
  if (is.null(FEMbasis$mesh$treelev)) { # if doesn't exist the tree information
    FEMbasis$mesh <- append(FEMbasis$mesh, tree_mesh)
  } # if already exist the tree information, don't append
  class(FEMbasis$mesh) <- mesh.class


  # Make Functional objects
  fit.FEM.mixed <- FEM.mixed(f, num_units, FEMbasis)
  PDEmisfit.FEM.mixed <- FEM.mixed(g, num_units, FEMbasis)


  # Save information of Barycenter
  if (is.null(bary.locations)) {
    bary.locations <- list(locations = locations, element_ids = bigsol[[11]], barycenters = bigsol[[12]])
  }
  class(bary.locations) <- "bary.locations"

  # iterative method part
  if (FLAG_ITERATIVE) {
    iterations <- bigsol[[13]]

    residuals <- list()
    start <- 0
    for (size in iterations) {
      residuals <- c(residuals, list(bigsol[[14]][(start + 1):(start + size + 1)]))
      start <- start + size + 1
    }
  }
  else {
    iterations <- NaN
    residuals <- NaN
  }

  # Prepare return list
  reslist <- NULL
  if (optim[3] == 1) {
    pure_obs_len <- length(which(!is.na(observations)))
    stderr <- sqrt(GCV_ * (pure_obs_len - dof) / pure_obs_len)

    reslist <- list(
      fit.FEM.mixed = fit.FEM.mixed, PDEmisfit.FEM.mixed = PDEmisfit.FEM.mixed,
      beta = beta, b_i = b_i, edf = dof, GCV = GCV_, stderr = stderr, bestlambda = bestlambda, bary.locations = bary.locations, iterations = iterations, residuals = residuals
    )
  } else {
    reslist <- list(fit.FEM.mixed = fit.FEM.mixed, PDEmisfit.FEM.mixed = PDEmisfit.FEM.mixed, beta = beta, b_i = b_i, bary.locations = bary.locations, iterations = iterations, residuals = residuals)
  }


  return(reslist)
}