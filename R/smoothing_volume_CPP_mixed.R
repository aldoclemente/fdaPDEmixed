CPP_smooth.volume.FEM.mixed<-function(locations, observations, FEMbasis, covariates, ndim, mydim, BC, num_units, random_effect, incidence_matrix, areal.data.avg, search, bary.locations, optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, FLAG_ITERATIVE, threshold, max.steps, threshold_residual, verbose, anderson_memory){
  
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  FEMbasis$mesh$tetrahedrons = FEMbasis$mesh$tetrahedrons - 1
  FEMbasis$mesh$faces = FEMbasis$mesh$faces - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1

  random_effect <- random_effect - 1
  
  if (is.null(covariates)) {
    covariates <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(DOF.matrix)) {
    DOF.matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(locations)) {
    locations <- matrix(nrow = 0, ncol = ndim)
  }
  
  
  if (is.null(incidence_matrix)) {
    incidence_matrix <- matrix(nrow = 0, ncol = 1)
  }
  
  if (is.null(BC$BC_indices)) {
    BC$BC_indices <- vector(length = 0)
  } else {
    BC$BC_indices <- as.vector(BC$BC_indices) - 1
  }
  
  if (is.null(BC$BC_values)) {
    BC$BC_values <- vector(length = 0)
  } else {
    BC$BC_values <- as.vector(BC$BC_values)
  }
  
  if (is.null(lambda)) {
    lambda <- vector(length = 0)
  } else {
    lambda <- as.vector(lambda)
  }

  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  observations <- as.vector(observations)
  storage.mode(locations) <- "double"
  storage.mode(observations) <- "double"
  storage.mode(random_effect) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$faces) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  storage.mode(num_units) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  areal.data.avg <- as.integer(areal.data.avg)
  storage.mode(areal.data.avg) <- "integer"
  storage.mode(search) <- "integer"
  
  storage.mode(optim) <- "integer"
  storage.mode(lambda) <- "double"
  DOF.matrix <- as.matrix(DOF.matrix)
  storage.mode(DOF.matrix) <- "double"
  storage.mode(DOF.stochastic.realizations) <- "integer"
  storage.mode(DOF.stochastic.seed) <- "integer"
  storage.mode(GCV.inflation.factor) <- "double"
  storage.mode(lambda.optimization.tolerance) <- "double"
  FLAG_ITERATIVE <- as.integer(FLAG_ITERATIVE)
  storage.mode(FLAG_ITERATIVE) <- "integer"
  verbose <- as.integer(verbose)
  storage.mode(verbose) <- "integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(threshold) <- "double"
  storage.mode(threshold_residual) <- "double"
  storage.mode(anderson_memory) <- "integer"
  ## Call C++ function
  bigsol <- .Call("regression_Laplace_mixed", locations, bary.locations, observations, num_units, random_effect, FEMbasis$mesh, FEMbasis$order, mydim, ndim, covariates, BC$BC_indices, BC$BC_values, incidence_matrix, areal.data.avg, search, optim, lambda, DOF.stochastic.realizations, DOF.stochastic.seed, DOF.matrix, GCV.inflation.factor, lambda.optimization.tolerance, PACKAGE = "fdaPDEmixed", 
                  FLAG_ITERATIVE, threshold, max.steps, threshold_residual, verbose, anderson_memory)
  return(bigsol)
}