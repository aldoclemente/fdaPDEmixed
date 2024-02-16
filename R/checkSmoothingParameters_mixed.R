# It's almost identical to the standard case, maybe can be made ugual?

checkSmoothingParameters_mixed <- function(locations = NULL, observations, FEMbasis, covariates, PDE_parameters = NULL, BC = NULL, incidence_matrix = NULL, areal.data.avg = TRUE, search = "tree", bary.locations = NULL, optim, lambda = NULL, DOF.stochastic.realizations = 100, DOF.stochastic.seed = 0, DOF.matrix = NULL, GCV.inflation.factor = 1, lambda.optimization.tolerance = 0.05, FLAG_ITERATIVE, threshold, max.steps, threshold_residual) {
  #################### Full Consistency Parameter Check #########################
  # Mesh type and methods
  if (is.null(FEMbasis)) {
    stop("FEMbasis required;  is NULL.")
  }
  if (class(FEMbasis) != "FEMbasis") {
    stop("'FEMbasis' is not class 'FEMbasis'")
  }

  if (class(FEMbasis$mesh) != "mesh.2D" & class(FEMbasis$mesh) != "mesh.2.5D" & class(FEMbasis$mesh) != "mesh.3D") {
    stop("Unknown mesh class")
  }

  if ((class(FEMbasis$mesh) == "mesh.2.5D" || class(FEMbasis$mesh) == "mesh.3D") & !is.null(PDE_parameters)) {
    stop("For mesh classes different from mesh.2D, anysotropic regularization is not yet implemented.
         Use Laplacian regularization instead")
  }

  # Locations & observations & areal
  if (!is.null(locations)) {
    if (any(is.na(locations))) {
      stop("Missing values not admitted in 'locations'.")
    }
  }

  if (is.null(observations)) {
    stop("observations required;  is NULL.")
  }

  if (is.null(covariates)) {
    stop("covariates required;  is NULL.")
  }


  if (!is.null(locations) & !is.null(incidence_matrix)) {
    stop("Both 'locations' and 'incidence_matrix' are given. In case of pointwise data, set 'incidence_matrix to NULL. In case of areal data, set 'locations' to NULL.")
  }

  if (any(incidence_matrix != 0 & incidence_matrix != 1)) {
    stop("Value different than 0 or 1 in 'incidence_matrix'.")
  }

  if (is.null(areal.data.avg)) {
    stop("'areal.data.avg' required; is NULL.")
  }
  if (!is.logical(areal.data.avg)) {
    stop("'areal.data.avg' is not logical")
  }

  # PDE_parameters
  if (!is.null(PDE_parameters)) {
    if (is.null(PDE_parameters$K)) {
      stop("'K' required in PDE_parameters;  is NULL.")
    }
    if (is.null(PDE_parameters$b)) {
      stop("'b' required in PDE_parameters;  is NULL.")
    }
    if (is.null(PDE_parameters$c)) {
      stop("'c' required in PDE_parameters;  is NULL.")
    }
  }

  space_varying <- FALSE

  if (!is.null(PDE_parameters$u)) {
    space_varying <- TRUE

    message("Smoothing: anysotropic and non-stationary case")

    if (!is.function(PDE_parameters$K)) {
      stop("'K' in 'PDE_parameters' is not a function")
    }
    if (!is.function(PDE_parameters$b)) {
      stop("'b' in 'PDE_parameters' is not a function")
    }
    if (!is.function(PDE_parameters$c)) {
      stop("'c' in 'PDE_parameters' is not a function")
    }
    if (!is.function(PDE_parameters$u)) {
      stop("'u' in 'PDE_parameters' is not a function")
    }
  }
  else if (!is.null(PDE_parameters)) {
    message("Smoothing: anysotropic and stationary case")
  }

  if (is.null(PDE_parameters)) {
    message("Smoothing: isotropic and stationary case")
  }

  # Boundary Conditions [BC]
  if (!is.null(BC)) {
    if (is.null(BC$BC_indices)) {
      stop("'BC_indices' required in BC;  is NULL.")
    }
    if (is.null(BC$BC_values)) {
      stop("'BC_indices' required in BC;  is NULL.")
    }
  }


  # Check the locations in 'bary.locations' and 'locations' are the same
  if (!is.null(bary.locations) & !is.null(locations)) {
    flag <- TRUE
    for (i in 1:nrow(locations)) {
      if (!(locations[i, 1] == bary.locations$locations[i, 1] & locations[i, 2] == bary.locations$locations[i, 2])) {
        flag <- FALSE
        break
      }
    }
    if (flag == FALSE) {
      stop("Locations are not same as the one in barycenter information.")
    }
  } # end of bary.locations

  # Optimization
  if (optim[1] == 1 & optim[2] == 1) {
    stop("Newton method can only be applied in a 'DOF.evaluation' = 'exact' context")
  }

  # --> Lambda related
  if (optim[1] == 0 & is.null(lambda)) {
    stop("'lambda' required for 'lambda.selection.criterion' = 'grid'; now is NULL.")
  }
  if (optim[1] != 0 & !is.null(lambda)) {
    if (length(lambda) > 1) {
      warning("In optimized methods 'lambda' is the initial value, all terms following the first will be discarded")
    }
  }

  # --> Stochastic related data
  if (!is.numeric(DOF.stochastic.realizations)) {
    stop("'DOF.stochastic.realizations' must be a positive integer")
  } else if (DOF.stochastic.realizations < 1) {
    stop("'DOF.stochastic.realizations' must be a positive integer")
  }

  if (!is.numeric(DOF.stochastic.seed)) {
    stop("'DOF.stochastic.seed' must be a non-negative integer")
  } else if (DOF.stochastic.seed < 0) {
    stop("'DOF.stochastic.seed' must be a non-negative integer")
  }

  if ((DOF.stochastic.realizations != 100 || DOF.stochastic.seed != 0) & optim[2] != 1) {
    warning("'DOF.stochastic.realizations' and 'DOF.stochastic.seed' are used just with 'DOF.evaluation' = 'stochastic'")
  }

  # --> GCV.inflation.factor related
  if (is.null(GCV.inflation.factor)) {
    stop("'GCV.inflation.factor' required;  is NULL.")
  } else if (!is.numeric(GCV.inflation.factor)) {
    stop("'GCV.inflation.factor' must be a non-negative real")
  } else if (GCV.inflation.factor < 0) {
    stop("'GCV.inflation.factor' must be a non-negative real")
  }
  if (GCV.inflation.factor != 1 & optim[3] != 1) {
    warning("'GCV' not selected as 'loss function', 'GCV.inflation.factor' unused")
  }

  # --> DOF.matrix related
  if (!is.null(DOF.matrix)) {
    if (optim[1] != 0) {
      stop("An optimization method needs DOF to be computed during the call, please set 'DOF.matrix' to 'NULL")
    }
    if (optim[2] != 0) {
      stop("'DOF.matrix' is passed to the function, 'DOF.evaluation' should be NULL")
    }
    if (optim[3] != 1) {
      warning("'GCV' is not the 'lambda.selection.lossfunction'. DOF.matrix is passed but GCV is not computed")
    }
  }
  if (is.null(DOF.matrix) & optim[2] == 0 & optim[3] == 1) {
    stop("Either 'DOF.matrix' different from NULL or 'DOF.evaluation' different from NULL, otherwise 'lambda.selection.lossfunction' = 'GCV' can't be computed")
  }

  # --> TOLERANCE
  if (!is.numeric(lambda.optimization.tolerance)) {
    stop("'stopping_criterion_tol' must be a numeric percentage between 0 and 1")
  } else if (lambda.optimization.tolerance >= 1 || lambda.optimization.tolerance <= 0) {
    stop("'stopping_criterion_tol' must be a numeric percentage between 0 and 1")
  }

  if (optim[1] == 0 & lambda.optimization.tolerance != 0.05) {
    warning("'lambda.optimization.tolerance' is not used in grid evaluation")
  }

  if (is.null(FLAG_ITERATIVE)) {
    stop("FLAG_ITERATIVE required;  is NULL.")
  }
  if (!is.logical(FLAG_ITERATIVE)) {
    stop("'FLAG_ITERATIVE' is not logical")
  }

  # Check max.steps and threshold for the iterative method
  if (!all.equal(max.steps, as.integer(max.steps)) || max.steps <= 0) {
    stop("'max.steps' must be a positive integer.")
  }

  #   if( !is.numeric(threshold) || threshold <= 0)
  #     stop("'threshold' must be a real positive")
  # print(" hhhhhhaaa?????")

  #   if( !is.numeric(threshold_residual) || threshold_residual <= 0)
  #     stop("'threshold_residual' must be a real positive")
  # print(" aaa?????")

  # Return information
  return(space_varying)
}