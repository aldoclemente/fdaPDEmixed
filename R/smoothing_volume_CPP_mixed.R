CPP_smooth.volume.FEM.mixed<-function(locations, observations, num_units, FEMbasis, lambda, covariates = NULL, incidence_matrix = NULL, ndim, mydim, BC = NULL, GCV, GCVMETHOD = 2, nrealizations = 100, DOF=TRUE, DOF_matrix=NULL, search, bary.locations, TESTFLAG)
{

  # C++ function for volumetric works with vectors not with matrices

  FEMbasis$mesh$tetrahedrons=c(t(FEMbasis$mesh$tetrahedrons))
  FEMbasis$mesh$nodes=c(t(FEMbasis$mesh$nodes))

  # Indexes in C++ starts from 0, in R from 1, opportune transformation

  FEMbasis$mesh$tetrahedrons=FEMbasis$mesh$tetrahedrons-1

  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(DOF_matrix))
  {
    DOF_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = ndim)
  }

  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)-1

  }

  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }

  ## Set propr type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  observations <- as.vector(observations)
  storage.mode(observations) <- "double"
  storage.mode(FEMbasis$mesh$order) <- "integer"
  storage.mode(FEMbasis$mesh$nnodes) <- "integer"
  storage.mode(FEMbasis$mesh$ntetrahedrons) <- "integer"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$tetrahedrons) <- "integer"
  storage.mode(num_units) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values)  <- "double"
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  DOF <- as.integer(DOF)
  storage.mode(DOF) <-"integer"

  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"
  storage.mode(search) <- "integer"

  TESTFLAG <- as.integer(TESTFLAG)
  storage.mode(TESTFLAG) <-"integer"
  
  ## Call C++ function
  bigsol <- .Call("regression_Laplace_mixed", locations, observations, num_units, FEMbasis$mesh, FEMbasis$mesh$order, mydim, ndim, lambda, covariates,
                  incidence_matrix, BC$BC_indices, BC$BC_values, GCV, GCVMETHOD, nrealizations,  DOF, DOF_matrix, search, bary.locations, TESTFLAG, PACKAGE = "fdaPDEmixed")

  return(bigsol)
}