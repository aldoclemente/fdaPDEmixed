f(!require(fdaPDEmixed)){
  devtools::install_github(repo ="aldoclemente/fdaPDEmixed")
}

if(!require(rstudioapi)) install.packages("rstudioapi")
# ------------------------------------------------------------------------------
library(fdaPDEmixed)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

set.seed(0)
foldername <- "data/"
if(!dir.exists(foldername)) dir.create(foldername)

foldername <- paste0(foldername,"test_1/")
if(!dir.exists(foldername)) dir.create(foldername)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

test.locations <- refine.mesh.2D(mesh, maximum_area = 0.0025, minimum_angle = 30)$nodes

Cov1 <- function(x,y){
  sin(2*pi*x) * cos(2*pi*y)
}

fs.test.time<-function(x,y,t=y) {
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  
  res=numeric(length =length(x))
  
  for(i in 1:length(x)) {
    if(x[i]>=0 && y[i]>0)
      res[i]=cos(t[i])*(0.25*pi+x[i])+(y[i]-0.5)^2
    
    if(x[i]>=0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    
    if(x[i]<0 && y[i]>0)
      res[i]=cos(t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    
    if(x[i]<0 && y[i]<=0)
      res[i]=cos(2*t[i])*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}
nnodes <- nrow(mesh$nodes)
betas <- as.matrix(c(3,0.5))
b <- as.matrix( c(-5., 0, 5.) )

n_obs <- c(100, 250, 500, 1000)

set.seed(0)
for(j in 1:length(n_obs)){
  path <- paste0(foldername, n_obs[j], "/")
  if(!dir.exists(path)) dir.create(path)
  
  locations <- sample_locations.mgcv(mesh, n_obs[j])
  
  write.csv(format(locations, digits=16), file=paste0(path, "locations_1.csv"))
  write.csv(format(rbind(locations,locations,locations), digits=16), 
            file=paste0(path, "locations.csv"))
  
  nlocs = dim(locations)[1]
  
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  
  write.csv(format(X1, digits=16), file=paste0(path, "covariates_1.csv"))
  write.csv(format(X2, digits=16), file=paste0(path, "covariates_2.csv"))
  write.csv(format(X3, digits=16), file=paste0(path, "covariates_3.csv"))
  write.csv(format(rbind(X1,X2,X3), digits=16), file=paste0(path, "covariates.csv"))
  
  func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
  
  obs1 <- X1%*% betas + func1 + X1[,1] * b[1] 
  
  obs2 <- X2%*% betas + func2 + X2[,1] * b[2]
  
  obs3 <- X3%*% betas + func3 + X3[,1] * b[3] 
  
  observations <- matrix(c(obs1, obs2, obs3))
  observations <- observations + rnorm(nlocs*3, mean=0, sd=0.05*(diff(range(c(func1, func2, func3)))))
  
  write.csv(format(observations[1:nlocs,], 
                   digits=16), file=paste0(path, "observations_1.csv"))
  write.csv(format(observations[(nlocs+1):(2*nlocs)], 
                   digits=16), file=paste0(path, "observations_2.csv"))
  write.csv(format(observations[(2*nlocs+1):(3*nlocs)], 
                   digits=16), file=paste0(path, "observations_3.csv"))
  write.csv(format(observations, digits=16), file=paste0(path, "observations.csv"))
  
  observations <- matrix(observations, nrow=nlocs, ncol=3)
  X = rbind(X1, X2, X3)
  
  write.csv(format(X1[,2], digits=16), file=paste0(path, "W_1.csv"))
  write.csv(format(X2[,2], digits=16), file=paste0(path, "W_2.csv"))
  write.csv(format(X3[,2], digits=16), file=paste0(path, "W_3.csv"))
  write.csv(format(X[,2], digits=16), 
            file=paste0(path, "W.csv"))
  
  write.csv(format(X1[,1], digits=16), file=paste0(path, "V_1.csv"))
  write.csv(format(X2[,1], digits=16), file=paste0(path, "V_2.csv"))
  write.csv(format(X3[,1], digits=16), file=paste0(path, "V_3.csv"))
  write.csv(format(X[,1], digits=16), 
            file=paste0(path, "V.csv"))
  
  # multi-domain "design_matrix"
  write.csv(format(
            cbind(as.matrix(X[,2]), 
                  as.matrix(c(X[1:nlocs,1],            rep(0, times=2*nlocs))),
                  as.matrix(c(rep(0, times=nlocs),X[(nlocs+1):(2*nlocs),1],  rep(0, times=nlocs))),
                  as.matrix(c(rep(0, times=2*nlocs), X[(2*nlocs+1):(3*nlocs),1]))), digits=16), 
            file=paste0(path, "X.csv"))
  
  
  # fdaPDE ---------------------------------------------------------------------
  invisible(capture.output(output_fdaPDE <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                          covariates = X, random_effect = c(1),
                                                                          FEMbasis = FEMbasis, lambda=1,
                                                                          FLAG_ITERATIVE = FALSE)))
  
  write.csv(format(output_fdaPDE$beta, digits=16), file = paste0(path,"beta_hat.csv"))
  write.csv(format(output_fdaPDE$b_i, digits=16), file = paste0(path,"b_hat.csv"))
  write.csv(format(output_fdaPDE$fit.FEM.mixed$coeff[1:nnodes,], 
                   digits=16), file = paste0(path,"f_1_hat.csv"))
  write.csv(format(output_fdaPDE$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),], 
                   digits=16), file = paste0(path,"f_2_hat.csv"))
  write.csv(format(output_fdaPDE$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),], 
                   digits=16), file = paste0(path,"f_3_hat.csv"))
  write.csv(format(output_fdaPDE$fit.FEM.mixed$coeff, 
                   digits=16), file = paste0(path,"f_hat.csv"))
  
  
  func1 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nlocs))
  
  write.csv(format(as.matrix(func1), digits=16), file = paste0(path, "f_1.csv")) # esatta !
  write.csv(format(as.matrix(func2), digits=16), file = paste0(path, "f_2.csv"))
  write.csv(format(as.matrix(func3), digits=16), file = paste0(path, "f_3.csv"))
  write.csv(format(as.matrix(c(func1, func2, func3)), 
                   digits=16), file = paste0(path, "f.csv"))     # esatta !
  
}
