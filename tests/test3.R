# Test 3 -----------------------------------------------------------------------
if(!require(fdaPDEmixed)){
  devtools::install_github(repo ="aldoclemente/fdaPDEmixed")
}

if(!require(mgcv)){
  install.packages("mgcv")
}

if(!require(rstudioapi)) install.packages("rstudioapi")
# ------------------------------------------------------------------------------
library(fdaPDEmixed)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

data("horseshoe3D", package = "fdaPDEmixed")
mesh=horseshoe3D
plot(mesh)
FEMbasis <- create.FEM.basis(mesh)

Cov1 <- function(x,y,z){
  sin(2*pi*x) * cos(2*pi*y) * z
}

Cov2 <- function(x,y,z){
  cos(2*pi*y) * z^2
}

nnodes <- nrow(mesh$nodes)
betas <- as.matrix(c(3, 0.5))
b_1 <- as.matrix( c(-5, 0.) )
b_2 <- as.matrix(c(5,0))
test.locations <- refine.by.splitting.mesh.3D(mesh)$nodes

test.X1 <- cbind( Cov1(test.locations[,1], test.locations[,2], test.locations[,3]),
                  Cov2(test.locations[,1], test.locations[,2], test.locations[,3]))
test.X2 <- test.X1 

n_sim <- 30
n_obs <- c(100, 250, 500, 1000)

rmse <- function(x,y){
  return(sqrt( mean( x-y )^2 ) ) 
}

# ------------------------------------------------------------------------------

lambda= 10^seq(-1,1,length=15)

# Building folders -------------------------------------------------------------
date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists("data/test_3/")){
  dir.create("data/test_3/")
}

folder.name = paste("data/test_3/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

for(j in 1:length(n_obs)){
  if(j!=4){
  idx <- sample(1:nnodes, size=n_obs[j])
  locations <- mesh$nodes[idx,]
  }else{
    idx <- sample(1:nrow(test.locations), size=n_obs[j])
    locations <- test.locations[idx,]
  }
  nlocs = dim(locations)[1]
  
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2], locations[,3]),
              Cov2(locations[,1], locations[,2], locations[,3]))
  X2 = cbind( Cov1(locations[,1], locations[,2], locations[,3]),
              Cov2(locations[,1], locations[,2], locations[,3]))
  
  func1 = fs.test.3D(x=locations[,1], y=locations[,2],  z = locations[,3])
  func2 = -func1
  
  test.func1 = fs.test.3D(x=test.locations[,1], y=test.locations[,2],  z=test.locations[,3])
  test.func2 = -test.func1
  
  results <- list(
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=2), # (b11, b12)
    b_2 = matrix(0, nrow=n_sim,ncol=2) # (b21, b22)
  )
  
  results <- list(fdaPDE=results, mgcv=results)
  
  errors <- list(
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=2),
    b_2 = matrix(0, nrow=n_sim,ncol=2),
    f_1 = matrix(0, nrow=n_sim,ncol=1),
    f_2 = matrix(0, nrow=n_sim,ncol=1),
    response = matrix(0, nrow=n_sim,ncol=1)
  )
  errors = list(fdaPDE = errors, mgcv = errors)
  
  for(i in 1:n_sim){
    # V == Cov1
    exact1 <- X1%*% betas + func1 + X1%*% b_1 
    exact2 <- X2%*% betas + func2 + X2%*% b_2
    
    exact <- c(exact1, exact2)
    observations <- exact
    observations <- observations + rnorm(nlocs*2, mean=0, sd=0.05*(diff(range(c(func1, func2)))))
    
    X = rbind(X1, X2)
    observations <- matrix(observations, nrow=nlocs, ncol=2)
    
    lambda.selection.criterion = "grid"
    lambda.selection.lossfunction = "GCV"
    DOF.evaluation = "exact"
    # fdaPDE ---------------------------------------------------------------------
    invisible(capture.output(output_fdaPDE <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                            covariates = X, random_effect = c(1,2),
                                                                            FEMbasis = FEMbasis, lambda = lambda, 
                                                                            lambda.selection.criterion = lambda.selection.criterion, 
                                                                            lambda.selection.lossfunction = lambda.selection.lossfunction,
                                                                            DOF.evaluation = DOF.evaluation, FLAG_ITERATIVE = TRUE)))
    
    best_lambda <- output_fdaPDE$bestlambda
    #results
    results$fdaPDE$beta_1[i] <- output_fdaPDE$beta[1, best_lambda]
    results$fdaPDE$beta_2[i] <- output_fdaPDE$beta[2, best_lambda]
    results$fdaPDE$b_1[i,] <- output_fdaPDE$b_i[1:2, best_lambda]
    results$fdaPDE$b_2[i,] <- output_fdaPDE$b_i[3:4, best_lambda]
    
    # errors
    errors$fdaPDE$beta_1[i] <- rmse(output_fdaPDE$beta[1, best_lambda], betas[1])
    errors$fdaPDE$beta_2[i] <- rmse(output_fdaPDE$beta[2, best_lambda], betas[2])
    errors$fdaPDE$b_1[i,] <- c(rmse(results$fdaPDE$b_1[i,1], b_1[1]), rmse(results$fdaPDE$b_1[1,2], b_1[2]))
    errors$fdaPDE$b_2[i,] <- c(rmse(results$fdaPDE$b_2[i,1], b_2[1]), rmse(results$fdaPDE$b_2[1,2], b_2[2]))
    errors$fdaPDE$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations),
                                 test.func1)
    errors$fdaPDE$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations),
                                 test.func2)
    
    
    y_hat1 <- test.X1%*% as.matrix(output_fdaPDE$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations) + 
      test.X1%*% as.matrix(output_fdaPDE$b_i[1:2, best_lambda])
    
    y_hat2 <- test.X2%*% as.matrix(output_fdaPDE$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations) + 
      test.X2%*%as.matrix(output_fdaPDE$b_i[3:4, best_lambda])
    
    test.observations <- c(test.X1%*% betas + test.func1 + test.X1%*% b_1 ,
                           test.X2%*% betas + test.func2 + test.X2%*% b_2 )
    errors$fdaPDE$response[i] <- rmse(c(y_hat1, y_hat2), test.observations)
    
  }
  save(errors, results,
       file = paste0(folder.name, "n_obs_", n_obs[j],".RData"))
  
}

# post-proc --------------------------------------------------------------------
# n_sim <- 30
# n_obs <- c(100, 250, 500, 1000)
# betas <- c(3, 0.5)
fill_col <- viridis::viridis(2, begin=0.25, end=0.95)[1]
# b_1 <- as.matrix( c(-5, 0.) )
# b_2 <- as.matrix(c(5,0))

estimates <- data.frame(beta1 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        beta2 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        b_11 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        b_12 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        b_21 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        b_22 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        beta1_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        beta2_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        f1_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        f2_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        response_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        n_obs = as.factor(rep(n_obs, each=n_sim))
)

#date_ = "..."
#folder.name = paste("data/test_3/",date_,"/",sep="")

for(i in 1:length(n_obs)){
  load(file = paste0(folder.name, "n_obs_", n_obs[i],".RData"))
  estimates$beta1[(1+(i-1)*n_sim):(i*n_sim)] <- results$fdaPDE$beta_1
  estimates$beta2[(1+(i-1)*n_sim):(i*n_sim)] <- results$fdaPDE$beta_2
  estimates$beta1_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$fdaPDE$beta_1
  estimates$beta2_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$fdaPDE$beta_2
  estimates$f1_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$fdaPDE$f_1
  estimates$f2_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$fdaPDE$f_2
  estimates$f3_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$fdaPDE$f_3
  estimates$response_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$fdaPDE$response
}

{
  mai_ = par("mai")
  mai_[2] = mai_[2] + 0.075
  pdf(paste0(folder.name, "test_3.pdf"), family = "serif", width = 7, height = 7)
  par(mai=mai_)
  boxplot(estimates$beta1_rmse ~ estimates$n_obs,ylab="RMSE", xlab="observations",
          main =expression(beta[1]),cex.lab = 2, cex.axis = 2, cex.main = 2,
          col=fill_col)
  par(mai=mai_)
  boxplot(estimates$beta2_rmse ~ estimates$n_obs,ylab="RMSE", xlab="observations",
          main =expression(beta[2]),cex.lab = 2, cex.axis = 2, cex.main = 2,
          col=fill_col)
  par(mai=mai_)
  boxplot(estimates$beta1 ~ estimates$n_obs,
          main =expression(hat(beta)[1]), ylab="", xlab="observations",
          cex.lab = 2, cex.axis = 2, cex.main = 2,
          col=fill_col)
  abline(h=betas[1], lty=2, lwd=3, col="red")
  par(mai=mai_)
  boxplot(estimates$beta2 ~ estimates$n_obs,
          main =expression(hat(beta)[2]), ylab="", xlab="observations",
          cex.lab = 2, cex.axis = 2, cex.main = 2,
          col=fill_col)
  abline(h=betas[2], lty=2, lwd=3, col="red")
  par(mai=mai_)
  boxplot(estimates$f1_rmse ~ estimates$n_obs,
          ylab="RMSE", xlab="observations",
          main =expression(f[1]),cex.lab = 2, cex.axis = 2, cex.main = 2,
          col=fill_col)
  par(mai=mai_)
  boxplot(estimates$f2_rmse ~ estimates$n_obs,
          ylab="RMSE", xlab="observations",
          main =expression(f[2]),cex.lab = 2, cex.axis = 2, cex.main = 2,
          col=fill_col)
  dev.off()
}

# plot mean estimates n = 250  -------------------------------------------------

j = 2
nnodes <- nrow(mesh$nodes)
idx <- sample(1:nnodes, size=n_obs[j])
locations <- mesh$nodes[idx,]
nlocs = dim(locations)[1]

nlocs <- nrow(locations)
X1 = cbind( Cov1(locations[,1], locations[,2], locations[,3]),
            Cov2(locations[,1], locations[,2], locations[,3]))
X2 = cbind( Cov1(locations[,1], locations[,2], locations[,3]),
            Cov2(locations[,1], locations[,2], locations[,3]))

func1 = fs.test.3D(x=locations[,1], y=locations[,2],  z = locations[,3])
func2 = -func1

test.func1 = fs.test.3D(x=test.locations[,1], y=test.locations[,2],  z=test.locations[,3])
test.func2 = -test.func1

estimates <- list(
  f_1 = matrix(0, nrow=nnodes,ncol=1),
  f_2 = matrix(0, nrow=nnodes,ncol=1)
)
lambda= 10^seq(-2,1,by=0.1)
for(i in 1:n_sim){
  exact1 <- X1%*% betas + func1 + X1%*% b_1 
  exact2 <- X2%*% betas + func2 + X2%*% b_2
  
  exact <- c(exact1, exact2)
  observations <- exact
  observations <- observations + rnorm(nlocs*2, mean=0, sd=0.05*(diff(range(c(func1, func2)))))
  
  X = rbind(X1, X2)
  observations <- matrix(observations, nrow=nlocs, ncol=2)
  
  # fdaPDE ---------------------------------------------------------------------
  invisible(capture.output(output_fdaPDE <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                          covariates = X, random_effect = c(1,2),
                                                                          FEMbasis = FEMbasis, lambda = lambda, 
                                                                          lambda.selection.criterion = "grid", 
                                                                          lambda.selection.lossfunction = "GCV",
                                                                          DOF.evaluation = "exact", FLAG_ITERATIVE = TRUE)))
  
  best_lambda <- output_fdaPDE$bestlambda
  estimates$f_1 <- estimates$f_1 + output_fdaPDE$fit.FEM.mixed$coeff[1:nnodes,best_lambda] / n_sim
  estimates$f_2 <- estimates$f_2 + output_fdaPDE$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda] / n_sim
}

save(estimates, file=paste0(folder.name, "n_obs_250_estimates.RData"))