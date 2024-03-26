# Test 2 -----------------------------------------------------------------------
if(!require(fdaPDEmixed)){
  devtools::install_github(repo ="aldoclemente/fdaPDEmixed")
}

if(!require(mgcv)){
  install.packages("mgcv")
}

if(!require(rstudioapi)) install.packages("rstudioapi")
# ------------------------------------------------------------------------------
library(fdaPDEmixed)
#library(mgcv)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

data(horseshoe2.5D)
mesh=horseshoe2.5D
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
test.locations <- refine.by.splitting.mesh.2.5D(mesh)$nodes
test.locations <- projection.points.2.5D(mesh, test.locations)

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

if( !dir.exists("data/test_2/")){
  dir.create("data/test_2/")
}

folder.name = paste("data/test_2/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}


for(j in 1:length(n_obs)){
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
    # fdaPDE ---------------------------------------------------------------------
    invisible(capture.output(output_fdaPDE <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                            covariates = X, random_effect = c(1,2),
                                                                            FEMbasis = FEMbasis, lambda = lambda, 
                                                                            lambda.selection.criterion = "grid", 
                                                                            lambda.selection.lossfunction = "GCV",
                                                                            DOF.evaluation = "exact", FLAG_ITERATIVE = TRUE)))
    
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
#n_sim <- 30
#n_obs <- c(100, 250, 500, 1000)
#betas <- ...
#b_1 <- ...
#b_2 <- ...
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

#date_ = ""
#folder.name = paste("data/test_2/",date_,"/",sep="")
{
  load(file = paste0(folder.name, "n_obs_100.RData"))
  estimates$beta1[1:n_sim] <- results$fdaPDE$beta_1
  estimates$beta2[1:n_sim] <- results$fdaPDE$beta_2
  estimates$b_11[1:n_sim] <- results$fdaPDE$b_1[,1]
  estimates$b_12[1:n_sim] <- results$fdaPDE$b_1[,2]
  estimates$b_21[1:n_sim] <- results$fdaPDE$b_2[,1]
  estimates$b_22[1:n_sim] <- results$fdaPDE$b_2[,2]
  estimates$beta1_rmse[1:n_sim] <- errors$fdaPDE$beta_1
  estimates$beta2_rmse[1:n_sim] <- errors$fdaPDE$beta_2
  estimates$f1_rmse[1:n_sim] <- errors$fdaPDE$f_1
  estimates$f2_rmse[1:n_sim] <- errors$fdaPDE$f_2
  estimates$response_rmse[1:n_sim] <- errors$fdaPDE$response
  
  load(file = paste0(folder.name,"n_obs_250.RData"))
  estimates$beta1[(n_sim+1):(2*n_sim)] <- results$fdaPDE$beta_1
  estimates$beta2[(n_sim+1):(2*n_sim)] <- results$fdaPDE$beta_2
  estimates$b_11[(n_sim+1):(2*n_sim)] <- results$fdaPDE$b_1[,1]
  estimates$b_12[(n_sim+1):(2*n_sim)] <- results$fdaPDE$b_1[,2]
  estimates$b_21[(n_sim+1):(2*n_sim)] <- results$fdaPDE$b_2[,1]
  estimates$b_22[(n_sim+1):(2*n_sim)] <- results$fdaPDE$b_2[,2]
  estimates$beta1_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$beta_1
  estimates$beta2_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$beta_2
  estimates$f1_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$f_1
  estimates$f2_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$f_2
  estimates$response_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$response
  
  load(file = paste0(folder.name, "n_obs_500.RData"))
  estimates$beta1[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$beta_1
  estimates$beta2[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$beta_2
  estimates$b_11[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$b_1[,1]
  estimates$b_12[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$b_1[,2]
  estimates$b_21[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$b_2[,1]
  estimates$b_22[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$b_2[,2]
  estimates$beta1_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$beta_1
  estimates$beta2_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$beta_2
  estimates$f1_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$f_1
  estimates$f2_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$f_2
  estimates$response_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$response
  
  load(file = paste0(folder.name, "n_obs_1000.RData"))
  estimates$beta1[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$beta_1
  estimates$beta2[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$beta_2
  estimates$b_11[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$b_1[,1]
  estimates$b_12[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$b_1[,2]
  estimates$b_21[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$b_2[,1]
  estimates$b_22[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$b_2[,2]
  estimates$beta1_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$beta_1
  estimates$beta2_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$beta_2
  estimates$f1_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$f_1
  estimates$f2_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$f_2
  estimates$response_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$response
}

{
  mai_ = par("mai")
  mai_[2] = mai_[2] + 0.075
  pdf(paste0(folder.name, "test_2.pdf"), family = "serif", width = 7, height = 7)
  par(mai=mai_)
  boxplot(estimates$beta1_rmse ~ estimates$n_obs,ylab="RMSE", xlab="observations",
          main =expression(beta[1]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
  par(mai=mai_)
  boxplot(estimates$beta2_rmse ~ estimates$n_obs,ylab="RMSE", xlab="observations",
          main =expression(beta[2]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
  par(mai=mai_)
  boxplot(estimates$beta1 ~ estimates$n_obs,
          main =expression(hat(beta)[1]), col="white", ylab="", xlab="observations",
          cex.lab = 2, cex.axis = 2, cex.main = 2)
  abline(h=betas[1], lty=2, lwd=3, col="red")
  par(mai=mai_)
  boxplot(estimates$beta2 ~ estimates$n_obs,
          main =expression(hat(beta)[2]), col="white", ylab="", xlab="observations",
          cex.lab = 2, cex.axis = 2, cex.main = 2)
  abline(h=betas[2], lty=2, lwd=3, col="red")
  par(mai=mai_)
  boxplot(estimates$f1_rmse ~ estimates$n_obs,
          ylab="RMSE", xlab="observations",
          main =expression(f[1]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
  par(mai=mai_)
  boxplot(estimates$f2_rmse ~ estimates$n_obs,
          ylab="RMSE", xlab="observations",
          main =expression(f[2]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
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

library(plotly)
options(warn=-1)

func1 = fs.test.3D(x=mesh$nodes[,1], y=mesh$nodes[,2],  z=mesh$nodes[,3])
func2 = -func1
# FEMbasis.ref <- create.FEM.basis(mesh.ref)

fig1 <- plot.FEM.2.5D(FEM(coeff = func1, FEMbasis))
fig2 <- plot.FEM.2.5D(FEM(coeff = func2, FEMbasis))

fig1 <- fig1 %>% layout(scene = list(
  camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9)))) # ,  dragmode="zoom"
fig2 <- fig2 %>% layout(scene = list(
  camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9))))

save_image(fig1, paste0(folder.name,"true_f1.pdf"))
save_image(fig2, paste0(folder.name,"true_f2.pdf"))

estimates_f1 <- FEM(coeff = estimates$f_1, FEMbasis)
estimates_f2 <- FEM(coeff = estimates$f_2, FEMbasis)

fig1 <- plot.FEM.2.5D(estimates_f1, 
                    limits = compute_limits(FEM(coeff = func1, FEMbasis)))
fig2 <- plot.FEM.2.5D(estimates_f2, 
                    limits = compute_limits(FEM(coeff = func2, FEMbasis)))

fig1 <- fig1 %>% layout(scene = list(
  camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9))))
fig2 <- fig2 %>% layout(scene = list(
  camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9))))
save_image(fig1, paste0(folder.name,"estimate_f1.pdf"))
save_image(fig2, paste0(folder.name,"estimate_f2.pdf"))

fig_cov1 <- plot.FEM.2.5D(FEM(coeff= Cov1(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis))
fig_cov1 <- fig_cov1 %>% layout(scene = list(
  camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9))))

fig_cov2 <- plot.FEM.2.5D(FEM(coeff= Cov2(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis))
fig_cov2 <- fig_cov2 %>% layout(scene = list(
  camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9))))
save_image(fig_cov1, paste0(folder.name,"cov_1.pdf"))
save_image(fig_cov2, paste0(folder.name,"cov_2.pdf"))

fig_mesh <- plot.mesh.2.5D(mesh) %>% 
  layout(scene = list(camera = list(
    eye = list(x = 1.5, y = -2.25,  z = 0.9))))

save_image(fig_mesh, paste0(folder.name, "Cshaped_surface.pdf"))
