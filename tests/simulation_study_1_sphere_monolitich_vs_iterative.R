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
library(plotly)
#library(mgcv)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

load("data/sphere2.5D.RData")
mesh=sphere2.5D
FEMbasis <- create.FEM.basis(mesh)

fdaPDE:::plot.mesh.2.5D(mesh)
plot.mesh.2.5D(mesh)

# set.seed(5847947)
# 
# a1 = rnorm(1,mean = 1, sd = 1)
# a2 = rnorm(1,mean = 1, sd = 1)
# a3 = rnorm(1,mean = 1, sd = 1)

# Evaluate exact solution on mesh nodes

# rho = sqrt(x^2 + y^2 + z^2)
# theta = 
# phi = acos(z/rho) = acos(z/sqrt(x^2 + y^2 + z^2))

set_theta <- function(x,y){
  result <- matrix(0, nrow=length(x), ncol=1)
  
  for(i in 1:length(x)){
    if(abs(x[i]) < 1e-13 & y[i] > 0){
      result[i] = pi/2
    }else if(abs(x[i]) < 1e-13 & y[i] < 0){
      result[i] = 3*pi/2
    }else if(x[i] > 0 & y[i] >= 0){
      result[i] = atan(y[i]/x[i])
    }else if( (x[i] > 0 & y[i] < 0) | (x[i] < 0 & y[i] > 0)){
      result[i] = atan(y[i]/x[i]) + 2*pi
    }else if((x[i] < 0 & y[i] <= 0)){
      result[i] = atan(y[i]/x[i]) + pi
    }else{
      result[i] = NA # x[i] == 0 & y[i] == 0
    }
  }
  return(result)
}

Cov1 = function(x, y, z){
  rho = sqrt(x^2 + y^2 + z^2)
  phi = acos(z/rho)
  theta = set_theta(x,y)
  
  return( (sin(theta)^2 * cos(2*phi) + cos(theta)^3) )
  
}

Cov1 = function(x, y, z){
  (1-z^2)*cos(2* ((x^2-y^2)/(x^2+y^2))) + z^3
}

Cov2 = function(x, y, z){
  1/4 * (3*z^2 - 1)
}

function_1 <- function(x,y,z){
  x^2- y^2 + z^2
}


function_2 <- function(x,y,z){
  sin(3*x)*cos(2*y) + cos(5*z)
}


func3 <- function(x,y,z){
  exp(x) * log(1+y^2+z^2)
}

plot.FEM.2.5D(FEM(coeff= Cov1(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis))
plot.colorbar(FEM(coeff= Cov1(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis), 
              file = "prova")
plot.FEM.2.5D(FEM(coeff= Cov2(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis))

plot.FEM.2.5D(FEM(coeff= function_1(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis))
plot.FEM.2.5D(FEM(coeff= function_2(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis)) # ci sta
plot.FEM.2.5D(FEM(coeff= func3(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis))



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
#n_obs <- c(500, 1000)
rmse <- function(x,y){
  return(sqrt( mean( (x-y)^2 ) ) ) 
}

# ------------------------------------------------------------------------------

lambda= 10^seq(-1,1,length=15)

# Building folders -------------------------------------------------------------
date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]

#date_ <- "2024-05-14-20_34_15"

if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists("data/simulation_1/")){
  dir.create("data/simulation_1/")
}

folder.name = paste("data/simulation_1/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

set.seed(0)
for(j in 1:length(n_obs)){
  idx <- sample(1:nnodes, size=n_obs[j])
  locations <- runif(n_obs, )
  nlocs = dim(locations)[1]
  
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2], locations[,3]),
              Cov2(locations[,1], locations[,2], locations[,3]))
  X2 = cbind( Cov1(locations[,1], locations[,2], locations[,3]),
              Cov2(locations[,1], locations[,2], locations[,3]))
  
  func1 = function_1(x=locations[,1], y=locations[,2],  z = locations[,3])
  func2 = function_2(x=locations[,1], y=locations[,2],  z = locations[,3])
  
  test.func1 = function_1(x=test.locations[,1], y=test.locations[,2],  z=test.locations[,3])
  test.func2 = function_2(x=test.locations[,1], y=test.locations[,2],  z=test.locations[,3])
  
  results <- list(
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=2), # (b11, b12)
    b_2 = matrix(0, nrow=n_sim,ncol=2), # (b21, b22)
    time = matrix(0, nrow=n_sim, ncol=1)
  )
  
  results <- list(monolitich=results, iterative=results)
  
  errors <- list(
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=2),
    b_2 = matrix(0, nrow=n_sim,ncol=2),
    f_1 = matrix(0, nrow=n_sim,ncol=1),
    f_2 = matrix(0, nrow=n_sim,ncol=1),
    response = matrix(0, nrow=n_sim,ncol=1)
  )
  errors = list(monolitich = errors, iterative = errors)
  
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
    start_ <- Sys.time()
    invisible(capture.output(output_iterative <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                            covariates = X, random_effect = c(1,2),
                                                                            FEMbasis = FEMbasis, lambda = lambda, 
                                                                            lambda.selection.criterion = "grid", 
                                                                            lambda.selection.lossfunction = "GCV",
                                                                            DOF.evaluation = "exact", FLAG_ITERATIVE = TRUE)))
    results$iterative$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    best_lambda <- output_iterative$bestlambda
    #results
    results$iterative$beta_1[i] <- output_iterative$beta[1, best_lambda]
    results$iterative$beta_2[i] <- output_iterative$beta[2, best_lambda]
    results$iterative$b_1[i,] <- output_iterative$b_i[1:2, best_lambda]
    results$iterative$b_2[i,] <- output_iterative$b_i[3:4, best_lambda]
    
    # errors
    errors$iterative$beta_1[i] <- rmse(output_iterative$beta[1, best_lambda], betas[1])
    errors$iterative$beta_2[i] <- rmse(output_iterative$beta[2, best_lambda], betas[2])
    errors$iterative$b_1[i,] <- c(rmse(results$fdaPDE$b_1[i,1], b_1[1]), rmse(results$fdaPDE$b_1[1,2], b_1[2]))
    errors$iterative$b_2[i,] <- c(rmse(results$fdaPDE$b_2[i,1], b_2[1]), rmse(results$fdaPDE$b_2[1,2], b_2[2]))
    errors$iterative$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_iterative$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations),
                                    test.func1)
    errors$iterative$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output_iterative$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations),
                                    test.func2)
    
    
    y_hat1 <- test.X1%*% as.matrix(output_iterative$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_iterative$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations) + 
      test.X1%*% as.matrix(output_iterative$b_i[1:2, best_lambda])
    
    y_hat2 <- test.X2%*% as.matrix(output_iterative$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_iterative$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations) + 
      test.X2%*%as.matrix(output_iterative$b_i[3:4, best_lambda])
    
    test.observations <- c(test.X1%*% betas + test.func1 + test.X1%*% b_1 ,
                           test.X2%*% betas + test.func2 + test.X2%*% b_2 )
    errors$iterative$response[i] <- rmse(c(y_hat1, y_hat2), test.observations)
    
    
    start_ <- Sys.time()
    invisible(capture.output(output_monolitich <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                               covariates = X, random_effect = c(1,2),
                                                                               FEMbasis = FEMbasis, lambda = lambda, 
                                                                               lambda.selection.criterion = "grid", 
                                                                               lambda.selection.lossfunction = "GCV",
                                                                               DOF.evaluation = "exact", FLAG_ITERATIVE = FALSE)))
    results$monolitich$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    best_lambda <- output_monolitich$bestlambda
    #results
    results$monolitich$beta_1[i] <- output_monolitich$beta[1, best_lambda]
    results$monolitich$beta_2[i] <- output_monolitich$beta[2, best_lambda]
    results$monolitich$b_1[i,] <- output_monolitich$b_i[1:2, best_lambda]
    results$monolitich$b_2[i,] <- output_monolitich$b_i[3:4, best_lambda]
    
    # errors
    errors$monolitich$beta_1[i] <- rmse(output_monolitich$beta[1, best_lambda], betas[1])
    errors$monolitich$beta_2[i] <- rmse(output_monolitich$beta[2, best_lambda], betas[2])
    errors$monolitich$b_1[i,] <- c(rmse(results$fdaPDE$b_1[i,1], b_1[1]), rmse(results$fdaPDE$b_1[1,2], b_1[2]))
    errors$monolitich$b_2[i,] <- c(rmse(results$fdaPDE$b_2[i,1], b_2[1]), rmse(results$fdaPDE$b_2[1,2], b_2[2]))
    errors$monolitich$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_monolitich$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations),
                                    test.func1)
    errors$monolitich$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output_monolitich$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations),
                                    test.func2)
    
    
    y_hat1 <- test.X1%*% as.matrix(output_monolitich$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_monolitich$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations) + 
      test.X1%*% as.matrix(output_monolitich$b_i[1:2, best_lambda])
    
    y_hat2 <- test.X2%*% as.matrix(output_monolitich$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_monolitich$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations) + 
      test.X2%*%as.matrix(output_monolitich$b_i[3:4, best_lambda])
    
    test.observations <- c(test.X1%*% betas + test.func1 + test.X1%*% b_1 ,
                           test.X2%*% betas + test.func2 + test.X2%*% b_2 )
    errors$monolitich$response[i] <- rmse(c(y_hat1, y_hat2), test.observations)
    
  }
  save(errors, results,
       file = paste0(folder.name, "n_obs_", n_obs[j],".RData"))
  
}

# post-proc --------------------------------------------------------------------
n_sim <- 30
n_obs <- c(100, 250, 500, 1000)
betas <- c(3, 0.5)

b_1 <- as.matrix( c(-5, 0.) )
b_2 <- as.matrix(c(5,0))

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
                        time = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        n_obs = as.factor(rep(n_obs, each=n_sim))
)

estimates2 <- estimates

date_ = "2024-05-14-20_34_15"
folder.name = paste("data/simulation_1/",date_,"/",sep="")

for(i in 1:length(n_obs)){
  load(file = paste0(folder.name, "n_obs_", n_obs[i],".RData"))
  
  estimates$beta1[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$beta_1
  estimates$beta2[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$beta_2
  estimates$beta1_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$monolitich$beta_1
  estimates$beta2_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$monolitich$beta_2
  estimates$f1_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$monolitich$f_1
  estimates$f2_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$monolitich$f_2
  estimates$response_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$monolitich$response
  estimates$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("monolithic", times=n_sim)
  estimates$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$time
  
  estimates2$beta1[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$beta_1
  estimates2$beta2[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$beta_2
  estimates2$beta1_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$iterative$beta_1
  estimates2$beta2_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$iterative$beta_2
  estimates2$f1_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$iterative$f_1
  estimates2$f2_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$iterative$f_2
  estimates2$f3_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$iterative$f_3
  estimates2$response_rmse[(1+(i-1)*n_sim):(i*n_sim)] <- errors$iterative$response
  estimates2$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("iterative", times=n_sim)
  estimates2$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$time
  
}


estimates <- rbind(estimates, estimates2)
estimates$method <- as.factor(estimates$method)
at_ <- c(1:2, 4:5, 7:8, 10:11)
fill_col <- viridis::viridis(3, begin=0.25, end=0.95)
fill_col <- fill_col[1:2]
legend <- levels(estimates$method)

{
  mar_ = par("mar")
  mar_[2] = mar_[2] + 0.25
  pdf(paste0(folder.name, "simulation_1.pdf"), family = "serif", width = 7, height = 7)
  par(mar=mar_)
  boxplot(estimates$beta1_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(beta[1]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$beta2_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(beta[2]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$beta1 ~estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(hat(beta)[1]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  abline(h=betas[1], lty=2, lwd=3, col="red")
  par(mar=mar_)
  boxplot(estimates$beta2 ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(hat(beta)[2]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  abline(h=betas[2], lty=2, lwd=3, col="red")
  par(mar=mar_)
  boxplot(estimates$f1_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(f[1]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$f2_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(f[2]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$time ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          ylim=c(55, 280),
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main ="execution time [s]")
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  
  # -- 8-14 no main
  
  par(mar=mar_)
  boxplot(estimates$beta1_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$beta2_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$beta1 ~estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  abline(h=betas[1], lty=2, lwd=3, col="red")
  par(mar=mar_)
  boxplot(estimates$beta2 ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  abline(h=betas[2], lty=2, lwd=3, col="red")
  par(mar=mar_)
  boxplot(estimates$f1_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$f2_rmse ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mar=mar_)
  boxplot(estimates$time ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="execution time [s]", xlab="observations", at = at_, xaxt="n", 
          ylim=c(55, 280),
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2)
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  
  
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

func1 = function_1(x=locations[,1], y=locations[,2],  z = locations[,3])
func2 = function_2(x=locations[,1], y=locations[,2],  z = locations[,3])

estimates_monolitich <- list(
  f_1 = matrix(0, nrow=nnodes,ncol=1),
  f_2 = matrix(0, nrow=nnodes,ncol=1)
)

estimates_iterative <- list(
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
  invisible(capture.output(output_iterative <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                          covariates = X, random_effect = c(1,2),
                                                                          FEMbasis = FEMbasis, lambda = lambda, 
                                                                          lambda.selection.criterion = "grid", 
                                                                          lambda.selection.lossfunction = "GCV",
                                                                          DOF.evaluation = "exact", FLAG_ITERATIVE = TRUE)))
  
  best_lambda <- output_iterative$bestlambda
  estimates_iterative$f_1 <- estimates_iterative$f_1 + output_iterative$fit.FEM.mixed$coeff[1:nnodes,best_lambda] / n_sim
  estimates_iterative$f_2 <- estimates_iterative$f_2 + output_iterative$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda] / n_sim
  
  invisible(capture.output(output_monolitich <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                             covariates = X, random_effect = c(1,2),
                                                                             FEMbasis = FEMbasis, lambda = lambda, 
                                                                             lambda.selection.criterion = "grid", 
                                                                             lambda.selection.lossfunction = "GCV",
                                                                             DOF.evaluation = "exact", FLAG_ITERATIVE = FALSE)))
  
  best_lambda <- output_monolitich$bestlambda
  estimates_monolitich$f_1 <- estimates_monolitich$f_1 + output_monolitich$fit.FEM.mixed$coeff[1:nnodes,best_lambda] / n_sim
  estimates_monolitich$f_2 <- estimates_monolitich$f_2 + output_monolitich$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda] / n_sim
}

save(estimates_monolitich, estimates_iterative,
     file=paste0(folder.name, "n_obs_250_estimates.RData"))

library(plotly)
options(warn=-1)

func1 = function_1(x=mesh$nodes[,1], y=mesh$nodes[,2],  z = mesh$nodes[,3])
func2 = function_2(x=mesh$nodes[,1], y=mesh$nodes[,2],  z = mesh$nodes[,3])
# FEMbasis.ref <- create.FEM.basis(mesh.ref)

{
  FEMobject <-  FEM(coeff = func1, FEMbasis)
  plot.FEM.2.5D(FEMobject, colorscale = viridis)
  snapshot3d(filename = paste0(folder.name,"true_f1.png"),
            fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
  plot.colorbar(FEMobject, colorscale =  viridis, 
              file = paste0(folder.name, "colorbar_f1"))
  
  plot.FEM.2.5D(FEM(coeff = estimates_iterative$f_1, FEMbasis), 
                colorscale = viridis, limits = compute_limits(FEMobject))
  snapshot3d(filename = paste0(folder.name,"iterative_f1.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
  plot.FEM.2.5D(FEM(coeff = estimates_monolitich$f_1, FEMbasis), 
                colorscale = viridis, limits = compute_limits(FEMobject))
  snapshot3d(filename = paste0(folder.name,"monolithic_f1.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()

}

{
  FEMobject <-  FEM(coeff = func2, FEMbasis)
  plot.FEM.2.5D(FEMobject, colorscale = viridis)
  snapshot3d(filename = paste0(folder.name,"true_f2.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
  plot.colorbar(FEMobject, colorscale =  viridis, 
                file = paste0(folder.name, "colorbar_f2"))
  
  plot.FEM.2.5D(FEM(coeff = estimates_iterative$f_2, FEMbasis), 
                colorscale = viridis, limits = compute_limits(FEMobject))
  snapshot3d(filename = paste0(folder.name,"iterative_f2.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
  plot.FEM.2.5D(FEM(coeff = estimates_monolitich$f_2, FEMbasis), 
                colorscale = viridis, limits = compute_limits(FEMobject))
  snapshot3d(filename = paste0(folder.name,"monolithic_f2.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
}

{
  FEMobject <- FEM(coeff= Cov1(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis)
  plot.FEM.2.5D(FEMobject, colorscale = viridis)
  snapshot3d(filename = paste0(folder.name,"cov_1.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
  plot.colorbar(FEMobject, colorscale =  viridis, 
                file = paste0(folder.name, "colorbar_cov1"))
}

{
  FEMobject <- FEM(coeff= Cov2(mesh$nodes[,1], mesh$nodes[,2], mesh$nodes[,3]), FEMbasis)
  plot.FEM.2.5D(FEMobject, colorscale = viridis)
  snapshot3d(filename = paste0(folder.name,"cov_2.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
  plot.colorbar(FEMobject, colorscale =  viridis, 
                file = paste0(folder.name, "colorbar_cov2"))
}

{
  plot.mesh.2.5D(mesh)
  snapshot3d(filename = paste0(folder.name,"unit_ball.png"),
             fmt = "png", width = 800, height = 750, webshot = rgl.useNULL())
  rgl.close()
  
}
