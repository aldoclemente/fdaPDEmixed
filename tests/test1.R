# Test 1 -----------------------------------------------------------------------
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

data(horseshoe2D)
mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

# #### uniform Grid Initialization
# lambda= 10^seq(-2,1,by=0.1)
# output_CPP1 = fdaPDE::smooth.FEM(
#   observations = as.matrix(sin(pi*mesh$nodes[,1])), #observations at mesh nodes
#   FEMbasis = FEMbasis,
#   lambda = lambda,
#   lambda.selection.criterion = "grid", 
#   lambda.selection.lossfunction = "GCV",
#   DOF.evaluation = "exact")
# 
# gridnum = 20
# minV1 = min(mesh$nodes[,1])
# maxV1 = max(mesh$nodes[,1])
# minV2 = min(mesh$nodes[,2])
# maxV2 = max(mesh$nodes[,2])
# 
# x = seq(minV1+0.125, maxV1-0.125, length.out = gridnum) ### SHOULD ADJUST THIS! (OR ELSE, ERROR IN SOAP FILM)
# y= seq(minV2+0.125, maxV2-0.125, length.out = gridnum)
# unif_grid <- expand.grid(x = x, y = y)
# 
# points1=fdaPDE::eval.FEM(output_CPP1$fit.FEM,
#                  locations = unif_grid)
# 
# loc = unif_grid[-which(is.na(points1)),]
# 
# # randomize
# set.seed(5847947)
# ind = sample.int(dim(loc)[1])[1:100] ##100 locations points as default

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
b <- as.matrix( c(-5, 0, 5) )

test.X1 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
test.X2 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
test.X3 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))

n_sim <- 30
n_obs <- c(100, 250, 500, 1000)

rmse <- function(x,y){
  return(sqrt( mean( x-y )^2 ) ) 
}

# ------------------------------------------------------------------------------

# mgcv boundary
bnd = mesh$nodes[as.logical(mesh$nodesmarkers), ]
bound <- list(list(x = bnd[,1], y = bnd[,2])) #, f = rep(0, nrow(crds))))
names(bound[[1]]) <- c("x", "y") #, "f")

# mgcv knots 
knots <- data.frame(v=rep(seq(-.5,3,by=.5),4),
                    w=rep(c(-.6,-.3,.3,.6),rep(8,4)))
names(knots) <- c("x", "y")

lambda= 10^seq(-2,1,by=0.1)

# Building folders -------------------------------------------------------------
date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]
if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists("data/test_1/")){
  dir.create("data/test_1/")
}

folder.name = paste("data/test_1/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}


for(j in 1:length(n_obs)){
  locations <- sample_locations(mesh, n_obs[j])
  nlocs = dim(locations)[1]
  
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  
  func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
  test.func1 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(0.5, nrow(test.locations)))
  test.func2 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1, nrow(test.locations)))
  test.func3 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1.5, nrow(test.locations)))
  
  
  results <- list(
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=2), # (b11, b12)
    b_2 = matrix(0, nrow=n_sim,ncol=2), # (b21, b22)
    b_3 = matrix(0, nrow=n_sim,ncol=2)  # (b31, b32)
  )
  
  results <- list(fdaPDE=results, mgcv=results)
  
  errors <- list(
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=2),
    b_2 = matrix(0, nrow=n_sim,ncol=2),
    b_3 = matrix(0, nrow=n_sim,ncol=2),
    f_1 = matrix(0, nrow=n_sim,ncol=1),
    f_2 = matrix(0, nrow=n_sim,ncol=1),
    f_3 = matrix(0, nrow=n_sim,ncol=1),
    response = matrix(0, nrow=n_sim,ncol=1)
  )
  errors = list(fdaPDE = errors, mgcv = errors)
  
  for(i in 1:n_sim){
    # V == Cov1
    exact1 <- X1%*% betas + func1 + X1[,1] * b[1] 
    obs1 <- exact1  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact1)))
    
    exact2 <- X2%*% betas + func2 + X2[,1] * b[2]
    obs2 <- exact2  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact2)))
    
    exact3 <- X3%*% betas + func3 + X3[,1] * b[3] 
    obs3 <- exact3  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact3)))
    
    exact <- c(exact1, exact2, exact3)
    observations <- c(obs1, obs2, obs3)
    observations <- observations + rnorm(nlocs*3, mean=0, sd=0.05*(diff(range(c(func1, func2, func3)))))
    observations <- matrix(observations, nrow=nlocs, ncol=3)
    X = rbind(X1, X2, X3)
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
    results$fdaPDE$b_3[i,] <- output_fdaPDE$b_i[5:6, best_lambda]
    
    # errors
    errors$fdaPDE$beta_1[i] <- rmse(output_fdaPDE$beta[1, best_lambda], betas[1])
    errors$fdaPDE$beta_2[i] <- rmse(output_fdaPDE$beta[2, best_lambda], betas[2])
    errors$fdaPDE$b_1[i,] <- c(rmse(results$fdaPDE$b_1[i,1], b[1]), rmse(results$fdaPDE$b_1[1,2], 0.))
    errors$fdaPDE$b_2[i,] <- c(rmse(results$fdaPDE$b_2[i,1], b[2]), rmse(results$fdaPDE$b_2[1,2], 0.))
    errors$fdaPDE$b_3[i,] <- c(rmse(results$fdaPDE$b_3[i,1], b[3]), rmse(results$fdaPDE$b_3[1,2], 0.))
    
    errors$fdaPDE$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations),
                          test.func1)
    errors$fdaPDE$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations),
                          test.func2)
    errors$fdaPDE$f_3[i] <- rmse(eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda]), FEMbasis), test.locations),
                          test.func3)
    
    y_hat1 <- test.X1%*% as.matrix(output_fdaPDE$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), FEMbasis), test.locations) + 
      test.X1%*% as.matrix(output_fdaPDE$b_i[1:2, best_lambda])
    
    y_hat2 <- test.X2%*% as.matrix(output_fdaPDE$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]), FEMbasis), test.locations) + 
      test.X2%*%as.matrix(output_fdaPDE$b_i[3:4, best_lambda])
    
    y_hat3 <- test.X3%*% as.matrix(output_fdaPDE$beta[, best_lambda]) + 
      eval.FEM(FEM(as.matrix(output_fdaPDE$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda]), FEMbasis), test.locations) + 
      test.X3%*% as.matrix(output_fdaPDE$b_i[5:6, best_lambda])
    
    test.observations <- c(test.X1%*% betas + test.func1 + test.X1[,1] * b[1],
                           test.X2%*% betas + test.func2 + test.X2[,1] * b[2],
                           test.X3%*% betas + test.func3 + test.X3[,1] * b[3])
    errors$fdaPDE$response[i] <- rmse(c(y_hat1, y_hat2, y_hat3), test.observations)
    
    # mgcv (soap film smoother) --------------------------------------------------
    
    # mgcv_loc = rbind(locations, locations, locations)
    # mgcv_X1 = c(X1[,1], X2[,1], X3[,1])
    # mgcv_X2 = c(X1[,2], X2[,2], X3[,2])
    # #mgcv_obs = observations
    # mgcv_fac <- as.factor(rep(c(1,2,3),c(nlocs,nlocs,nlocs)))
    # 
    # mgcv_data <- data.frame(x = mgcv_loc[,1], 
    #                   y = mgcv_loc[,2], 
    #                   obs = as.vector(observations), 
    #                   cov1 = mgcv_X1, 
    #                   cov2 = mgcv_X2, 
    #                   fac=mgcv_fac)
    # 
    # output_mgcv <- mgcv::gam(obs ~ cov1+ cov2 +
    #                            s(cov1, fac, bs='re')+ #, bs='re')+ #allow penalization on random-effects
    #                            s(cov2, fac, bs='re')+ #, bs='re') + #s(fac, bs="re")+
    #                            s(x, y, by=fac, bs = "sf", # soap-film smoother #, id=1,
    #                              xt = list(bnd = bound))-1, 
    #                          data = mgcv_data, knots = knots, method="REML") 
  
  
    # output_mgcv <- gamm(obs ~ cov1+ cov2 +
    #                      s(x, y, by=fac, bs = "sf", # soap-film smoother #, id=1, 
    #                        xt = list(bnd = bound))-1, random=list(cov1=~1, cov2=~1),
    #                    data = mgcv_data)
    # output_mgcv <- output_mgcv$gam
    # 
    # results$mgcv$beta_1[i] <- output_mgcv$coefficients[1]
    # results$mgcv$beta_2[i] <- output_mgcv$coefficients[2]
    # 
    # errors$mgcv$beta_1[i] <- rmse(results$mgcv$beta_1[i], betas[1])
    # errors$mgcv$beta_2[i] <- rmse(results$mgcv$beta_2[i], betas[2])
    # errors$mgcv$response[i] <- rmse(as.vector(exact), output_mgcv$fitted.values)
    
    # fittare su test? boh 
    #errors$mgcv$response[i] <- rmse(  )
    
    # mgcv_loc = rbind(test.locations, test.locations, test.locations)
    # mgcv_X1 = c(test.X1[,1], test.X2[,1], test.X3[,1])
    # mgcv_X2 = c(test.X1[,2], test.X2[,2], test.X3[,2])
    # mgcv_fac <- as.factor(rep(c(1,2,3),c(nrow(test.locations),nrow(test.locations),nrow(test.locations))))
    # mgcv_nlocs = nrow(test.locations)*3
    # 
    # mgcv_test <- data.frame(x = mgcv_loc[,1],
    #                         y = mgcv_loc[,2],
    #                         obs = test.observations,
    #                         cov1 = mgcv_X1,
    #                         cov2 = mgcv_X2,
    #                         fac=mgcv_fac)
    # 
    # tmp <- predict.gam(output_mgcv,  newdata = mgcv_test[1:nrow(test.locations),],
    #                    type="terms", terms="s(x,y):fac1")
    # idx_na <- which(is.na(tmp))
    # if(length(idx_na) != 0){
    # errors$mgcv$f_1[i] <-  rmse(test.func1[-idx_na], tmp[-idx_na])
    # }else{
    #   errors$mgcv$f_1[i] <- rmse(test.func1, tmp)
    # }
    # tmp <- predict.gam(output_mgcv, newdata = mgcv_test[((nrow(test.locations)+1) : (2*nrow(test.locations))), ],
    #                    type="terms", terms="s(x,y):fac2")
    # idx_na <- which(is.na(tmp))
    # if(length(idx_na) != 0){
    #   errors$mgcv$f_2[i] <-  rmse(test.func2[-idx_na], tmp[-idx_na])
    # }else{
    #   errors$mgcv$f_2[i] <-  rmse(test.func2, tmp)
    # }
    # 
    # tmp <- predict.gam(output_mgcv,  newdata = mgcv_test[((2*nrow(test.locations)+1) : (3*nrow(test.locations))), ], 
    #                    type="terms", terms="s(x,y):fac3")
    # if(length(idx_na) != 0){
    #   errors$mgcv$f_3[i] <-  rmse(test.func3[-idx_na], tmp[-idx_na])
    # }else{
    #   errors$mgcv$f_3[i] <-  rmse(test.func3, tmp)
    # }
    # 
    # tmp <- predict.gam(output_mgcv, newdata =  mgcv_test, type = "response")
    # idx_na <- which(is.na(tmp))
    # if(length(idx_na) != 0){
    #   errors$mgcv$response[i] <- rmse(test.observations[-idx_na], tmp[-idx_na])
    # }else{
    #   errors$mgcv$response[i] <- rmse(test.observations, tmp)
    # }
  }
  save(errors, results,
       file = paste0(folder.name, "n_obs_", n_obs[j],".RData"))
  
}

# post-proc --------------------------------------------------------------------
n_sim <- 30
n_obs <- c(100, 250, 500, 1000)
betas <- as.matrix(c(3,0.5))
b <- as.matrix( c(-5, 0, 5) )
estimates <- data.frame(beta1 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        beta2 = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        beta1_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        beta2_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        f1_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        f2_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        f3_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        response_rmse = matrix(nrow = n_sim * length(n_obs), ncol = 1),
                        n_obs = as.factor(rep(n_obs, each=n_sim))
                       )

date_ = "2024-03-06-14_47_06"
folder.name = paste("data/test_1/",date_,"/",sep="")
{
load(file = paste0("data/test_1/", date_, "n_obs_100.RData"))
estimates$beta1[1:n_sim] <- results$fdaPDE$beta_1
estimates$beta2[1:n_sim] <- results$fdaPDE$beta_2
estimates$beta1_rmse[1:n_sim] <- errors$fdaPDE$beta_1
estimates$beta2_rmse[1:n_sim] <- errors$fdaPDE$beta_2
estimates$f1_rmse[1:n_sim] <- errors$fdaPDE$f_1
estimates$f2_rmse[1:n_sim] <- errors$fdaPDE$f_2
estimates$f3_rmse[1:n_sim] <- errors$fdaPDE$f_3
estimates$response_rmse[1:n_sim] <- errors$fdaPDE$response

load(file = paste0("data/test_1/", date_, "n_obs_250.RData"))
estimates$beta1[(n_sim+1):(2*n_sim)] <- results$fdaPDE$beta_1
estimates$beta2[(n_sim+1):(2*n_sim)] <- results$fdaPDE$beta_2
estimates$beta1_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$beta_1
estimates$beta2_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$beta_2
estimates$f1_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$f_1
estimates$f2_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$f_2
estimates$f3_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$f_3
estimates$response_rmse[(n_sim+1):(2*n_sim)] <- errors$fdaPDE$response

load(file = paste0("data/test_1/", date_, "n_obs_500.RData"))
estimates$beta1[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$beta_1
estimates$beta2[(2*n_sim+1):(3*n_sim)] <- results$fdaPDE$beta_2
estimates$beta1_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$beta_1
estimates$beta2_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$beta_2
estimates$f1_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$f_1
estimates$f2_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$f_2
estimates$f3_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$f_3
estimates$response_rmse[(2*n_sim+1):(3*n_sim)] <- errors$fdaPDE$response

load(file = paste0("data/test_1/", date_, "n_obs_1000.RData"))
estimates$beta1[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$beta_1
estimates$beta2[(3*n_sim+1):(4*n_sim)] <- results$fdaPDE$beta_2
estimates$beta1_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$beta_1
estimates$beta2_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$beta_2
estimates$f1_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$f_1
estimates$f2_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$f_2
estimates$f3_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$f_3
estimates$response_rmse[(3*n_sim+1):(4*n_sim)] <- errors$fdaPDE$response
}

{
mai_ = par("mai")
mai_[2] = mai_[2] + 0.075
pdf(paste0(folder.name, "test_1.pdf"), family = "serif", width = 7, height = 7)
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
par(mai=mai_)
boxplot(estimates$f3_rmse ~ estimates$n_obs,
        ylab="RMSE", xlab="observations",
        main =expression(f[3]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
dev.off()
}
# ------------------------------------------------------------------------------

boxplot(data.frame(fdaPDE=errors$fdaPDE$f_2,
                   mgcv=errors$mgcv$f_2), 
        ylab="RMSE", xlab="observations",
        main =expression(f[2]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
boxplot(data.frame(fdaPDE=errors$fdaPDE$f_3,
                   mgcv=errors$mgcv$f_3), 
        ylab="RMSE", xlab="observations",
        main =expression(f[3]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)

boxplot(data.frame(fdaPDE=errors$fdaPDE$beta_1,
                   mgcv=errors$mgcv$beta_1), 
        ylab="RMSE", xlab="observations",
        main =expression(beta[1]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)

boxplot(data.frame(fdaPDE=errors$fdaPDE$beta_2,
                   mgcv=errors$mgcv$beta_2), 
        ylab="RMSE",
        main =expression(beta[2]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)

boxplot(data.frame(fdaPDE=results$fdaPDE$beta_2,
                   mgcv=results$mgcv$beta_2),
        main =expression(hat(beta)[2]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
abline(h=betas[2], lty=2, lwd=3, col="red")

boxplot(data.frame(fdaPDE=errors$fdaPDE$f_1,
                   mgcv=errors$mgcv$f_1), 
        ylab="RMSE",
        main =expression(f[1]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
boxplot(data.frame(fdaPDE=errors$fdaPDE$f_2,
                   mgcv=errors$mgcv$f_2), 
        ylab="RMSE",
        main =expression(f[2]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)
boxplot(data.frame(fdaPDE=errors$fdaPDE$f_3,
                   mgcv=errors$mgcv$f_3), 
        ylab="RMSE",
        main =expression(f[3]), col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)

boxplot(data.frame(fdaPDE=errors$fdaPDE$response,
                   mgcv=errors$mgcv$response), 
        ylab="RMSE",
        main = "fitted", col="white",cex.lab = 2, cex.axis = 2, cex.main = 2)


# plot mean estimates n = 250  -------------------------------------------------

date_ = "2024-03-06-14_47_06"
folder.name = paste("data/test_1/",date_,"/",sep="")
j = 2
nnodes <- nrow(mesh$nodes)
locations <- sample_locations(mesh, n_obs[j])
nlocs <- nrow(locations)
X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))

func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
test.func1 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(0.5, nrow(test.locations)))
test.func2 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1, nrow(test.locations)))
test.func3 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1.5, nrow(test.locations)))

estimates <- list(
  f_1 = matrix(0, nrow=nnodes,ncol=1),
  f_2 = matrix(0, nrow=nnodes,ncol=1),
  f_3 = matrix(0, nrow=nnodes,ncol=1)
)
lambda= 10^seq(-2,1,by=0.1)
for(i in 1:n_sim){
  # V == Cov1
  exact1 <- X1%*% betas + func1 + X1[,1] * b[1] 
  obs1 <- exact1  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact1)))
  
  exact2 <- X2%*% betas + func2 + X2[,1] * b[2]
  obs2 <- exact2  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact2)))
  
  exact3 <- X3%*% betas + func3 + X3[,1] * b[3] 
  obs3 <- exact3  # + rnorm(nlocs,mean=0,sd=0.05*diff(range(exact3)))
  
  exact <- c(exact1, exact2, exact3)
  observations <- c(obs1, obs2, obs3)
  observations <- observations + rnorm(nlocs*3, mean=0, sd=0.05*(diff(range(c(func1, func2, func3)))))
  observations <- matrix(observations, nrow=nlocs, ncol=3)
  X = rbind(X1, X2, X3)
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
  estimates$f_3 <- estimates$f_3 + output_fdaPDE$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda] / n_sim
}

save(estimates, file=paste0(folder.name, "n_obs_250_estimates.RData"))
     
library(plotly)
options(warn=-1)
mesh.ref = refine.mesh.2D(mesh, minimum_angle = 30, maximum_area = 0.0025)
nnodes <- nrow(mesh.ref$nodes)
func1 = fs.test.time(x=mesh.ref$nodes[,1], y=mesh.ref$nodes[,2],  t=rep(0.5, nnodes))
func2 = fs.test.time(x=mesh.ref$nodes[,1], y=mesh.ref$nodes[,2],  t=rep(1, nnodes))
func3 = fs.test.time(x=mesh.ref$nodes[,1], y=mesh.ref$nodes[,2],  t=rep(1.5, nnodes))
FEMbasis.ref <- create.FEM.basis(mesh.ref)
fig1 <- contour.FEM(FEM(coeff = func1, FEMbasis.ref))
fig2 <- contour.FEM(FEM(coeff = func2, FEMbasis.ref))
fig3 <- contour.FEM(FEM(coeff = func3, FEMbasis.ref))
fig1 <- fig1 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
fig2 <- fig2 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
fig3 <- fig3 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
save_image(fig1, paste0(folder.name,"true_f1.pdf"))
save_image(fig2, paste0(folder.name,"true_f2.pdf"))
save_image(fig3, paste0(folder.name,"true_f3.pdf"))

estimates_f1 <- eval.FEM(FEM(coeff = estimates$f_1, FEMbasis), locations = mesh.ref$nodes)
estimates_f2 <- eval.FEM(FEM(coeff = estimates$f_2, FEMbasis), locations = mesh.ref$nodes)
estimates_f3 <- eval.FEM(FEM(coeff = estimates$f_3, FEMbasis), locations = mesh.ref$nodes)

fig1 <- contour.FEM(FEM(coeff = estimates_f1, FEMbasis.ref), 
                    limits = compute_limits(FEM(coeff = func1, FEMbasis.ref)))
fig2 <- contour.FEM(FEM(coeff = estimates_f2, FEMbasis.ref), 
                    limits = compute_limits(FEM(coeff = func2, FEMbasis.ref)))
fig3 <- contour.FEM(FEM(coeff = estimates_f3, FEMbasis.ref),
                    limits = compute_limits(FEM(coeff = func3, FEMbasis.ref)))
fig1 <- fig1 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
fig2 <- fig2 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
fig3 <- fig3 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
save_image(fig1, paste0(folder.name,"estimate_f1.pdf"))
save_image(fig2, paste0(folder.name,"estimate_f2.pdf"))
save_image(fig3, paste0(folder.name,"estimate_f3.pdf"))

fig_cov1 <- contour.FEM(FEM(coeff= Cov1(mesh.ref$nodes[,1], mesh.ref$nodes[,2]), FEMbasis.ref))
fig_cov1 <- fig_cov1 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))

fig_cov2 <- contour.FEM(FEM(coeff= rnorm(nrow(mesh.ref$nodes), mean = 0, sd = 2), FEMbasis.ref))
fig_cov2 <- fig_cov2 %>% layout(scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2))))
save_image(fig_cov1, paste0(folder.name,"cov_1.pdf"))
save_image(fig_cov2, paste0(folder.name,"cov_2.pdf"))

fig_mesh <- plot(mesh)

save_image(fig_mesh, paste0(folder.name, "Cshaped.pdf"))
