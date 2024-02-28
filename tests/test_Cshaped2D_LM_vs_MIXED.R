if(!require(fdaPDEmixed)){
  devtools::install_github(repo ="aldoclemente/fdaPDEmixed")
}
if(!require(rstudioapi)) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#-------------------------------------------------------------------------------
library(fdaPDEmixed)

data(horseshoe2D)
mesh=create.mesh.2D(nodes=horseshoe2D$boundary_nodes, 
                    segments = horseshoe2D$boundary_segments)
mesh = refine.mesh.2D(mesh, maximum_area = 0.025, minimum_angle = 30)
locations <- mesh$nodes[!mesh$nodesmarkers, ]

mesh = refine.mesh.2D(mesh, maximum_area = 0.0075, minimum_angle = 30)
FEMbasis <- create.FEM.basis(mesh)

test.locations <- refine.mesh.2D(mesh, maximum_area = 0.0025, minimum_angle = 30)$nodes
plot(mesh)
points(locations, pch=16, col="blue")
points(test.locations, pch=16, col="red")

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

nlocs <- nrow(locations)
nnodes <- nrow(mesh$nodes)
betas <- as.matrix(c(3,0.5))

X = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))

# image(FEM(Cov1(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis))
# image(FEM(rnorm(nnodes, sd=2), FEMbasis))
# 
# image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes)), FEMbasis))
# image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(2, nnodes)), FEMbasis))
# image(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes)), FEMbasis))

func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
# func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
# func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
test.func1 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(0.5, nrow(test.locations)))
# test.func2 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1, nrow(test.locations)))
# test.func3 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1.5, nrow(test.locations)))

n_sim <- 30
results_mixed_dummy <- list(
  beta_1 = matrix(0, nrow=n_sim,ncol=1),
  beta_2 = matrix(0, nrow=n_sim,ncol=1),
  time = matrix(0, nrow=n_sim,ncol=1)
)

results_lm <- results_mixed_dummy

errors_mixed_dummy <- list(
  beta_1 = matrix(0, nrow=n_sim,ncol=1),
  beta_2 = matrix(0, nrow=n_sim,ncol=1),
  f_1 = matrix(0, nrow=n_sim,ncol=1),
  response = matrix(0, nrow=n_sim,ncol=1)
)

errors_lm <- errors_mixed_dummy

results_mixed <- list(
  beta_1 = matrix(0, nrow=n_sim,ncol=1),
  beta_2 = matrix(0, nrow=n_sim,ncol=1),
  b_1 = matrix(0, nrow=n_sim,ncol=2), # (b11, b12)
  time = matrix(0, nrow=n_sim,ncol=1)
)

errors_mixed <- list(
  beta_1 = matrix(0, nrow=n_sim,ncol=1),
  beta_2 = matrix(0, nrow=n_sim,ncol=1),
  b_1 = matrix(0, nrow=n_sim, ncol=2),
  f_1 = matrix(0, nrow=n_sim,ncol=1),
  response = matrix(0, nrow=n_sim,ncol=1)
)

rmse <- function(x,y){
  return(sqrt( mean( x-y )^2 ) ) 
}

for(i in 1:n_sim){
  exact1 <- X%*% betas + func1 
  observations <- exact1
  observations <- observations + rnorm(nlocs, mean=0, sd=0.05*(diff(range(func1))))
  observations <- as.matrix(observations)
  
  #lambda= 10^seq(-2,1,by=0.1) # 31
  lambda= 10^seq(-1,1,length=10)
  start <- Sys.time()
  output_mixed_dummy = smooth.FEM.mixed(observations = observations, locations = locations,
                                        covariates = X,
                                        FEMbasis = FEMbasis, lambda = lambda, 
                                        lambda.selection.criterion = "grid", 
                                        lambda.selection.lossfunction = "GCV",
                                        DOF.evaluation = "exact")
  results_mixed_dummy$time[i] <- difftime(Sys.time(), start, units="secs")
  
  start <- Sys.time()
  output_mixed = smooth.FEM.mixed(observations = observations, locations = locations,
                                  covariates = X, random_effect = c(1,2),
                                  FEMbasis = FEMbasis, lambda = lambda, 
                                  lambda.selection.criterion = "grid", 
                                  lambda.selection.lossfunction = "GCV",
                                  DOF.evaluation = "exact")
  results_mixed$time[i] <- difftime(Sys.time(), start, units="secs")
  
  start <- Sys.time()
  output_lm = fdaPDE::smooth.FEM(observations = observations, locations = locations,
                                 covariates = X,
                                 FEMbasis = FEMbasis, lambda = lambda, 
                                 lambda.selection.criterion = "grid", 
                                 lambda.selection.lossfunction = "GCV",
                                 DOF.evaluation = "exact")
  results_lm$time[i] <- difftime(Sys.time(), start, units="secs")
  
  best_lambda_mixed_dummy <- output_mixed_dummy$bestlambda
  best_lambda_mixed <- output_mixed$bestlambda
  best_lambda_lm <- output_lm$optimization$lambda_position
  
  #results DUMMY mixed ---------------------------------------------------------
  results_mixed_dummy$beta_1[i] <- output_mixed_dummy$beta[1, best_lambda_mixed_dummy]
  results_mixed_dummy$beta_2[i] <- output_mixed_dummy$beta[2, best_lambda_mixed_dummy]
  
  #results mixed ---------------------------------------------------------------
  results_mixed$beta_1[i] <- output_mixed$beta[1, best_lambda_mixed]
  results_mixed$beta_2[i] <- output_mixed$beta[2, best_lambda_mixed]
  results_mixed$b_1[i,] <- output_mixed$b_i[, best_lambda_mixed]
  
  #results lm ------------------------------------------------------------------
  results_lm$beta_1[i] <- output_lm$solution$beta[1]
  results_lm$beta_2[i] <- output_lm$solution$beta[2]
  
  # errors DUMMY mixed ---------------------------------------------------------
  errors_mixed_dummy$beta_1[i] <- rmse(output_mixed_dummy$beta[1, best_lambda_mixed_dummy], betas[1])
  errors_mixed_dummy$beta_2[i] <- rmse(output_mixed_dummy$beta[2, best_lambda_mixed_dummy], betas[2])
  
  errors_mixed_dummy$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_mixed_dummy$fit.FEM.mixed$coeff[1:nnodes,best_lambda_mixed_dummy]), FEMbasis), test.locations),
                                    test.func1)
  
  y_hat1 <- X%*% as.matrix(output_mixed_dummy$beta[, best_lambda_mixed_dummy]) + 
    eval.FEM(FEM(as.matrix(output_mixed_dummy$fit.FEM.mixed$coeff[1:nnodes,best_lambda_mixed_dummy]), FEMbasis), locations)
  
  errors_mixed_dummy$response[i] <- rmse(y_hat1, observations)
  
  # errors DUMMY mixed ---------------------------------------------------------
  errors_mixed$beta_1[i] <- rmse(output_mixed$beta[1, best_lambda_mixed], betas[1])
  errors_mixed$beta_2[i] <- rmse(output_mixed$beta[2, best_lambda_mixed], betas[2])
  
  errors_mixed$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_mixed$fit.FEM.mixed$coeff[1:nnodes,best_lambda_mixed]), FEMbasis), test.locations),
                              test.func1)
  
  y_hat1 <- X%*% as.matrix(output_mixed$beta[, best_lambda_mixed]) + 
    eval.FEM(FEM(as.matrix(output_mixed$fit.FEM.mixed$coeff[1:nnodes,best_lambda_mixed]), FEMbasis), locations) +
    X%*% as.matrix(output_mixed$b_i[, best_lambda_mixed])
  
  errors_mixed$response[i] <- rmse(y_hat1, observations)
  
  # errors lm ------------------------------------------------------------------
  errors_lm$beta_1[i] <- rmse(output_lm$solution$beta[1], betas[1])
  errors_lm$beta_2[i] <- rmse(output_lm$solution$beta[2], betas[2])
  
  errors_lm$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output_lm$fit.FEM$coeff[1:nnodes,]), FEMbasis), test.locations),
                           test.func1)
  
  y_hat1 <- X%*% as.matrix(output_lm$solution$beta) + 
    eval.FEM(FEM(as.matrix(output_lm$fit.FEM$coeff[1:nnodes,]), FEMbasis), locations)
  
  errors_lm$response[i] <- rmse(y_hat1, observations)
}


errors_dummy <- results_mixed_dummy
errors_dummy$beta_1 <- errors_dummy$beta_1 - results_lm$beta_1
errors_dummy$beta_2 <- errors_dummy$beta_2 - results_lm$beta_2


errors <- results_mixed
errors$beta_1 <- errors$beta_1 - results_lm$beta_1
errors$beta_2 <- errors$beta_2 - results_lm$beta_2

if(!dir.exists("data/")) {
  dir.create("data/")
}

folder.name = "data/test_Cshaped2D/"

if( !dir.exists(folder.name)){
  dir.create(folder.name)
}

# store results
save(errors_lm , errors_mixed_dummy, errors_mixed, 
     results_lm, results_mixed_dummy, results_mixed , 
     errors_dummy, errors,
     file = paste0(folder.name, "data.RData"))

# imgs -------------------------------------------------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
folder.name = "data/test_Cshaped2D/"
if(!dir.exists(folder.name)) cat("Warning! Data not available\n")
  
pdf(paste0(folder.name, "results.pdf"))
boxplot(data.frame(beta_1 = errors_dummy$beta_1, beta_2 = errors_dummy$beta_2), main="Errore \n(Rand. Effect ==NULL)")
abline(h=0., lty=2, col="red")

boxplot(data.frame(beta_1 = errors$beta_1, beta_2 = errors$beta_2), main="Errore \n (Rand. Effect !=NULL)")
abline(h=0., lty=2, col="red")

boxplot(data.frame(f_1_mixed_dummy = errors_mixed_dummy$f_1,
                   f_1_mixed = errors_mixed$f_1,
                   f_1_lm = errors_lm$f_1),
        main="RMSE\nf")

boxplot(data.frame(beta_1_mixed_dummy = results_mixed_dummy$beta_1, 
                   beta_1_mixed = results_mixed$beta_1,
                   beta_1_lm = results_lm$beta_1), main =" Beta 1")
abline(h=betas[1], lty=2, col="red")

boxplot(data.frame(beta_2_mixed_dummy = results_mixed_dummy$beta_2, 
                   beta_2_mixed = results_mixed$beta_2,
                   beta_2_lm = results_lm$beta_2), main =" Beta 2")
abline(h=betas[2], lty=2, col="red")

boxplot(data.frame(b_11 = results_mixed$b_1[,1], b_12 = results_mixed$b_1[,2]), main="b \n(Rand. Effect !=NULL)")
abline(h=0., lty=2, col="red")

boxplot(data.frame(mixed_dummy=results_mixed_dummy$time,
                   mixed = results_mixed$time,
                   lm = results_lm$time), main="Time [s]")

dev.off()


source("utils.R")

png(paste0(folder.name,"Cshaped_mesh.png"))
plot(mesh, pch=".", asp=1)
dev.off()

png(paste0(folder.name,"Cshaped_locs.png"))
plot(mesh, pch=".", asp=1)
points(locations, pch=16, col="red3")
dev.off()

if(!require(plotly)) install.packages("plotly")

fig <- contour.FEM(FEM(Cov1(mesh$nodes[,1], mesh$nodes[,2]), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))

fig <- contour.FEM(FEM(rnorm(nnodes, sd=2), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))


fig <-contour.FEM(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(0.5, nnodes)), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))

fig <- contour.FEM(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(2, nnodes)), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))


fig <-contour.FEM(FEM(fs.test.time(x=mesh$nodes[,1], y=mesh$nodes[,2],  t=rep(1.5, nnodes)), FEMbasis))
fig %>% layout( scene = list(
  camera = list(
    eye = list(x = 0, y = -0.01,  z = 2.25))))
