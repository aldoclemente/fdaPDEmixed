# Test 1 -----------------------------------------------------------------------
if(!require(fdaPDEmixed)){
  devtools::install_github(repo ="aldoclemente/fdaPDEmixed")
}

if(!require(mgcv)){
  install.packages("mgcv")
}

if(.Platform$GUI == "RStudio"){
  if(!require(rstudioapi)) install.packages("rstudioapi")
}
# ------------------------------------------------------------------------------
rm(list=ls())

library(fdaPDEmixed)
library(mgcv)
if(.Platform$GUI == "RStudio"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
source("utils.R")
source("../utils/utils.R")
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

test.X1 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
test.X2 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
test.X3 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))

n_sim <- 30
n_obs <- c(100, 250, 500, 1000)

rmse <- function(x,y){
  return(sqrt(mean( (x-y)^2 ) ) ) 
}

# ------------------------------------------------------------------------------
# nodi 1 , 55
# mgcv boundary
bnd = mesh$nodes[as.logical(mesh$nodesmarkers), ]

outer <- matrix(nrow=nrow(bnd), ncol=2)
R <- 0.05

for(i in 1:nrow(bnd)){
  if(i!=1 & i!=nrow(bnd)){
    slope <- ifelse( abs(bnd[i-1,1] - bnd[i+1,1]) < 1e-10, NA,
                     (bnd[i-1,2] - bnd[i+1,2]) / (bnd[i-1,1] - bnd[i+1,1]))
  }else if(i==1){
    slope <- ifelse( abs(bnd[nrow(bnd),1] - bnd[i+1,1]) < 1e-10, NA,
                     (bnd[nrow(bnd),2] - bnd[i+1,2]) / (bnd[nrow(bnd),1] - bnd[i+1,1]))
  }else if(i==nrow(bnd)){
    slope <- ifelse( abs(bnd[i-1,1] - bnd[1,1]) < 1e-10, NA,
                     (bnd[i-1,2] - bnd[1,2]) / (bnd[i-1,1] - bnd[1,1]))
  }
  slope <- ifelse(is.na(slope), 0 , -1/slope)
  if(!is.infinite(slope)){
    sinTheta <- sin(atan(slope))
    cosTheta <- cos(atan(slope))
  }else{
    sinTheta <- 1
    cosTheta <- 0
  }
  
  pt1 <- t(as.matrix(c(bnd[i,1] + cosTheta*R, 
                       bnd[i,2] + sinTheta*R)))
  pt2 <- t(as.matrix(c(bnd[i,1] - cosTheta*R, 
                       bnd[i,2] - sinTheta*R)))
  pts <- rbind(pt1, pt2)
  idx <- fdaPDE:::CPP_search_points(mesh, pts)
  outer[i,] <- pts[which(idx == 0),]
}

bound <- list(list(x = outer[,1], y = outer[,2])) #, f = rep(0, nrow(crds))))
names(bound[[1]]) <- c("x", "y") #, "f")

x11()
fdaPDE:::plot.mesh.2D(mesh)
points(outer, pch=16, col="forestgreen")
points(bound[[1]]$x, bound[[1]]$y, col="red", pch=1)

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

# method  | n_obs | hat{beta_1} | err beta_1 | hat{beta_2} | err beta_2 |
# alpha_1 | err alpha_1 | alpha_2 | err alpha_2 | alpha_3 | err alpha_3 |
# rmse f_1 | rmse f_2 | rmse f_3
results = data.frame(matrix(NA,nrow=0, ncol=15))
names(results) <- c("method", "n_obs", 
                    "beta_1", "err_beta_1",
                    "beta_2", "err_beta_2",
                    "alpha_1", "err_alpha_1",
                    "alpha_2", "err_alpha_2",
                    "alpha_3", "err_alpha_3",
                    "err_f_1", "err_f_2", "err_f_3")


#estimates <- list(rep(data.frame(matrix(NA, nrow=nnodes, ncol=n_sim)), times=length(n_obs)))
estimates <- list(fdaPDE=list(), mgcv=list())
#data.frame(matrix(NA, nrow=nnodes, ncol=nsim))

mask_smooth <- function(data, k, num_nodes=nnodes){ # k stat idx
  if(is.vector(data)) return(data[((k-1)*num_nodes+1):(k*num_nodes)])
    return(data[((k-1)*num_nodes+1):(k*num_nodes),])
}

set.seed(0)
for(j in 1:length(n_obs)){
  locations <- sample_locations.mgcv(mesh, n_obs[j])
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
  s
  
  estimates$fdaPDE[[as.character(n_obs[j])]] = matrix(NA, nrow=3*nnodes,ncol=n_sim)
  estimates$mgcv[[as.character(n_obs[j])]] = matrix(NA, nrow=3*nnodes,ncol=n_sim)
  
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
    lambda.selection.criterion = "grid" 
    lambda.selection.lossfunction = "GCV"
    DOF.evaluation = "exact"
    # fdaPDE ---------------------------------------------------------------------
    invisible(capture.output(output_fdaPDE <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                            covariates = X, random_effect = c(1),
                                                                            FEMbasis = FEMbasis, lambda = lambda,
                                                                            lambda.selection.criterion = lambda.selection.criterion,
                                                                            lambda.selection.lossfunction = lambda.selection.lossfunction,
                                                                            DOF.evaluation = DOF.evaluation, FLAG_ITERATIVE = TRUE)))
    #source("lalalal.R")
    best_lambda <- output_fdaPDE$bestlambda
    #results
    results <- rbind(results,
                      data.frame(method="fdaPDE", n_obs = n_obs[j], 
                      beta_1=output_fdaPDE$beta[1, best_lambda], err_beta_1=rmse(output_fdaPDE$beta[1, best_lambda], betas[1]),
                      beta_2=output_fdaPDE$beta[2, best_lambda], err_beta_2=rmse(output_fdaPDE$beta[2, best_lambda], betas[2]),
                      alpha_1=output_fdaPDE$b_i[1, best_lambda], err_alpha_1=rmse(output_fdaPDE$b_i[1, best_lambda], b[1]),
                      alpha_2=output_fdaPDE$b_i[2, best_lambda], err_alpha_2=rmse(output_fdaPDE$b_i[2, best_lambda], b[2]),
                      alpha_3=output_fdaPDE$b_i[3, best_lambda], err_alpha_3=rmse(output_fdaPDE$b_i[3, best_lambda], b[3]),
                      err_f_1=rmse(eval.FEM(FEM(mask_smooth(output_fdaPDE$fit.FEM.mixed$coeff[,best_lambda],1), 
                                                FEMbasis), test.locations), test.func1),
                      err_f_2=rmse(eval.FEM(FEM(mask_smooth(output_fdaPDE$fit.FEM.mixed$coeff[,best_lambda],2), 
                                                FEMbasis), test.locations), test.func2),
                      err_f_3=rmse(eval.FEM(FEM(mask_smooth(output_fdaPDE$fit.FEM.mixed$coeff[,best_lambda],3), 
                                                FEMbasis), test.locations), test.func3)))
    
    estimates$fdaPDE[[as.character(n_obs[j])]][,i] <- c(mask_smooth(output_fdaPDE$fit.FEM.mixed$coeff[,best_lambda],1),
                                                        mask_smooth(output_fdaPDE$fit.FEM.mixed$coeff[,best_lambda],2),
                                                        mask_smooth(output_fdaPDE$fit.FEM.mixed$coeff[,best_lambda],3))
    
    # mgcv (soap film smoother) --------------------------------------------------
    
    mgcv_loc = rbind(locations, locations, locations)
    mgcv_X1 = c(X1[,1], X2[,1], X3[,1])
    mgcv_X2 = c(X1[,2], X2[,2], X3[,2])
    mgcv_obs = observations
    unit1 <- rep(c(1,0,0),c(nlocs,nlocs,nlocs))
    unit2 <- rep(c(0,1,0),c(nlocs,nlocs,nlocs))
    unit3 <- rep(c(0,0,1),c(nlocs,nlocs,nlocs))
    fac <- as.factor(rep(c(1,2,3), c(nlocs, nlocs, nlocs)))
    fac <- factor(rep(c(1,2,3), c(nlocs, nlocs, nlocs)), levels = c(1,2,3))
    
    mgcv_data <- data.frame(x = mgcv_loc[,1],
                            y = mgcv_loc[,2],
                            obs = as.vector(observations),
                            cov1 = mgcv_X1,
                            cov2 = mgcv_X2,
                            unit1 = unit1, unit2 = unit2, unit3 = unit3,
                            fac=fac)
    
    # output_mgcv <- gam(obs ~ cov1+ cov2 + #s(cov1, by=fac) + s(cov2,by=fac)+
    #                      ti(cov1, by=unit1, bs="re") + ti(cov2,by=unit1, bs="re") + # s(cov1,by=fac,bs="sz") + s(cov2,by=fac,bs="sz")
    #                      ti(cov1, by=unit2, bs="re") + ti(cov2,by=unit2, bs="re") +
    #                      ti(cov1, by=unit3, bs="re") + ti(cov2,by=unit3, bs="re") +
    #                      #s(cov1,by=fac,bs="sz") + s(cov2,by=fac,bs="sz")+
    #                      s(x, y, by=unit1, bs = "sf", # soap-film smoother #, id=1,
    #                        xt = list(bnd = bound)) +
    #                      s(x, y, by=unit2, bs = "sf", # soap-film smoother #, id=1,
    #                        xt = list(bnd = bound)) +
    #                      s(x, y, by=unit3, bs = "sf", # soap-film smoother #, id=1,
    #                        xt = list(bnd = bound)) -1, #random=list(cov1=~1, cov2=~1),
    #                    data = mgcv_data)
    
    # NB s(x,y, by=fac) è corretta ma DEVI aggiungere fac come parametro perché le s_i(x,y) hanno media nulla !
    output_mgcv <- gam(obs ~ cov1+ cov2 + fac + cov1:fac + #cov1:unit1 + cov1:unit2 + cov1:unit3 +
                         s(x, y, by=fac, bs = "sf", id=1,# soap-film smoother #, id=1,
                           xt = list(bnd = bound)) -1 ,
                       data = mgcv_data)
    
    # guarda quaderno
    A = matrix(c(1,1,0,0, 
                 0,-1,1,0,
                 0,-1,0,1,
                 0,1,1,1), nrow = 4, ncol=4, byrow = T) 
    
    rhs = c(output_mgcv$coefficients[c(1,6:7)], 0)
    mgcv_coeffs <- solve(A,rhs)
  
    mgcv_loc = rbind(test.locations, test.locations, test.locations)
    mgcv_X1 = c(test.X1[,1], test.X2[,1], test.X3[,1])
    mgcv_X2 = c(test.X1[,2], test.X2[,2], test.X3[,2])
    unit1 <- rep(c(1,0,0),c(nrow(test.locations),nrow(test.locations),nrow(test.locations)))
    unit2 <- rep(c(0,1,0),c(nrow(test.locations),nrow(test.locations),nrow(test.locations)))
    unit3 <- rep(c(0,0,1),c(nrow(test.locations),nrow(test.locations),nrow(test.locations)))
    mgcv_nlocs = nrow(test.locations)*3
    fac <- as.factor(rep(c(1,2,3), c(nrow(test.locations), nrow(test.locations), nrow(test.locations))))
    
    test.observations <- c(test.X1%*% betas + test.func1 + test.X1[,1] * b[1],
                           test.X2%*% betas + test.func2 + test.X2[,1] * b[2],
                           test.X3%*% betas + test.func3 + test.X3[,1] * b[3])
    
    mgcv_test <- data.frame(x = mgcv_loc[,1],
                            y = mgcv_loc[,2],
                            obs = test.observations,
                            cov1 = mgcv_X1,
                            cov2 = mgcv_X2,
                            unit1=as.numeric(unit1), unit2=as.numeric(unit2), unit3=as.numeric(unit3),
                            fac=fac)
    
    tmp1 <- predict.gam(output_mgcv,  newdata = mgcv_test[1:nrow(test.locations),],
                       type="terms", terms="s(x,y):fac1")
    tmp1 <- tmp1 + output_mgcv$coefficients[3]
    
    tmp2 <- predict.gam(output_mgcv, newdata = mgcv_test[((nrow(test.locations)+1) : (2*nrow(test.locations))), ],
                       type="terms", terms="s(x,y):fac2")
    tmp2 <- tmp2 + output_mgcv$coefficients[4]
        tmp3 <- predict.gam(output_mgcv,  newdata = mgcv_test[((2*nrow(test.locations)+1) : (3*nrow(test.locations))), ],
                       type="terms", terms="s(x,y):fac3")
    tmp3 <- tmp3 + output_mgcv$coefficients[5]
  
    results =  rbind(results, 
                     data.frame(method="mgcv", n_obs = n_obs[j], 
                      beta_1=mgcv_coeffs[1], err_beta_1=rmse(mgcv_coeffs[1], betas[1]),
                      beta_2=output_mgcv$coefficients[2], err_beta_2=rmse(output_mgcv$coefficients[2], betas[2]),
                      alpha_1=mgcv_coeffs[2], err_alpha_1=rmse(mgcv_coeffs[2], b[1]),
                      alpha_2=mgcv_coeffs[3], err_alpha_2=rmse(mgcv_coeffs[3], b[2]),
                      alpha_3=mgcv_coeffs[4], err_alpha_3=rmse(mgcv_coeffs[4], b[3]),
                      err_f_1=rmse(tmp1,test.func1),
                      err_f_2=rmse(tmp2,test.func2),
                      err_f_3=rmse(tmp3,test.func3)))
    
    # mgcv estimates 
    mgcv_loc = rbind(mesh$nodes, mesh$nodes, mesh$nodes)
    mgcv_X1 = c(Cov1(mesh$nodes[,1], mesh$nodes[,2]), Cov1(mesh$nodes[,1], mesh$nodes[,2]), Cov1(mesh$nodes[,1], mesh$nodes[,2]))
    mgcv_X2 = rnorm(3*nnodes, mean = 0, sd = 2)
    unit1 <- rep(c(1,0,0),c(nnodes,nnodes,nnodes))
    unit2 <- rep(c(0,1,0),c(nnodes,nnodes,nnodes))
    unit3 <- rep(c(0,0,1),c(nnodes,nnodes,nnodes))
    fac <- factor(rep(c(1,2,3), each=nnodes), levels = c(1,2,3))
    mgcv_test <- data.frame(x = mgcv_loc[,1],
                            y = mgcv_loc[,2], cov1 = mgcv_X1, cov2 = mgcv_X2,
                            unit1=as.numeric(unit1), unit2=as.numeric(unit2), unit3=as.numeric(unit3),
                            fac=fac)
    
    tmp1 <- predict.gam(output_mgcv,  newdata = mask_smooth(mgcv_test, 1),
                        type="terms", terms="s(x,y):fac1")
    tmp1 <- tmp1 + output_mgcv$coefficients[3]
    tmp2 <- predict.gam(output_mgcv, newdata = mask_smooth(mgcv_test, 2),
                        type="terms", terms="s(x,y):fac2")
    tmp2 <- tmp2 + output_mgcv$coefficients[4]
    tmp3 <- predict.gam(output_mgcv,  newdata = mask_smooth(mgcv_test, 3),
                        type="terms", terms="s(x,y):fac3")
    tmp3 <- tmp3 + output_mgcv$coefficients[5]
    estimates$mgcv[[as.character(n_obs[j])]][,i] <- c(tmp1, tmp2, tmp3)
  }
}

results$method <- as.factor(results$method)
results$n_obs <- as.factor(results$n_obs)
save(results, file = paste0(folder.name, "results.RData"))
save(estimates, file = paste0(folder.name, "estimates.RData"))

plot_boxplot(results, n_obs="n_obs", method="method", 
             filename = paste0(folder.name,"boxplots.pdf"))


source("../utils/plot_smooth_2D.R")
# 100 obs
plot_smooth_2D(FEM(mask_smooth(rowMeans(estimates$fdaPDE[[1]]),1), FEMbasis)) # :)
plot_smooth_2D(FEM(mask_smooth(rowMeans(estimates$fdaPDE[[1]]),2), FEMbasis))
plot_smooth_2D(FEM(mask_smooth(rowMeans(estimates$fdaPDE[[1]]),3), FEMbasis))
# post-proc --------------------------------------------------------------------

at_ <- c(1:2, 4:5, 7:8, 10:11)
fill_col <- viridis::viridis(2, begin=0.25, end=0.95)
legend <- c("fdaPDE", "mgcv")
{
  mai_ = par("mai")
  mai_[2] = mai_[2] + 0.075
  pdf(paste0(folder.name, "test_1.pdf"), family = "serif", width = 7, height = 7)
  par(mai=mai_)
  boxplot(results$err_beta_1 ~ results$method + as.numeric(results$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(beta[1]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mai=mai_)
  boxplot(results$err_beta_2 ~ results$method + as.numeric(results$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(beta[2]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mai=mai_)
  boxplot(results$beta_1 ~results$method + as.numeric(results$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(hat(beta)[1]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  abline(h=betas[1], lty=2, lwd=3, col="red")
  par(mai=mai_)
  boxplot(results$beta_2 ~ results$method + as.numeric(results$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(hat(beta)[2]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  abline(h=betas[2], lty=2, lwd=3, col="red")
  par(mai=mai_)
  boxplot(results$err_f_1 ~ results$method + as.numeric(results$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(f[1]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mai=mai_)
  boxplot(results$err_f_2 ~ results$method + as.numeric(results$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(f[2]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  par(mai=mai_)
  boxplot(results$err_f_3 ~ results$method + as.numeric(results$n_obs),
          ylab="RMSE", xlab="observations", at = at_, xaxt="n",
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main =expression(f[3]))
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  dev.off()
}
