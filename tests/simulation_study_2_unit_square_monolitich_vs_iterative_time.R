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

betas <- as.matrix(c(3,0.5))
b <- as.matrix( c(-5., 0, 5.) )

n_sim <- 30

rmse <- function(x,y){
  return(sqrt(mean( (x-y)^2 ) ) ) 
}

n <- c(16, 32, 64, 128)  

# Building folders -------------------------------------------------------------
date_ = unlist(strsplit(as.character(gsub(":","_",gsub(" ","-",Sys.time()))), 
                        split = "[.]"))[1]

#date_ <- "2024-05-14-20_34_15"

if(!dir.exists("data/")) {
  dir.create("data/")
}

if( !dir.exists("data/simulation_2/")){
  dir.create("data/simulation_2/")
}

folder.name = paste("data/simulation_2/",date_,"/",sep="")

if(!dir.exists(folder.name)) {
  dir.create(folder.name)
}

# locations at nodes -----------------------------------------------------------
set.seed(0)
for(j in 1:length(n)){
  datafolder = paste0("data/unit_square_", n[j], "/")
  nodes <- read.csv(paste0(datafolder, "points.csv"))[,2:3]
  elements <- read.csv(paste0(datafolder, "elements.csv"))[,2:4]
  mesh <- create.mesh.2D(nodes=nodes, triangles=elements)
  FEMbasis <- create.FEM.basis(mesh)
  
  locations = mesh$nodes
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  
  func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
  
  results <- list(
    time = matrix(0, nrow=n_sim,ncol=1)
  )
  
  results <- list(monolithic=results, richardson=results, iterative=results)
  
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
    start_ <- Sys.time()
    invisible(capture.output(output_iterative <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, lambda = 1, 
                                                                               #lambda.selection.criterion = "grid", 
                                                                               #lambda.selection.lossfunction = "GCV",
                                                                               #DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE)))
    results$iterative$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_richardson <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, lambda = 1, 
                                                                               #lambda.selection.criterion = "grid", 
                                                                               #lambda.selection.lossfunction = "GCV",
                                                                               #DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE, anderson_memory = 1L)))
    results$richardson$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_monolitich <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                                covariates = X, random_effect = c(1),
                                                                                FEMbasis = FEMbasis, lambda = 1, 
                                                                                #lambda.selection.criterion = "grid", 
                                                                                #lambda.selection.lossfunction = "GCV",
                                                                                #DOF.evaluation = "exact", 
                                                                                FLAG_ITERATIVE = FALSE)))
    results$monolitich$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
  
  }
  save(results,
       file = paste0(folder.name, "n_", n[j],".RData"))
}

# post-proc 
n_obs <- vector(mode="integer", length = 4)
h <- 1 / n

for(j in 1:length(n)){
  datafolder = paste0("data/unit_square_", n[j], "/")
  nodes <- read.csv(paste0(datafolder, "points.csv"))[,2:3]
  n_obs[j] <- nrow(nodes)  
}


estimates <- data.frame(
  time = matrix(nrow = n_sim * length(n_obs), ncol = 1),
  n_obs = as.factor(rep(n_obs, each=n_sim))
)

estimates2 <- estimates
estimates3 <- estimates
for(i in 1:length(n)){
  load(file = paste0(folder.name, "n_", n[i],".RData"))
  
  estimates$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("monolithic", times=n_sim)
  estimates$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$time
  
  estimates2$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("iterative", times=n_sim)
  estimates2$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$time
  
  estimates3$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("richardson", times=n_sim)
  estimates3$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$richardson$time
}

estimates <- rbind(estimates, estimates2, estimates3)
estimates$method <- as.factor(estimates$method)

#at_ <- c(1:2, 4:5, 7:8, 10:11)
at_ <- c(1:3, 5:7, 9:11, 13:15)
fill_col <- viridis::viridis((length(levels(estimates$method))+1), begin=0.25, end=0.95)
fill_col <- fill_col[1:length(levels(estimates$method))]
legend <- levels(estimates$method)

{
  mar_ = par("mar")
  mar_[2] = mar_[2] + 0.25
  pdf(paste0(folder.name, "simulation_2.pdf"), family = "serif", width = 7, height = 7)
  par(mar=mar_)
  boxplot(estimates$time ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          ylim=c(min(estimates$time), max(estimates$time)),
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main ="execution time [s]")
  axis(side = 1, at = c(1.5, 4.5, 7.5, 10.5), labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  
  monolithic <- estimates[estimates$method == "monolithic", ]
  iterative <- estimates[estimates$method == "iterative", ]
  
  plot(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylab="", xlab="log(nodes)", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time [s]")
  points(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
         pch=16, cex=2, col=fill_col[2])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         pch=16, cex=2, col=fill_col[1])
  points(log(n_obs), log(n_obs)-10,
         type="l", lty=3, lwd=3, col="black")
  points(log(n_obs), 2*log(n_obs)-16,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = log(n_obs), labels = round(log(n_obs), digits = 2), cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  plot(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylab="", xlab="nodes", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time [s]")
  points(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
         pch=16, cex=2, col=fill_col[2])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         pch=16, cex=2, col=fill_col[1])
  points(n_obs, n_obs,
         type="l", lty=3, lwd=3, col="black")
  points(n_obs, n_obs,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  
  dev.off()
  
}

# nnodes fixed -----------------------------------------------------------------
n <- 128
n_obs <- c(1000, 2000, 4000, 8000, 16000)

datafolder = paste0("data/unit_square_", n, "/")
nodes <- read.csv(paste0(datafolder, "points.csv"))[,2:3]
elements <- read.csv(paste0(datafolder, "elements.csv"))[,2:4]
mesh <- create.mesh.2D(nodes=nodes, triangles=elements)
FEMbasis <- create.FEM.basis(mesh)
nnodes <- nrow(mesh$nodes)

set.seed(0)
for(j in 1:length(n_obs)){
  nlocs <- n_obs[j]
  locations <- cbind(runif(nlocs), runif(nlocs))
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  
  func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
  
  results <- list(
    time = matrix(0, nrow=n_sim,ncol=1)
  )
  
  results <- list(monolithic=results, iterative=results, richardson=results)
  
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
    start_ <- Sys.time()
    invisible(capture.output(output_iterative <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                               locations = locations,
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, lambda = 1, 
                                                                               #lambda.selection.criterion = "grid", 
                                                                               #lambda.selection.lossfunction = "GCV",
                                                                               #DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE)))
    results$iterative$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_richardson <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                               locations = locations,
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, lambda = 1, 
                                                                               #lambda.selection.criterion = "grid", 
                                                                               #lambda.selection.lossfunction = "GCV",
                                                                               #DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE, anderson_memory = 1L)))
    results$richardson$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    
    start_ <- Sys.time()
    invisible(capture.output(output_monolitich <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                                locations = locations,
                                                                                covariates = X, random_effect = c(1),
                                                                                FEMbasis = FEMbasis, lambda = 1, 
                                                                                #lambda.selection.criterion = "grid", 
                                                                                #lambda.selection.lossfunction = "GCV",
                                                                                #DOF.evaluation = "exact", 
                                                                                FLAG_ITERATIVE = FALSE)))
    results$monolitich$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
  }
  save(results,
       file = paste0(folder.name, "n_obs_", n_obs[j],".RData"))
}


estimates <- data.frame(
  time = matrix(nrow = n_sim * length(n_obs), ncol = 1),
  n_obs = as.factor(rep(n_obs, each=n_sim))
)

estimates2 <- estimates

for(i in 1:length(n_obs)){
  load(file = paste0(folder.name, "n_obs_", n_obs[i],".RData"))
  
  estimates$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("monolithic", times=n_sim)
  estimates$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$time
  
  estimates2$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("iterative", times=n_sim)
  estimates2$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$time
}

estimates <- rbind(estimates, estimates2)
estimates$method <- as.factor(estimates$method)
at_ <- c(1:2, 4:5, 7:8, 10:11, 12:13)
fill_col <- viridis::viridis(3, begin=0.25, end=0.95)
fill_col <- fill_col[1:2]
legend <- levels(estimates$method)

{
  mar_ = par("mar")
  mar_[2] = mar_[2] + 0.25
  pdf(paste0(folder.name, "simulation_2_pt2.pdf"), family = "serif", width = 7, height = 7)
  par(mar=mar_)
  boxplot(estimates$time ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          ylim=c(min(estimates$time), max(estimates$time)),
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main ="execution time [s]")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  
  monolithic <- estimates[estimates$method == "monolithic", ]
  iterative <- estimates[estimates$method == "iterative", ]
  
  plot(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylim= c(min(log(estimates$time)), max(log(estimates$time))),
       ylab="", xlab="log(nodes)", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time [s]")
  points(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
         pch=16, cex=2, col=fill_col[2])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         pch=16, cex=2, col=fill_col[1])
  points(log(n_obs), log(n_obs)-10,
         type="l", lty=3, lwd=3, col="black")
  points(log(n_obs), 2*log(n_obs)-16,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = log(n_obs), labels = round(log(n_obs), digits = 2), cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  plot(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
       ylim= c(min(estimates$time), max(estimates$time)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylab="", xlab="nodes", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time [s]")
  points(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
         pch=16, cex=2, col=fill_col[2])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         pch=16, cex=2, col=fill_col[1])
  points(n_obs, log(n_obs),
         type="l", lty=3, lwd=3, col="black")
  points(n_obs, n_obs,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  
  dev.off()
  
}

# mesh nodes fixed - grid GCV --------------------------------------------------
n <- 32
n_obs <- c(250, 500, 1000, 2000, 4000)

datafolder = paste0("data/unit_square_", n, "/")
nodes <- read.csv(paste0(datafolder, "points.csv"))[,2:3]
elements <- read.csv(paste0(datafolder, "elements.csv"))[,2:4]
mesh <- create.mesh.2D(nodes=nodes, triangles=elements)
FEMbasis <- create.FEM.basis(mesh)
nnodes <- nrow(mesh$nodes)

lambda <- 10^seq(-3,0,length=10)

set.seed(0)
for(j in 1:length(n_obs)){
  nlocs <- n_obs[j]
  locations <- cbind(runif(nlocs), runif(nlocs))
  nlocs <- nrow(locations)
  X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  
  func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
  
  results <- list(
    time = matrix(0, nrow=n_sim,ncol=1)
  )
  
  results <- list(monolithic=results, iterative=results, richardson=results)
  
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
    start_ <- Sys.time()
    invisible(capture.output(output_iterative <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                               locations = locations,
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, lambda = lambda, 
                                                                               lambda.selection.criterion = "grid", 
                                                                               lambda.selection.lossfunction = "GCV",
                                                                               DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE)))
    results$iterative$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_richardson <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                               locations = locations,
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, lambda = lambda, 
                                                                               lambda.selection.criterion = "grid", 
                                                                               lambda.selection.lossfunction = "GCV",
                                                                               DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE, anderson_memory = 1L)))
    results$richardson$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_monolitich <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, 
                                                                                locations = locations,
                                                                                covariates = X, random_effect = c(1),
                                                                                FEMbasis = FEMbasis, lambda = lambda, 
                                                                                lambda.selection.criterion = "grid", 
                                                                                lambda.selection.lossfunction = "GCV",
                                                                                DOF.evaluation = "exact", 
                                                                                FLAG_ITERATIVE = FALSE)))
    results$monolitich$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
  }
  save(results,
       file = paste0(folder.name, "GCV_n_obs_", n_obs[j],".RData"))
}

estimates <- data.frame(
  time = matrix(nrow = n_sim * length(n_obs), ncol = 1),
  n_obs = as.factor(rep(n_obs, each=n_sim))
)


for(i in 1:length(n_obs)){
  load(file = paste0(folder.name, "GCV_n_obs_", n_obs[i],".RData"))
  
  estimates$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("monolithic", times=n_sim)
  estimates$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$time
  
  estimates2$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("iterative", times=n_sim)
  estimates2$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$time
}

estimates <- rbind(estimates, estimates2)

at_ <- c(1:2, 4:5, 7:8, 10:11, 13:14)
fill_col <- viridis::viridis(3, begin=0.25, end=0.95)
fill_col <- fill_col[1:2]
legend <- levels(estimates$method)

{
  mar_ = par("mar")
  mar_[2] = mar_[2] + 0.25
  pdf(paste0(folder.name, "simulation_2_GCV.pdf"), family = "serif", width = 7, height = 7)
  par(mar=mar_)
  boxplot(estimates$time ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          ylim=c(min(estimates$time), max(estimates$time)),
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main ="execution time [s]")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  
  monolithic <- estimates[estimates$method == "monolithic", ]
  iterative <- estimates[estimates$method == "iterative", ]
  
  plot(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylim= c(min(log(estimates$time)), max(log(estimates$time))),
       ylab="", xlab="log(nodes)", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time [s]")
  points(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
         pch=16, cex=2, col=fill_col[2])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         pch=16, cex=2, col=fill_col[1])
  points(log(n_obs), log(n_obs)-10,
         type="l", lty=3, lwd=3, col="black")
  points(log(n_obs), 2*log(n_obs)-16,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = log(n_obs), labels = round(log(n_obs), digits = 2), cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  plot(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
       ylim= c(min(estimates$time), max(estimates$time)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylab="", xlab="nodes", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time [s]")
  points(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
         pch=16, cex=2, col=fill_col[2])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         pch=16, cex=2, col=fill_col[1])
  points(n_obs, log(n_obs),
         type="l", lty=3, lwd=3, col="black")
  points(n_obs, n_obs,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  
  dev.off()
  
}

# n_obs fixed ------------------------------------------------------------------
# GCV newton exact
n <- c(20,30,40,50)#,60)
nlocs <- 1000
locations <- cbind(runif(nlocs), runif(nlocs))
nlocs <- nrow(locations)

scheme = c("iterative", "monolithic", "richardson")
set.seed(0)
for(j in 1:length(n)){
  x <- seq(0,1,length.out=n[j])
  y <- x
  nodes <- expand.grid(x,y)
  mesh <- create.mesh.2D(nodes=nodes)
  nnodes = nrow(mesh$nodes)
  FEMbasis <- create.FEM.basis(mesh)
  X1 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X2 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  X3 = cbind( Cov1(locations[,1], locations[,2]), rnorm(nlocs, mean = 0, sd = 2))
  
  func1 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(0.5, nlocs))
  func2 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1, nlocs))
  func3 = fs.test.time(x=locations[,1], y=locations[,2],  t=rep(1.5, nlocs))
  
  test.locations = refine.by.splitting.mesh.2D(mesh)$nodes
  test.X1 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
  test.X2 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
  test.X3 <- cbind( Cov1(test.locations[,1], test.locations[,2]), rnorm(nrow(test.locations), mean = 0, sd = 2))
  test.func1 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(0.5, nrow(test.locations)))
  test.func2 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1, nrow(test.locations)))
  test.func3 = fs.test.time(x=test.locations[,1], y=test.locations[,2],  t=rep(1.5, nrow(test.locations)))
  
  results <- list(
    time = matrix(0, nrow=n_sim,ncol=1),
    beta_1 = matrix(0, nrow=n_sim,ncol=1),
    beta_2 = matrix(0, nrow=n_sim,ncol=1),
    b_1 = matrix(0, nrow=n_sim,ncol=1), 
    b_2 = matrix(0, nrow=n_sim,ncol=1),
    b_3 = matrix(0, nrow=n_sim,ncol=1),
    f_1 = matrix(0, nrow=n_sim,ncol=1),
    f_2 = matrix(0, nrow=n_sim,ncol=1),
    f_3 = matrix(0, nrow=n_sim,ncol=1)
  )
  
  results <- list(monolithic=results, richardson=results, iterative=results)
  
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
    start_ <- Sys.time()
    invisible(capture.output(output_iterative <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                               covariates = X, random_effect = c(1),
                                                                               FEMbasis = FEMbasis, #lambda = 1, 
                                                                               lambda.selection.criterion = "newton", 
                                                                               lambda.selection.lossfunction = "GCV",
                                                                               DOF.evaluation = "exact",
                                                                               FLAG_ITERATIVE = TRUE)))
    results$iterative$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_richardson <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                                covariates = X, random_effect = c(1),
                                                                                FEMbasis = FEMbasis,# lambda = 1, 
                                                                                lambda.selection.criterion = "newton", 
                                                                                lambda.selection.lossfunction = "GCV",
                                                                                DOF.evaluation = "exact",
                                                                                FLAG_ITERATIVE = TRUE, anderson_memory = 1L)))
    results$richardson$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    start_ <- Sys.time()
    invisible(capture.output(output_monolithic <- fdaPDEmixed::smooth.FEM.mixed(observations = observations, locations = locations,
                                                                                covariates = X, random_effect = c(1),
                                                                                FEMbasis = FEMbasis, #lambda = 1, 
                                                                                lambda.selection.criterion = "newton", 
                                                                                lambda.selection.lossfunction = "GCV",
                                                                                DOF.evaluation = "exact", 
                                                                                FLAG_ITERATIVE = FALSE)))
    results$monolitich$time[i] <- as.numeric(difftime(Sys.time(), start_, units = "secs"))
    
    output <- list(iterative = output_iterative, richardson = output_richardson, monolithic = output_monolithic)
    for(k in scheme){
      best_lambda = output[[k]]$bestlambda
      results[[k]]$beta_1[i,] <- rmse(output[[k]]$beta[1], betas[1])
      results[[k]]$beta_2[i,] <- rmse(output[[k]]$beta[2], betas[2])
      results[[k]]$b_1[i,] <- rmse(output[[k]]$b_i[1], b[1])
      results[[k]]$b_2[i,] <- rmse(output[[k]]$b_i[2], b[2])
      results[[k]]$b_3[i,] <- rmse(output[[k]]$b_i[3], b[3])
      results[[k]]$f_1[i] <- rmse(eval.FEM(FEM(as.matrix(output[[k]]$fit.FEM.mixed$coeff[1:nnodes,best_lambda]), 
                                              FEMbasis), test.locations),
                                       test.func1)
      results[[k]]$f_2[i] <- rmse(eval.FEM(FEM(as.matrix(output[[k]]$fit.FEM.mixed$coeff[(nnodes+1):(2*nnodes),best_lambda]),
                                              FEMbasis), test.locations),
                                       test.func2)
      results[[k]]$f_3[i] <- rmse(eval.FEM(FEM(as.matrix(output[[k]]$fit.FEM.mixed$coeff[(2*nnodes+1):(3*nnodes),best_lambda]), 
                                              FEMbasis), test.locations),
                                       test.func2)
      
    }
    
  }
  save(results,
       file = paste0(folder.name, "GCV_newton_n_", n[j],".RData"))
}

n <- c(20,30,40,50) #,60)
nlocs <- 1000
n_sim= 30
#folder.name <- "data/simulation_2/2024-05-20-09_39_03/"

estimates <- data.frame(
  time = matrix(nrow = n_sim * length(n), ncol = 1),
  n_obs = as.factor(rep(n, each=n_sim))
)

estimates2 <- estimates
estimates3 <- estimates
for(i in 1:length(n)){
  load(file = paste0(folder.name, "GCV_newton_n_", n[i],".RData"))
  
  estimates$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("monolithic", times=n_sim)
  estimates$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$monolitich$time
  
  estimates2$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("iterative", times=n_sim)
  estimates2$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$iterative$time
  
  estimates3$method[(1+(i-1)*n_sim):(i*n_sim)] <- rep("richardson", times=n_sim)
  estimates3$time[(1+(i-1)*n_sim):(i*n_sim)] <- results$richardson$time
}

estimates <- rbind(estimates, estimates2, estimates3)
estimates$method <- as.factor(estimates$method)

#at_ <- c(1:2, 4:5, 7:8, 10:11)
at_ <- c(1:3, 5:7, 9:11, 13:15)
fill_col <- viridis::viridis((length(levels(estimates$method))+1), begin=0.25, end=0.95)
fill_col <- fill_col[1:length(levels(estimates$method))]
legend <- levels(estimates$method)

estimates <- rbind(estimates, estimates2)
legend <- levels(estimates$method)
n_obs <- n
{
  mar_ = par("mar")
  mar_[2] = mar_[2] + 0.25
  pdf(paste0(folder.name, "simulation_2_GCV_newton.pdf"), family = "serif", width = 7, height = 7)
  par(mar=mar_)
  boxplot(estimates$time ~ estimates$method + as.numeric(estimates$n_obs),
          ylab="", xlab="observations", at = at_, xaxt="n",
          ylim=c(min(estimates$time), max(estimates$time)),
          col=fill_col,cex.lab = 2, cex.axis = 2, cex.main = 2,
          main ="execution time [s]")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topright",legend=legend, fill=fill_col, horiz=T, cex=1.5, inset=0.0125, 
         bty="n")
  
  monolithic <- estimates[estimates$method == "monolithic", ]
  iterative <- estimates[estimates$method == "iterative", ]
  richardson <- estimates[estimates$method == "richardson", ]
  
  plot(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylim= c(min(log(estimates$time)), max(log(estimates$time))),
       ylab="log(times)", xlab="log(nodes)", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time")
  points(log(n_obs), log(tapply(monolithic$time, monolithic$n_obs, mean)),
         pch=16, cex=2, col=fill_col[2])
  
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(log(n_obs), log(tapply(iterative$time, iterative$n_obs, mean)),
         pch=16, cex=2, col=fill_col[1])
  
  points(log(n_obs), log(tapply(richardson$time, richardson$n_obs, mean)),
         type="l", lty=4, lwd=4, col=fill_col[3])
  points(log(n_obs), log(tapply(richardson$time, richardson$n_obs, mean)),
         pch="*", cex=2, col=fill_col[3])
  
  points(log(n_obs), log(n_obs)-10,
         type="l", lty=3, lwd=3, col="black")
  points(log(n_obs), 2*log(n_obs)-16,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = log(n_obs), labels = round(log(n_obs), digits = 2), cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  plot(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
       ylim= c(min(estimates$time), max(estimates$time)),
       type="l", lty=2, lwd=4, col=fill_col[2], 
       ylab="time [s]", xlab="nodes", xaxt="n",
       cex.lab = 2, cex.axis = 2, cex.main = 2,
       main ="execution time")
  points(n_obs, tapply(monolithic$time, monolithic$n_obs, mean),
         pch=16, cex=2, col=fill_col[2])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         type="l", lty=2, lwd=4, col=fill_col[1])
  points(n_obs, tapply(iterative$time, iterative$n_obs, mean),
         pch=16, cex=2, col=fill_col[1])
  points(n_obs, log(tapply(richardson$time, richardson$n_obs, mean)),
         type="l", lty=4, lwd=4, col=fill_col[3])
  points(n_obs, log(tapply(richardson$time, richardson$n_obs, mean)),
         pch="*", cex=2, col=fill_col[3])
  points(n_obs, log(n_obs),
         type="l", lty=3, lwd=3, col="black")
  points(n_obs, n_obs,
         type="l", lty=3, lwd=3, col="red")
  axis(side = 1, at = n_obs, labels = n_obs, cex.lab = 2, cex.axis = 2)
  legend("topleft", legend=legend, fill=fill_col, horiz=F, cex=1.5, inset=0.0125, 
         bty="n")
  
  
  dev.off()
  
}

