library("plotrix"); library("latex2exp"); library("RColorBrewer");
library("viridis")

contour.FEM <- function(x, limits=NULL, ...){
  mesh <- x$FEMbasis$mesh
  plot_data <- data.frame(  X=mesh$nodes[,1], 
                            Y=mesh$nodes[,2],
                            Z=rep(0, times=nrow(mesh$nodes)),
                            coeff=x$coeff[1:nrow(mesh$nodes)])
  
  coeff <- apply(mesh$triangles, MARGIN=1, FUN = function(edge){
    mean(x$coeff[edge,])
  })
  
  if(is.null(limits)) limits = c(min(coeff), max(coeff))
  cmin = limits[1]; cmax=limits[2]
  
  I=mesh$triangles[,1]-1
  J=mesh$triangles[,2]-1
  K=mesh$triangles[,3]-1
  fig<- plot_ly(plot_data, x=~X, y=~Y, z=~Z,
                i = I, j = J, k = K, cmin = limits[1], cmax=limits[2],
                intensity=~coeff, color=~coeff, type="mesh3d", 
                colorbar=list(title=""), ...) %>%
    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      yaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      zaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      camera = list(
        eye = list(x = 0, y = -0.01,  z = 1.25))),  dragmode="zoom") %>%
    colorbar(len = 1, title="")
}

# sample_locations <- function(mesh, N){
#   
#   gridnum = 2 * floor( sqrt(N) )
#   minX = min(mesh$nodes[,1])
#   maxX = max(mesh$nodes[,1])
#   minY = min(mesh$nodes[,2])
#   maxY = max(mesh$nodes[,2])
#   
#   x = seq(minX+0.125, maxX-0.125, length.out = gridnum)
#   y = seq(minY+0.125, maxY-0.125, length.out = gridnum)
#   unif_grid <- expand.grid(x = x, y = y)
#   idx <- fdaPDE:::CPP_search_points(mesh, as.matrix(unif_grid))
#   locations <- as.matrix( unif_grid[which(idx != 0), ] )
#   
#   pts_list <- lapply(as.list(as.data.frame(t(locations))), sf::st_point)
#   sf_pts <- sf::st_sfc(pts_list)
#   
#   pts_list <- rbind(mesh$nodes[as.logical(mesh$nodesmarkers),],
#                     mesh$nodes[which(as.logical(mesh$nodesmarkers)==TRUE)[1],])
#   sf_polygon <- sf::st_polygon(list(pts_list))
#   
#   idx <- sf::st_intersects(sf_pts, sf::st_boundary(sf_polygon), sparse = FALSE)
#   if( length( which( idx == TRUE ) ) != 0 ){
#     locations = locations[-which( idx == TRUE ), ]
#   }
#   return(locations[sample(1:nrow(locations), size=N),])
# }

sample_locations.mgcv <- function(mesh, N){
  
  gridnum = 2 * floor( sqrt(N) )
  minX = min(mesh$nodes[,1])
  maxX = max(mesh$nodes[,1])
  minY = min(mesh$nodes[,2])
  maxY = max(mesh$nodes[,2])
  
  x = seq(minX+0.125, maxX-0.125, length.out = gridnum)
  y = seq(minY+0.125, maxY-0.125, length.out = gridnum)
  unif_grid <- expand.grid(x = x, y = y)
  
  idx <- mgcv::fs.test(unif_grid[,1],unif_grid[,2])
  idx <- which(!is.na(idx))
  
  locations <- as.matrix( unif_grid[idx, ] )
  locations <- locations[sample(1:nrow(locations), size=N),]
  rownames(locations) <- rep("", N)
  colnames(locations) <- c("","")
  return(locations)
}

compute_limits <- function(x){
  mesh <- x$FEMbasis$mesh
  coeff <- apply(mesh$triangles, MARGIN=1, FUN = function(edge){
    mean(x$coeff[edge,])
  })
  limits = c(1e10, -1e10)
  limits[1] = min(coeff, limits[1])
  limits[2] = max(coeff, limits[2])
  return(limits)
}

plot.mesh.2D <- function(x, ...){
  plot_data <- data.frame(  X=x$nodes[,1], 
                            Y=x$nodes[,2],
                            Z=rep(0,times=nrow(x$nodes)))
  I=x$triangles[,1]-1
  J=x$triangles[,2]-1
  K=x$triangles[,3]-1
  fig <- plot_ly(...) %>% 
    add_markers(x = x$nodes[,1],
                y = x$nodes[,2],
                color = I('black'), size = I(1)) %>%
    add_segments(x = x$nodes[x$edges[,1],1],
                 y = x$nodes[x$edges[,1],2],
                 xend = x$nodes[x$edges[,2],1],
                 yend = x$nodes[x$edges[,2],2], 
                 color = I('black'), size = I(1),
                 showlegend = F) %>%
    layout( 
      xaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      yaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F)
    )
  fig
}

plot.FEM.2.5D <- function(x, limits=NULL, ...){
  mesh <- x$FEMbasis$mesh
  plot_data <- data.frame(  X=mesh$nodes[,1], 
                            Y=mesh$nodes[,2],
                            Z=mesh$nodes[,3],
                            coeff =x$coeff[1:nrow(mesh$nodes)])
  
  coeff <- apply(mesh$triangles, MARGIN=1, FUN = function(edge){
    mean(x$coeff[edge,])
  })
  I=x$triangles[,1]-1
  J=x$triangles[,2]-1
  K=x$triangles[,3]-1
  if(is.null(limits)) limits = c(min(coeff), max(coeff))
  cmin = limits[1]; cmax=limits[2]
  
  fig<- plot_ly(plot_data, x=~X, y=~Y, z=~Z,
                i = I, j = J, k = K, cmin = limits[1], cmax=limits[2],
                intensity=~coeff, color=~coeff, type="mesh3d", 
                colorbar=list(title=""), ...) %>%
    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      yaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      zaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      camera = list(
        eye = list(x = 0, y = -0.01,  z = 1.25)))) %>%
    colorbar(len = 1, title="")
}

plot.mesh.2.5D <- function(x, ...){
  plot_data <- data.frame(  X=x$nodes[,1], 
                            Y=x$nodes[,2],
                            Z=x$nodes[,3])
  I=x$triangles[,1]-1
  J=x$triangles[,2]-1
  K=x$triangles[,3]-1
  
  plot_ly(plot_data, x = ~X, y = ~Y, z = ~Z, 
          i = I, j = J, k = K, type = 'mesh3d', facecolor = "lightgray", vertexcolor="black") %>%
    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      yaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      zaxis = list(
        title = '', showgrid = F, zeroline = F, showticklabels = F),
      camera = list(
        eye = list(x = 0, y = -0.01,  z = 1.25)))) %>%
    colorbar(len = 1, title="")
}

plot_mesh.2.5D <- function(mesh, ...){
  data_plot <- data.frame(X=mesh$nodes[,1], 
                          Y=mesh$nodes[,2],
                          Z=mesh$nodes[,3])
  
  data_edge <- data.frame(X= mesh$nodes[mesh$edges[,1]][,1],
                          Y= mesh$nodes[mesh$edges[,1]][,2],
                          Z= mesh$nodes[mesh$edges[,1]][,3])
  
  I=(mesh$triangles[,1]-1); J=(mesh$triangles[,2]-1); K=(mesh$triangles[,3]-1)
  fig <- plot_ly(data_plot,
                 type = 'mesh3d', x = ~X, y = ~Y, z = ~Z,
                 i = I, j = J, k = K,
                 hoverinfo = 'none', facecolor="lightgray")%>%
    
    layout(scene = list(
      aspectmode = "data", 
      xaxis = list(title = '',showgrid = F,zeroline = F,showticklabels = F),
      yaxis = list(title = '',showgrid = F,zeroline = F,showticklabels = F),
      zaxis = list(title = '',showgrid = F,zeroline = F,showticklabels = F)))
}

# ------------------------------------------------------------------------------

## GRAPHICAL SETTINGS ----------------------------------------------------------
zoom = 0.7657689 # 0.6
userMatrix = rbind(c(  0.02091786,  0.99873853, -0.04564825,    0),
                   c( -0.13139695,  0.04800860,  0.99016660,    0),
                   c(  0.99110913, -0.01471432,  0.13223548,    0),
                   c(  0.00000000,  0.0000000,   0.0000000,     1))
windowRect = c(70,  106, 1920, 1117)

#pp <- par3d(no.readonly = TRUE)

plot.mesh.2.5D <- function(mesh, M = NULL, m = NULL, ROI=NULL, NA_ = NULL,...){
  
  FEM = FEM(rep(0, nrow(mesh$nodes)), create.FEM.basis(mesh))
  FEM$coeff[ROI,] <- 1 
  FEM$coeff[NA_,] <- 2
  
  
  if (is.null(m)) { m = min(FEM$coeff)}
  if (is.null(M)) { M = max(FEM$coeff)}
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  #p = jet.col(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  p = c("lightgray", "red3", "blue3")
  palette(p)
  
  ncolor = length(p)
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
    diffrange = M - m
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf]))/diffrange*(ncolor-1)+1
    
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col,...)
    rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
              z = nodes[edges,3],
              color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

plot.FEM.2.5D <- function(FEM, limits=NULL, colorscale = jet.col,...){
  
  if (is.null(limits)){
    limits <- c(0,0)
    limits[1] = min(FEM$coeff, na.rm = T)
    limits[2] = max(FEM$coeff, na.rm = T)
  }
  
  triangles = c(t(FEM$FEMbasis$mesh$triangles))
  ntriangles = nrow(FEM$FEMbasis$mesh$triangles)
  order = FEM$FEMbasis$mesh$order
  nodes = FEM$FEMbasis$mesh$nodes
  edges = matrix(rep(0, 6*ntriangles), ncol = 2)
  for(i in 0:(ntriangles-1)){
    edges[3*i+1,] = c(triangles[3*order*i+1], triangles[3*order*i+2])
    edges[3*i+2,] = c(triangles[3*order*i+1], triangles[3*order*i+3])
    edges[3*i+3,] = c(triangles[3*order*i+2], triangles[3*order*i+3])
  }
  edges = edges[!duplicated(edges),]
  edges <- as.vector(t(edges))
  
  coeff = FEM$coeff
  
  FEMbasis = FEM$FEMbasis
  
  mesh = FEMbasis$mesh
  
  p = colorscale(n = 1000, alpha = 0.8)
  # alternative color palette: p <- colorRampPalette(c("#0E1E44", "#3E6DD8", "#68D061", "#ECAF53", "#EB5F5F", "#E11F1C"))(1000)
  
  grDevices::palette(p)
  
  ncolor = length(p)
  
  nsurf = dim(coeff)[[2]]
  for (isurf in 1:nsurf)
  {
    open3d(zoom = zoom, userMatrix = userMatrix, windowRect = windowRect)
    rgl.pop("lights") 
    light3d(specular = "black") 
    
    diffrange = diff(range(limits)) 
    
    col = coeff[triangles,isurf]
    col = (col - min(coeff[,isurf], na.rm =T))/diffrange*(ncolor-1)+1
    if(abs(diffrange) < 1e-10) col = rep(M, times=length(col)) # costanti
    rgl.triangles(x = nodes[triangles ,1], y = nodes[triangles ,2],
                  z = nodes[triangles,3],
                  color = col) #,...)
    # rgl.lines(x = nodes[edges ,1], y = nodes[edges ,2],
    #           z = nodes[edges,3],
    #           color = "black",...)
    aspect3d("iso")
    
    if (nsurf > 1 && isurf<nsurf)
    {readline("Press a button for the next plot...")}
  }
}

plot.colorbar <- function(x, limits= NULL, colorscale = jet.col, 
                          horizontal = FALSE, cex.axis = 2,
                          file = "plot.pdf"){
  
  mesh <- x$FEMbasis$mesh
  
  coeff <- apply(mesh$triangles, MARGIN=1, FUN = function(edge){
    mean(x$coeff[edge,])
  })
  
  if(is.null(limits)) limits = c(min(coeff, na.rm = T), max(coeff, na.rm = T))
  cmin = limits[1]; cmax=limits[2]
  
  exps <- -15:15
  range_ <- round( diff(range(limits)), digits = 2)
  cmin_exp <- which(floor(log10(abs(signif(signif(cmin, digits = 2) / 10^exps, digits = 0))))==0)
  cmax_exp <- which(floor(log10(abs(signif(signif(cmax, digits = 2) / 10^exps, digits = 0))))==0)
  k <- exps[max(cmin_exp, cmax_exp)]
  
  at = seq(0, 100, length.out=5)
  labels = as.character(round(seq(cmin*10^(-k), cmax*10^(-k), length.out=5), 2))
  text_ <- ifelse(k != 0, paste0("$\\times 10^{", k,"}$"), "")
  
  #x11(width = 11, height = 3)
  ncolor = 1000
  diffrange = cmax - cmin 
  if(abs(diffrange) < 1e-10) ncolor = 1 # costanti
  
  labels_grad = ifelse(text_ == "", "", TeX(text_))
  
  pdf(paste0(file, "_horiziontal.pdf"), family = "serif", width = 11, height = 3)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 112.5), c(0, 15), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 100, 5, col = colorscale(ncolor), border = "black")
  axis(1, at = at, labels = labels, # at = c(0,33.33,66.66,100)
       cex.axis = cex.axis, lwd.ticks = 0, lwd = 0) # lwd.ticks = 2, 
  text(107,2, labels_grad, cex = 2)
  dev.off()
  
  pdf(paste0(file, "_vertical.pdf"), family = "serif", width = 3, height = 11)
  #x11(width = 3, height = 11)
  par(mai = c(1,0.75,0,0))
  plot(c(0, 15), c(0, 112.5), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", frame.plot = F)
  gradient.rect(0, 0, 5, 100, col = colorscale(ncolor), border = "black", gradient = "y")
  axis(4, at = at, labels = labels, # at = c(0,33.33,66.66,100)
       cex.axis = cex.axis, lwd.ticks = 0, lwd = 0,  line=-7.5) # lwd.ticks = 2, 
  text(2.5, 107, labels_grad, cex = 2)
  dev.off()
  
}
