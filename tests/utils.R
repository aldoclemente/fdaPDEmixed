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

sample_locations <- function(mesh, N){
  
  gridnum = 2 * floor( sqrt(N) )
  minX = min(mesh$nodes[,1])
  maxX = max(mesh$nodes[,1])
  minY = min(mesh$nodes[,2])
  maxY = max(mesh$nodes[,2])
  
  x = seq(minX+0.125, maxX-0.125, length.out = gridnum)
  y = seq(minY+0.125, maxY-0.125, length.out = gridnum)
  unif_grid <- expand.grid(x = x, y = y)
  idx <- fdaPDE:::CPP_search_points(mesh, as.matrix(unif_grid))
  locations <- as.matrix( unif_grid[which(idx != 0), ] )
  
  pts_list <- lapply(as.list(as.data.frame(t(locations))), sf::st_point)
  sf_pts <- sf::st_sfc(pts_list)
  
  pts_list <- rbind(mesh$nodes[as.logical(mesh$nodesmarkers),],
                    mesh$nodes[which(as.logical(mesh$nodesmarkers)==TRUE)[1],])
  sf_polygon <- sf::st_polygon(list(pts_list))
  # plot(sf_polygon, col="red")  
  # plot(sf_pts, col="black", pch=16, add=T)
  # plot(sf::st_boundary(sf_polygon), col="blue", lwd=3, add = T)
  
  idx <- sf::st_intersects(sf_pts, sf::st_boundary(sf_polygon), sparse = FALSE)
  if( length( which( idx == TRUE ) ) != 0 ){
    locations = locations[-which( idx == TRUE ), ]
  }
  return(locations[sample(1:nrow(locations), size=N),])
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
    layout(scene = list(
      aspectratio=list(x=1,y=1)), 
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
