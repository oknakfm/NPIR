## \{ y_i | |x_{ik}-x| < \epsilon \}
Empirical_level_set <- function(X = X, Y = Y, y, k = NULL, epsilon){
  if(k>d){
    warning("k exceeds the dimension d so we use k=NULL instead")
    k = NULL
  } 
  if(is.null(k)){
    id = which( pdist(y,Y)@dist < epsilon )
  }else{
    id = which( abs(y-Y[,k]) < epsilon )
  }
  return(list(id=id, X=X[id,], Y=Y[id,]))
}

## omega : B(\infty) \to B(2)
omega <- function(Y){
  omega_entrywise <- function(y){
    if(sum(y^2)==0) c(0,0) else y * max(abs(y)) / sqrt(sum(y^2))
  }
  if(is.null(nrow(Y))){
    omega_entrywise(Y)
  }else{
    t(apply(Y,1,omega_entrywise))
  }
} 

## omega^{-1} : B(2) \to B(\infty)
omega_inverse <- function(Y){
  omega_inv_entrywise <- function(y){
    if(sum(y^2)==0) c(0,0) else y * sqrt(sum(y^2)) / max(abs(y))
  }
  if(is.null(nrow(Y))){
    omega_inv_entrywise(Y)
  }else{
    t(apply(Y,1,omega_inv_entrywise))
  }
}

## angles : R^d to [0,2\pi) from the vector e=(0,1)
angles <- function(Y){
  theta <- function(y) (-Arg(complex(real = y[1], imaginary = y[2])) + pi/2) %% (2*pi)
  if(is.null(nrow(Y))) theta(Y) else apply(Y,1,theta) 
} 

## rotate vectors theta degrees
rotation <- function(Y,theta) Y %*%  matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)

## Polar to Cartesian: vectors(r=||x|| , theta=angles(x)) = y
vectors <- function(r,theta)  c(r * sin(theta), r * cos(theta))

## visualizing the data points
visualizer <- function(type = "square",  ## "square" or "sphere"
                       data, ## data[[1]], data[[2]], ..., data[[n]] are visualized differently for the list "data"
                       colors = NULL, pch = NULL,
                       legend = NULL, ## names of data[[1]], data[[2]], ..., and ignored if NULL 
                       title = NULL ## plot title
){
  L = length(data)
  xl = yl = c(-1.2,1.2)
  plot(0, type="n", xlim=xl, ylim=yl, xlab="x1", ylab="x2", main = title)
  
  if(is.null(colors)) colors = 1:L
  if(is.null(pch)) pch = 1:L
  
  if(type=="square"){
    segments(-1, -1, 1, -1); segments(1, -1, 1, 1); segments(1, 1, -1, 1); segments(-1, 1, -1, -1)
  }else if(type=="sphere"){
    theta = seq(0, 2 * pi, length = 200) # angles for drawing points around the circle
    lines(x = cos(theta), y = sin(theta))
  }else{
    warnings("no matched type")
  }
  
  for(l in 1:L){
    par(new=T)
    plot(data[[l]], xlim=xl, ylim=yl, xlab=" ", ylab=" ", xaxt="n", yaxt="n", col=colors[l], pch=pch[l], type="p")
  }
  
  if(!is.null(legend)) legend("topright", legend=legend, col=colors, pch=pch, bg="white")
}

visualizer_demo <- function(type = "square", X = X, Y = Y, k = 1, epsilon = epsilon, title = "demo"){
  
  sn1  = Empirical_level_set(X = X, Y = Y, y = -1,   k = k, epsilon = epsilon) ## f_1^{-1}(-1)
  sp1  = Empirical_level_set(X = X, Y = Y, y =  1,   k = k, epsilon = epsilon) ## f_1^{-1}(+1)
  
  ## visualization on sphere
  sn05 = Empirical_level_set(X = X, Y = Y, y = -0.5, k = k, epsilon = epsilon) ## f_1^{-1}(-0.5)
  s0   = Empirical_level_set(X = X, Y = Y, y = 0,    k = k, epsilon = epsilon) ## f_1^{-1}(0)
  sp05 = Empirical_level_set(X = X, Y = Y, y = 0.5,  k = k, epsilon = epsilon) ## f_1^{-1}(0.5)
  
  visualizer(type=type, 
             data=list(sn1=sn1$X, sn05=sn05$X, s0=s0$X, sp05=sp05$X, sp1=sp1$X), 
             legend=c("-1","-0.5","0","0.5","1"), colors=1:5, pch=c(1,1,2,3,3),title=title)
  
}

## k-nearest neighbour estimator
kNN <- function(X_obs = X, Y_obs = Y, X_points = c(0,0), k = 10, truncate = FALSE){
  D = as.matrix(pdist(X = X_obs, Y = X_points))
  W = Matrix(D*0,sparse=T)
  ind <- function(i) order(D[,i])[1:k]
  for(i in 1:ncol(D)) W[ind(i),i]=1
  
  Y_pred = 
    if(is.null(nrow(X_points))){
      apply(Y_obs[which(W[,1]==1),],2,mean)
    }else{
      t(sapply(1:nrow(X_points), function(j) apply(Y_obs[which(W[,j]==1),],2,mean)))
    }
  
  if(truncate){
    Y_pred[which(Y_pred < -1, arr.ind = T)] = -1
    Y_pred[which(Y_pred > 1, arr.ind = T)] = 1
  }
  
  return(Y_pred)
}
