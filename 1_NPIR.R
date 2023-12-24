DIR_NAME = getwd()
OUTPUT_DIR_NAME = paste0(DIR_NAME,"/Output")
dir.create(OUTPUT_DIR_NAME, showWarnings=FALSE)
setwd(DIR_NAME)

## ===============
## Dependencies
## ===============
require(pdist)
require(Matrix)

# install functions
setwd(DIR_NAME)
source("A_functions.R")

## ================
## Synthetic data
## ================

n = 100^2 # num of samples
d = 2     # covariate dimension
v = 0.001  # V[Y|X] = v * I
# v = 0.1

TWIST = TRUE

for(TWIST in c(TRUE,FALSE)){
  # herein, stay in the directory for output
  if(TWIST){
    tmp_DIR_NAME = paste0(OUTPUT_DIR_NAME,"/TWIST")
  } else{
    tmp_DIR_NAME = paste0(OUTPUT_DIR_NAME,"/NOT_TWIST")
  }
  dir.create(tmp_DIR_NAME, showWarnings=FALSE)
  setwd(tmp_DIR_NAME)

  ## f(x) = E[y|x] : I^d \to I^d
  if(!TWIST){
    deg = pi * 0.3
    f <- function(x){
      z = omega(x)
      z = vectors( r = sqrt(sum(z^2))^(abs(sin(angles(z)))), theta = angles(z) + deg )
      return(omega_inverse(z))
    }
  }else{
    deg = pi * 0.4
    f <- function(x){
      xi <- function(x) if(x <= - 3/4) 4*x + 3 else (4/7)*x + (3/7)
      z = omega(x); z = rotation(z, deg); z = omega_inverse(z)
      z[,1] = sapply(z[,1],xi)
      z = omega(z); z = rotation(z,deg); z = omega_inverse(z)
      z[,1] = sapply(z[,1],xi)
      z = omega(z); z = rotation(z,deg); z = omega_inverse(z)
      return(z)
    }
  }

  X = matrix(runif(n*d, min=-1, max=1),n,d) # covariates
  E = matrix(rnorm(n*d, mean=0, sd=sqrt(v)),n,d)  # error
  Y = t(apply(X,1,f)) + E # outcomes
  
  save(file = "D1_Synthetic_Data.RData",X,E,Y,deg,f)
  
  ### =============================================
  ### Estimating a homeomorphism \rho: I^2 \to I^2
  ### =============================================
  
  ## [== Estimation ==]
  
  # vertices of the square I^2
  X_vertices = rbind(c(1,1),c(1,-1),c(-1,-1),c(-1,1))
  
  # first-step estimator f^1(x) for the vertices x
  f1_vertices = kNN(X_obs = X, Y_obs = Y, X_points = X_vertices, k = 10, truncate = FALSE)
  
  # f^1(x) \in I^2 mapped to D^2
  zeta = omega(f1_vertices)
  
  # [[z]] = x \in [0,2\pi] satisfying z-x = 2*pi*m
  md <- function(z) z %% (2*pi)
  
  # theta is selected so that 0 < [[angle(f(x1)) + theta]] < [[angle(f(x4)) + theta]] < 2*pi
  theta = md( md(2*pi - angles(zeta[1,])) + md(angles(zeta[1,]) - angles(zeta[4,])) / 2 )
  
  # eta = v(1, [[angle(zeta)+theta]])
  eta_angle = (angles(zeta) + theta) %% (2*pi)
  eta = NULL; for(j in 1:4) eta[[j]] = vectors(1,eta_angle[j])
  
  # tau:[0,2\pi] to [0,2pi]
  tau <- function(x){
    t1 = eta_angle[1]; t2 = eta_angle[2]; t3 = eta_angle[3]; t4 = eta_angle[4]
    tau_entrywise <- function(x){
      t = x
      id = sum(eta_angle < t)
      if(id==0){
        res = t/(4*t1)
      }else if(id==1){
        res = t/(2*t2 - 2*t1) + (t2 - 3*t1)/(4*t2 - 4*t1)
      }else if(id==2){
        res = t/(2*t3 - 2*t2) + (3*t3 - 5*t2)/(4*t3 - 4*t2)
      }else if(id==3){
        res = t/(2*t4 - 2*t3) + (5*t4 - 7*t3)/(4*t4 - 4*t3)
      }else if(id==4){
        res = t/(8*pi - 4*t4) + (7*pi - 4*t4)/(4*pi - 2*t4)
      }
      return(res * pi)
    }
    if(length(x)>1){
      sapply(x,tau_entrywise)
    }else{
      tau_entrywise(x)
    }
  }
  
  # tau^{-1}:[0,2\pi] to [0,2\pi]
  tau_inverse <- function(x){
    t1 = eta_angle[1]; t2 = eta_angle[2]; t3 = eta_angle[3]; t4 = eta_angle[4]
    p1 = pi/4; p2 = 3*pi/4; p3 = 5*pi/4; p4 = 7*pi/4
    tau_inverse_entrywise <- function(x){
      t = x
      id = sum(c(p1,p2,p3,p4) < t)
      if(id==0){
        res = t1/(pi/4) * t
      }else if(id==1){
        res = t1 + (t-p1) * (t2-t1)/(pi/2)
      }else if(id==2){
        res = t2 + (t-p2) * (t3-t2)/(pi/2)
      }else if(id==3){
        res = t3 + (t-p3) * (t4-t3)/(pi/2)
      }else if(id==4){
        res = t4 + (t-p4) * (2*pi-t4)/(pi/4)
      }
      return(res)
    }
    if(length(x)>1){
      sapply(x,tau_inverse_entrywise)
    }else{
      tau_inverse_entrywise(x)
    }
  } 
  
  
  # R:D^2 \to D^2
  R <- function(Z){
    n = nrow(Z)
    if(is.null(n)){
      vectors(r = sqrt(sum(Z^2)), theta = tau( md(angles(Z) + theta) ))
    }else{
      rt = cbind(apply(Z,1,function(x) sqrt(sum(x^2))), tau( md(angles(Z) + theta) ) )
      t(sapply(1:n, function(i) vectors(r = rt[i,1], theta = rt[i,2])))
    }
  }
  
  # R^{-1}:D^2 \to D^2
  R_inverse <- function(Z){
    n = nrow(Z)
    if(is.null(n)){
      vectors(r = sqrt(sum(Z^2)), theta = md(tau_inverse(angles(Z)) - theta))
    }else{
      rt = cbind(apply(Z,1,function(x) sqrt(sum(x^2))), md(tau_inverse(angles(Z)) - theta) )
      t(sapply(1:n, function(i) vectors(r = rt[i,1], theta = rt[i,2])))
    }
  }
  
  # rho: I^2 \to I^2
  rho <- function(Z) omega_inverse(R(omega(Z)))
  
  # rho^{-1}: I^2 \to I^2
  rho_inverse <- function(Z) omega_inverse(R_inverse(omega(Z)))
  
  Z = rho(Y)
  
  
  ## [== Plot ==]
  
  # visualization of the homeomorphism \rho
  epsilon = 0.02 ## margin parameter for empirical level set
  for(k in 1:2){
    pdf(file=paste0("P1_rotation_k=",k,".pdf"), width=10, height=5)
    par(mfrow=c(1,2)); par(mar = c(4.2, 4.2, 2, 1)); par(oma = c(0,0,0,0))
    visualizer_demo(X=X, Y=Y, k=k, epsilon=epsilon, title=expression({1: L[y[1]](X)}))
    visualizer_demo(X=X, Y=Z, k=k, epsilon=epsilon, title=expression({2: L[rho(Y)[1]](X)}))
    dev.off()
  }
  
  save(file = "D2_Estimating_rho.RData",f1_vertices,zeta,theta,eta,tau,tau_inverse,R,R_inverse,rho,rho_inverse,Z)
  
  
  ## =============================================
  ## Estimating a second-step estimator f^2(x)
  ## =============================================
  
  ## [== Estimation ==]
  
  # parameters
  M = 100; n_test = M^2 # number of grid for visualization (not used for estimation)
  # t_list = c(1,3,5) # list of t, specifying a grid \hat{I} = {0, \pm 1/t, \pm 2/t, ..., \pm 1} (used for estimation)
  t_list = 2^(2:4)
  
  # New covariates for visualizing groundtruth/estimators
  X_test = expand.grid(seq(-1,1,length.out=sqrt(n_test)),seq(-1,1,length.out=M))
  G = expand.grid(1:M,1:M)
  
  save(file = "D3_settings.RData", M, n_test, t_list, X_test, G)
  
  
  # groundtruth (formed as matrix, to plot heatmap)
  groundtruth_M1 = groundtruth_M2 = matrix(0,M,M); groundtruth = t(apply(X_test,1,f))
  for(j in 1:(M^2)){
    groundtruth_M1[G[j,1],G[j,2]] = groundtruth[j,1]
    groundtruth_M2[G[j,1],G[j,2]] = groundtruth[j,2]
  }
  
  # first-step estimator f^1(x)
  f1_M1 = f1_M2 = matrix(0,M,M); f1 = kNN(X_obs = X, Y_obs = Y, X_points = X_test, k = 10, truncate = TRUE)
  for(j in 1:(M^2)){
    f1_M1[G[j,1],G[j,2]] = f1[j,1]
    f1_M2[G[j,1],G[j,2]] = f1[j,2]
  }
  
  # empirically rotated set g = rho(f^1(x)) (before smooth interpolation)
  g_M1 = g_M2 = matrix(0,M,M); g = rho(f1)
  for(j in 1:(M^2)){
    g_M1[G[j,1],G[j,2]] = g[j,1]
    g_M2[G[j,1],G[j,2]] = g[j,2]
  }
  
  save(file = "D3_Groundtruth.RData",groundtruth_M1, groundtruth_M2, f1_M1, f1_M2, g_M1, g_M2)
  
  
  for(t in t_list){
    div = 2 * t + 1 ## number of division (for each direction) in the grid
    I = seq(-1,1,length.out = div)
    delta = I[2]-I[1] ## interval between each grid
    
    # justify the boundaries of empirically rotated g 
    coherent_g <- array(0,dim=c(div,div,2))
    for(i in 1:div){
      for(j in 1:div){
        tmp <- rho(kNN(X_obs = X, Y_obs = Y, X_points = c(I[i],I[j]), k = 20, truncate = TRUE))
        if(i==1) tmp[1] = -1 else if(i==div) tmp[1] = 1
        if(j==1) tmp[2] = -1 else if(j==div) tmp[2] = 1
        coherent_g[i,j,] = tmp
      }
    }
    
    # smoothly interpolate the coherent_g)
    g_dagger_function <- function(x){
      i = sum(I <= x[1]); j = sum(I <= x[2])
      alpha = (x[1] - I[i])/delta; beta = (x[2] - I[j])/delta
      
      if(i == length(I) && j == length(I)){
        tmp = coherent_g[i,j,]
      }else if(i == length(I)){
        u3 = coherent_g[i,j,]; u4 = coherent_g[i,j+1,]
        tmp = beta * u4 + (1-beta) * u3
      }else if(j == length(I)){
        u2 = coherent_g[i+1,j,]; u3 = coherent_g[i,j,]
        tmp = alpha * u2 + (1-alpha) * u3
      }else{
        u1 = coherent_g[i+1,j+1,]; u2 = coherent_g[i+1,j,]; u3 = coherent_g[i,j,]; u4 = coherent_g[i,j+1,]
        v1 = beta * u1 + (1-beta) * u2
        v2 = beta * u4 + (1-beta) * u3
        tmp = alpha * v1 + (1-alpha) * v2
      }
      return(tmp)
    }
    g_dagger = t(apply(X_test,1,g_dagger_function))
    
    # second step estimator (proposed)
    second_step <- function(x) rho_inverse(g_dagger_function(x))
    f2 = t(apply(X_test,1,second_step))
    
    g_dagger_M1 = g_dagger_M2 = f2_M1 = f2_M2 = matrix(0,M,M)
    
    for(j in 1:(M^2)){
      g_dagger_M1[G[j,1],G[j,2]] = g_dagger[j,1]
      g_dagger_M2[G[j,1],G[j,2]] = g_dagger[j,2]
      
      f2_M1[G[j,1],G[j,2]] = f2[j,1]
      f2_M2[G[j,1],G[j,2]] = f2[j,2]  
    }
    
    save(file = paste0("D3_second_step_t=",t,".RData"), t, div, I, delta, coherent_g, g_dagger_function, second_step, 
         g_dagger, f2, g_dagger_M1, g_dagger_M2, f2_M1, f2_M2, 
         coherent_g)
    
  }
  
  
  ## [== Plot ==]
  
  axis_labels <- function(){
    axis(side=1, at=seq(0,1,0.25), labels=seq(-1,1,0.5))
    axis(side=2, at=seq(0,1,0.25), labels=seq(-1,1,0.5))
    return(NULL)
  }
  
  for(k in 1:2){
    pdf(file=paste0("P2_effect_of_smoothness_k=",k,".pdf"), width=5*(1+length(t_list)), height=5)
    par(mfrow=c(1,1+length(t_list)),cex=1.5,cex.axis=1.3); par(mar = c(4.2, 4.2, 2, 1)); par(oma = c(0,0,0,0))
    if(k==1){
      image(groundtruth_M1, main=expression(Groundtruth: {f['*,1']}), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
      for(t in t_list){
        load(file=paste0("D3_second_step_t=",t,".RData"))
        image(f2_M1, main=bquote(Estimator: {hat(f)}['n,1'] ~ '( t' == .(t) ~ ')'), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
      }
    }else{
      image(groundtruth_M2, main=expression(Groundtruth: {f['*,2']}), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
      for(t in t_list){
        load(file=paste0("D3_second_step_t=",t,".RData"))
        image(f2_M2, main=bquote(Estimator: {hat(f)}['n,2'] ~ '( t' == .(t) ~ ')'), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
      }
    }
    dev.off()  
  }
  
  
  pdf(file=paste0("P3_1_groundtruth.pdf"), width=5, height=5)
  par(cex=1.5,cex.axis=1.3); par(mar = c(2, 2, 2, 1)); par(oma = c(0,0,0,0))
  image(groundtruth_M1, main=expression(Groundtruth: {f['*,1']}), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
  dev.off()
  
  pdf(file=paste0("P3_2_first_step.pdf"), width=5, height=5)
  par(cex=1.5,cex.axis=1.3); par(mar = c(2, 2, 2, 1)); par(oma = c(0,0,0,0))
  image(f1_M1, main=expression('first step:'~ {hat(f)['n,1']^{(1)}}), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
  dev.off()
  
  pdf(file=paste0("P3_3_g.pdf"), width=5, height=5)
  par(cex=1.5,cex.axis=1.3); par(mar = c(2, 2, 2, 1)); par(oma = c(0,0,0,0))
  image(g_M1, main=expression({hat(g)['n,1']}), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
  dev.off()
  
  for(t in t_list){
    
    load(file=paste0("D3_second_step_t=",t,".RData"))
    
    pdf(file=paste0("P3_4_g_dagger_t=",t,".pdf"), width=5, height=5)
    par(cex=1.5,cex.axis=1.3); par(mar = c(2, 2, 2, 1)); par(oma = c(0,0,0,0))
    image(g_dagger_M1, main=bquote({hat(g)['n,1']^{"\u2020"}} ~ '( t' == .(t) ~ ')'), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
    dev.off()
    
    pdf(file=paste0("P3_5_second_step_t=",t,".pdf"), width=5, height=5)
    par(cex=1.5,cex.axis=1.3); par(mar = c(2, 2, 2, 1)); par(oma = c(0,0,0,0))
    image(f2_M1, main=bquote({hat(f)['n,1']^{(2)}} ~ '( t' == .(t) ~ ')'), zlim=c(-1.1,1.1), xaxt="n", yaxt="n"); axis_labels()
    dev.off()
    
  }
  
  
  pdf(file="P4_Lg_grid.pdf", height=5, width=5 * length(t_list))
  par(mfrow=c(1,length(t_list)), cex=1.5, cex.axis=1); 
  par(mar = c(4.2, 4.2, 2, 1)); par(oma = c(0,0,0,0))
  
  for(t in t_list){
    
    load(file=paste0("D3_second_step_t=",t,".RData"))
    
    plot(0,0,type="n",xlim=c(-1,1),ylim=c(-1,1),xlab=expression(x[1]),ylab=expression(x[2]),
         main=paste0("t=",t))
    
    draw_seg <- function(v1,v2) segments(x0 = v1[1], y0 = v1[2], 
                                         x1 = v2[1], y1 = v2[2])
    
    for(i in 1:(div-1)){
      for(j in 1:(div-1)){
        draw_seg(coherent_g[i,j,], coherent_g[i+1,j,])
        draw_seg(coherent_g[i,j,], coherent_g[i,j+1,])
      }
    }
    for(k in 1:(div-1)){
      draw_seg(coherent_g[div,k,],coherent_g[div,k+1,])
      draw_seg(coherent_g[k,div,],coherent_g[k+1,div,])
    }
    
  }
  dev.off()
}

setwd(DIR_NAME)
