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

t_list = 2^(2:5)

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
  
  pdf("P5_levelset_of_f.pdf",height=5,width=5*length(t_list))
  par(mfrow=c(1,length(t_list)),cex=1.5); par(mar = c(4.2, 4.2, 2, 1)); par(oma = c(0,0,0,0))
  
  for(t in t_list){
    
    z1_grid = seq(-1,1,1/t)
    z2_grid = seq(-1,1,length.out=1000)
    
    f1_M = f2_M = array(0,dim=c(length(z1_grid),length(z2_grid),2))
    for(i in 1:length(z1_grid)){
      for(j in 1:length(z2_grid)){
        f1_M[i,j,] = f(c(z1_grid[i],z2_grid[j]))
        f2_M[i,j,] = f(c(z2_grid[j],z1_grid[i]))
      }
    }
    
    plot(0,type="n",xlim=c(-1,1),ylim=c(-1,1),
         xlab=expression(x[1]),ylab=expression(x[2]), main=paste0("t=",t))
    draw_seg <- function(v1,v2,col="black") segments(x0 = v1[1], y0 = v1[2], 
                                                     x1 = v2[1], y1 = v2[2], col=col, lwd=0.5)
    for(i in 1:length(z1_grid)){
      for(j in 1:(length(z2_grid)-1)){
        draw_seg(f1_M[i,j,],f1_M[i,j+1,],col="red")
        draw_seg(f2_M[i,j,],f2_M[i,j+1,],col="blue")
      }
    }
    
  }
  
  dev.off()
  
}

setwd(DIR_NAME)