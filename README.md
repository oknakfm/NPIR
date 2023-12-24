# Overview
This repository provides R Source codes to reproduce numerical experiments in the following arXiv preprint:

```
@article{okuno2023NPIR,
    year      = {2023},
    publisher = {},
    volume    = {},
    number    = {},
    pages     = {},
    author    = {Akifumi Okuno and Masaaki Imaizumi},
    title     = {Minimax Analysis for Inverse Risk in Nonparametric Planer Invertible Regression},
    journal   = {Electronic Journal of Statistics},
    note      = {accepted}
}
```

## 1_NPIR.R
Computes nonparametric planer invertible regression. The results are saved into (automatically-generated) output folder. 

### This script estimates the invertible function by the following procedure：
- For I=[-1,1], generate x1,x2,...,xn uniformly randomly from I^2
- Generate y1,y2,...,yn \in I^2 fandomly from a normal distribution N(f(xi),vI)
- Compute a first step estimator f^{(1)}(x) using k-nearest neighbour
- Empirically estimate the homeomorphism \rho (via estimation of the rotation of the four-vertices of I)
- After estimating g = \rho \circ f^{(1)}, define coherent_g which aligns four endpoints
- Using the coherent_g, define g_dagger which smoothly interpolate internal region of each triangles connecting vertices
- Compute the proposed estimator f^{(2)} = \rho^{-1} \circ g_\dagger

Experimental results for a function (whose level set includes twists) are saved into OUTPUT/TWIST. 
Those for a function whose level set does not include twists (i.e., its estimation is easier) are saved into OUTPUT/NON_TWIST. 

### Output 1 (P1): Rotation that matches the four vertices of I^2
<img src="/ForReadMe/1_rotation.png" width="400">

### Output 2 (P2): Influence of the precision parameter t to the estimation of the function f
おおよそうまく推定できている  
<img src="/ForReadMe/2_effect_of_t.png" width="800">

### 出力3 (P3): Visualization of each step of the estimation
<img src="/ForReadMe/3_1_groundtruth.png" width="200"> <img src="/ForReadMe/3_3_g.png" width="200"> <img src="/ForReadMe/3_5_second_step_t=3.png" width="200">

### 出力4 (P4): Grid used in the function estimation
<img src="/ForReadMe/4_Lg_grid.png" width="600">

## 2_TrueLevelSets.R

### 出力5 (P5): Visualization of the level-set for the underlying true function f*
<img src="/ForReadMe/5_levelset_of_f.png" width=800>

## A_functions.R
This script contains several functions to be called from the above scripts.


