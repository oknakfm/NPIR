# Overview
This repository provides R Source codes to reproduce numerical experiments in the following paper:

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

## <a href="https://github.com/oknakfm/NPIR/blob/main/1_NPIR.R">1_NPIR.R</a>
Computes nonparametric planer invertible regression. The results are saved into (automatically-generated) output folder. 

### This script estimates the invertible function by the following procedure:
- For I=[-1,1], generate x1,x2,...,xn uniformly randomly from I^2, and also generate y1,y2,...,yn \in I^2 fandomly from a normal distribution N(f(xi),vI). These (x_i,y_i) are used for experiments.
- Compute a first step estimator f^{(1)}(x) using k-nearest neighbour
- Empirically estimate the homeomorphism \rho (via estimation of the rotation of the four-vertices of I)
- After estimating g = \rho \circ f^{(1)}, define coherent_g which slightly fixes the four vertices
- Using the coherent_g, define g_dagger which smoothly interpolates the internal region (of each triangle connecting vertices)
- Compute the proposed estimator f^{(2)} = \rho^{-1} \circ g_\dagger

Experimental results for a function (whose level set includes twists) are saved into OUTPUT/TWIST. 
Those for a function whose level set does not include twists (i.e., its estimation is easier) are saved into OUTPUT/NON_TWIST. 

### Output 1 (P1): Rotation that matches the four vertices of I^2
<img src="/ForReadMe/1_rotation.png" width="400">

### Output 2 (P2): Influence of the precision parameter t to the estimation of the function f
<img src="/ForReadMe/2_effect_of_t.png" width="800">

### Output 3 (P3): Visualization of each step of the estimation
<img src="/ForReadMe/3_1_groundtruth.png" width="200"> <img src="/ForReadMe/3_3_g.png" width="200"> <img src="/ForReadMe/3_5_second_step_t=3.png" width="200">

### Output 4 (P4): Grid used in the function estimation
<img src="/ForReadMe/4_Lg_grid.png" width="600">

## <a href="https://github.com/oknakfm/NPIR/blob/main/2_true_levelsets.R">2_TrueLevelSets.R</a>

### Output 5 (P5): Visualization of the level-set for the underlying true function f*
<img src="/ForReadMe/5_levelset_of_f.png" width=800>

## <a href="https://github.com/oknakfm/NPIR/blob/main/A_functions.R">A_functions.R</a>
This script contains several functions to be called from the above scripts.


