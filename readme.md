# GRFics
![](https://i.imgur.com/1Ptm5NI.png)

GRFics is an R library for efficient generation of approximate realizations from two-dimensional Gaussian random fields (GRFs).

The package aims to be simple to use, even for those without a background in statistics. 

## Technical details
The approach utilized was first described in Lindgren et al. (2011), where a stochastic partial differential equation (SPDE) is discretized to obtain a Gaussian Markov random field (GMRF). 
This package implements the variation of the approximation derived in  Fuglstad et al. (2015), where the SPDE is defined on a rectangular region and discretized on a regular grid. A periodic boundary condition is used, leading to realizations that are periodic along the vertical and horizontal boundaries.

## Installation
The package can be installed by using the function `install_github` from the `devtools` package:

``` r
devtools::install_github("mathiasisaksen/GRFics")
```

## Usage
Each function contains documentation with examples of usage. In RStudio, this documentation can be accessed by writing, for example, `?generate.grf.object`.

## References
Finn Lindgren, Håvard Rue, and Johan Lindström. "An explicit link between Gaussian fields and Gaussian markov random fields: the stochastic partial differential equation approach." *Journal of the Royal Statistical Society: Series B
(Statistical Methodology)*, 73(4):423–498, 2011.

Geir-Arne Fuglstad, Finn Lindgren, Daniel Simpson, and Håvard Rue. "Exploring a new class of non-stationary spatial Gaussian random fields with varying local anisotropy." *Statistica Sinica*, pages 115–133, 2015.
