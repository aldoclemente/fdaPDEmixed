## Installation
Make sure to have the following dependencies installed on your system:

- a C++17 compliant compiler
- the  following R packages: `Rcpp`, `RcppEigen`, `rgl`, `Matrix`, `plot3D`, `plot3Drgl`, `shiny`, `MASS`, `testthat`.

then, to install the latest stable version of the package, you can either:

1.  use the `devtools` package. From the R console, execute

<!-- -->

    devtools::install_github("aldoclemente/fdaPDEmixed", ref="main") 

2.  clone this repository and install. From a terminal, execute

<!-- -->

    git clone --recurse-submodules git@github.com:aldoclemente/fdaPDEmixed.git 
    cd fdaPDEmixed/ 

and install the package from the R console

    install.packages(".", type="source", repos=NULL) 
    
    
