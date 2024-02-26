## Installation
Make sure to have the following dependencies installed on your system:

- a C++17 compliant compiler
- the  following R packages: `Rcpp`, `RcppEigen`, `rgl`, `Matrix`, `plot3D`, `plot3Drgl`, `shiny`, `MASS`, `testthat`.

1.  use the `devtools` package. From the R console, execute

<!-- -->

    devtools::install_github("aldoclemente/fdaPDEISCHIA", ref="main") 

2.  clone this repository and install. From a terminal, execute

<!-- -->

    git clone -b main git@github.com:aldoclemente/fdaPDEISCHIA.git 
    cd path/to/fdaPDEISCHIA 

and install the package from the R console

    install.packages(".", type="source", repos=NULL) 
    
    
