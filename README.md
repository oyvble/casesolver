
# CaseSolver

Software for comparing DNA profiles within cases (Expert system)

## Requiring EuroForMix and EFMex

It is recommended to have euroformix R-package installed (at least v4.0)
before installing CaseSolver. To enable the EFMex option you need to
have the EFMex R-package installed (at least v0.8.1).

## Installation of other required packages

``` r
install.packages(c('gWidgets2tcltk','readxl')) #required
install.packages(c('plotly','igraph','plotrix','R2HTML','officer','flextable')) #useful
```

## Installation of CaseSolver

Installation directly from source through GitHub (does not require
compilation):

``` r
install.packages("remotes")
remotes::install_github("oyvble/casesolver")
```

## Run CaseSolver

``` r
library(casesolver);gui() 
```
