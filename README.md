# RANN2
<!-- badges: start -->
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://jefferis.github.io/RANN2/)
[![Build Status](https://travis-ci.org/jefferis/RANN2.svg)](https://travis-ci.org/jefferis/RANN2)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

This package is an updated version of the [RANN](https://cran.r-project.org/package=RANN) 
package, making use of the Rcpp package. 
For basic use, there is little difference with original RANN 
package although there are some small (typically 5-10%) speedups for certain 
query/target size combinations. **RANN2** also includes experimental 
functionality via `WANN` objects to:
  * keep ANN points in memory to avoid repeated copying
  * keep the ANN k-d tree in memory to avoid repeated building
  * separate building the k-d tree from allocating the points
  * permit very fast self queries
  * permit queries of the points from one ANN tree against a second tree

## Installation
Currently there isn't a released version on [CRAN](https://cran.r-project.org/),
although we are considering a submission when the package develops sufficiently
distinct functionality from the original RANN package.

### Development version
You can use the **devtools** package to install the development version:

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("jefferis/RANN2")
```

Note: Windows users need [Rtools](http://www.murdoch-sutherland.com/Rtools/) and [devtools](https://cran.r-project.org/package=devtools) to install this way.

## Use
### Basic use
The expectation is that for 90% of users the `nn2` function should be the only
way that the library is used. This takes a target matrix of R points, copies them into
an array used by ANN and builds a k-d tree. It then iterates over the query
points, searching the tree one at a time.

### Advanced use
**RANN2** adds `WANN` objects, which allow fine control of when the k-d tree is
built and removed.

```
data(kcpoints)
w1=WANN(kcpoints[[1]])
library(microbenchmark)
microbenchmark(w1sq<-w1$selfQuery(k=1,eps=0))
microbenchmark(nn2(kcpoints[[1]],k=1))
w2=WANN(kcpoints[[2]])
# NB must pass the Cpp object not the reference class object
w1$queryWANN(w2$.CppObject)
```

WANN objects will primarily be useful if you make repeated queries. You can also
delay building the k-d tree:

```
w1=WANN(kcpoints[[1]])
w1$querySelf(k=1,eps=0)
w1$build_tree()
w1$delete_tree()
```
if only a fraction of the objects will need to be searched; the tree will
automatically be built when it is queried. You can also explicitly control
when the tree is built or deleted (for memory management). The tree is wrapped
in an R reference class (R5) object which imposes a significant performance
penalty for building small trees (< ~ 1000 points).

### Changing ANN data types
By default ANN uses `double`s for both points and returned distances. You can
save space by changing this if you want. To do to this you must recompile after
setting either `ANN_COORD_TYPE` or `ANN_DIST_TYPE` in `src/MAKEVARS` or 
`MAKEVARS.win` as appropriate. e.g. 
```
PKG_CPPFLAGS=-I. -IANN -DRANN -DANN_COORD_TYPE=float
```
would switch to the use of floats for the main ANN coordinate type. Note however
that the k-d tree itself appears to occupy ~ 2x the space of the underlying
double coordinates.

### Linking and using ANN library
This package compiles the static library for the ANN library and provides the headers for it, so developers can directly include them in their cpp code. 

#### Instructions

1. `DESCRIPTION` file: 
    ```
    LinkingTo: RANN2
    ```
1. `src/Makevars` file: 
    ```R
    PKG_IMPORT=RANN2
    PKG_HOME=`${R_HOME}/bin/Rscript -e 'cat(system.file(package=\"$(PKG_IMPORT)\"))'`
    PKG_LIBS=-L$(PKG_HOME)/lib -l$(PKG_IMPORT)
    ```
1. `src/Makevars.win` file:
    ```R
    PKG_IMPORT=RANN2
    PKG_HOME=`${R_HOME}/bin/Rscript -e 'cat(system.file(package=\"$(PKG_IMPORT)\"))'`
    PKG_LIBS+=-L$(PKG_HOME)/lib -l$(PKG_IMPORT)

    PKG_CPPFLAGS+=-DDLL_EXPORTS
    ```
1. Your `cpp` file: 
    ```C
    #include <ANN.h>
    ```
    
For usage example: [src/ann.cpp](https://github.com/jefferis/RANN2/blob/master/src/nn.cpp)
