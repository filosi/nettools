# nettools
This package is available on [CRAN](https://cran.r-project.org/web/packages/nettools/index.html).

## R package for differential network analysis
branches:
* master 		-> version freeze to 1.0.4
* dev 		-> development version only unstable and new functions


## Installation instructions

### Package Dependecy
**From CRAN**
```
install.packages(c("infotheo", "rootSolve", "dtw", "minerva", "combinat", "WGCNA", "R.utils"))
```

**From bioconductor**
```
source("https://bioconductor.org/biocLite.R")
biocLite(c("impute", "minet", "AnnotationDbi", "GO.db"))
```

**Optional**
```
install.packages(c("igraph", "Matrix", "Rcpp"))
```

### Option 1 

**Download the zip file from github**
```
https://github.com/filosi/nettools/archive/master.zip
```
**Move to the download directory**
``` 
cd ~/Dowloads
```

**Unpack it with:**
``` 
unzip nettools-master.zip
```

**Build the binary for R**
``` 
R CMD build nettools-master
```

**Install it using (change the version according to the latest version):**
``` 
R CMD INSTALL nettools_X.Y.Z.tar.gz
```

### Option 2

Directly from R with the devtools library (https://github.com/hadley/devtools) installed 

``` R
install.packages("devtools")
library(devtools)
install_github("filosi/nettools")
```
