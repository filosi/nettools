# nettools

## R package for network analysis

branches:
* master 		-> version freeze to 0.9.4
* dev 		-> development version for adding network stability functions merged into master

## Installation instructions

### Option 1 

**Download the zip file**

**Move to the download directory**
``` 
cd ~/Dowloads
```

**Unpack it with:**
``` 
unzip nettools-dev.zip
```

**Build the binary for R**
``` CMD build nettools-dev
R CMD build nettoo-dev

**Install it using (change the version according to the latest version):**
``` 
R CMD INSTALL nettools_1.0.4.tar.gz
```

### Option 2

Directly from R with the devtools library (https://github.com/hadley/devtools) installed 

``` R
library(devtools)
install_github("nettools","MPBA")
```
