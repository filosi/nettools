nettools

R package for network analysis

branches:
master 		-> version freeze to 0.9.4
dev 		-> development version for adding network stability functions merged into master

## Download the zip file

## Move to the download directory
cd ~/Dowloads

## Unpack it with:
unzip nettools-dev.zip

## Build the binary for R
R CMD build nettools-dev

## Install it using (change the version according to the latest version):
R CMD INSTALL nettools_0.9.5.tar.gz
