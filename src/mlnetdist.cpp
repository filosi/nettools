/*  
	nettools R package
    
    This code is written by Michele Filosi <michele.filosi@gmail.com>.
    Copyright (C) 2012 Michele Filosi, Copyright (C) 2012 Fondazione
    Bruno Kessler.

    References: 
	a) A. Sole-Ribalta, M. De Domenico, N.E. Kouvaris, A. Diaz-Guuilera, 
	   S. Gomez, A. Arenas.
	   Spectral properties of the Laplacian of multiplex networks.
	   arXiv: 1307.2090v1
	
    b) M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman, C. Furlanello
       Stability Indicators in Network Reconstruction, PLOSONE 2014
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
#include "mlnetdist.h"

using namespace Rcpp;


/* Read the array */
array3d *readarray(double *mod, int *dim, int n){
  
  array3d *myarray;
  int i, j, s, offset;
  
  // Initialize array3d structure
  myarray = (array3d *) malloc (sizeof(array3d));

  // Array dimensions
  myarray->dim = (int *) malloc (n * sizeof(int));
  
  // Initialize array dimension
  for (i=0; i<n; i++){
	myarray->dim[i] = dim[(n - (i+1))];
  }
  
  // Allocate space for the numeric array
  myarray->A = (double ***) malloc (myarray->dim[0] * sizeof(double **));
  
  // Initialize the array with the input values
  for (i=0; i<myarray->dim[0]; i++){
  	myarray->A[i] = (double **) malloc (myarray->dim[1] * sizeof(double *));
  	offset = myarray->dim[1] * myarray->dim[2] * i;

  	for (j=0; j<myarray->dim[1]; j++){
  	  myarray->A[i][j] = (double *) malloc (myarray->dim[2] * sizeof(double));
	  
  	  for (s=0; s<myarray->dim[2]; s++){
  		myarray->A[i][j][s] = mod[(offset + (myarray->dim[1] * j) + s)];
  	  }
  	}
  }
  return myarray;
}


/* Compute the degrees of each matrix in the array3d. Namely for the
   matrix defined by slicing the 3 dimension array3d -> dim[2] 

   Returns a bi-dimensional array
*/
double **arrdegree(array3d *myarray){
  
  double **degree;
  double mysum;
  int i, j, s;
  
  /* Alloc degree vector */
  degree = (double **) malloc (myarray->dim[0]*sizeof(double*));
  for (s=0; s<myarray->dim[0]; s++){
	degree[s] = (double *) malloc (myarray->dim[1]*sizeof(double));
  }

  /* Compute degree by layer */
  for (s=0; s<myarray->dim[0]; s++){
	for (i=0; i<myarray->dim[1]; i++){
	  mysum = 0.0;
	  for (j=0; j<myarray->dim[2]; j++){
		mysum = mysum + myarray->A[s][i][j];
	  }
	  degree[s][i] = mysum;
	}
  }
  return degree;
}


/* 
  This function copy the original array and 
  compute the laplacian for each adjacency matrix.
*/
void arrLap (array3d **myarray){
  
  int i, j, s;
  double **mydegree;
  array3d *myarrtmp = *myarray;
  
  /* Compute the degree array for each layer */
  mydegree = arrdegree(myarrtmp);

  for (i=0; i<myarrtmp->dim[0]; i++){
  	for (j=0; j<myarrtmp->dim[1]; j++){
  	  for (s=0; s<myarrtmp->dim[2]; s++){
		if (j == s){
		  /* Assign degree onto the major diagonal */
		  myarrtmp->A[i][j][s] = mydegree[i][s];
		} else {
		  /* Fill the off-diagonal positions of the Laplacian by slide */
		  myarrtmp->A[i][j][s] = (-1 * myarrtmp->A[i][j][s]);
		}
  	  }
  	}
  }

  // Free degree array
  for (i=0; i<myarrtmp->dim[0]; i++){
	free(mydegree[i]);
  }
  free(mydegree);
}

/*
  This function computes the intra layer laplacian matrix defined by:
  LL = L_1    0    0    0    ....
        0     L_2  0    0    ....
		.     
		.        0    .....
		.
		0       .....         L_l
  
  where l = 1, ..., #layers and L_l is the Laplacian of the network on the layer l
  
  Inputs: 
  NumericVector -> a R array mapped as numeric vector in Rcpp.
  A multilayer network is represented by a 3D multidimensional
  array. The firsst 2 dimensions are the dimension of the network on
  each layer. The third dimension the layer dimension.

  Outputs:
  NumericMatrix -> a NumericMatrix of n x p  | n = 1, ..., #nodes and p=1, ..., #layers
*/


// [[Rcpp::export]]
NumericMatrix intraLaplacian(NumericVector mod) {
  IntegerVector mydim;
  mydim = as<IntegerVector>(mod.attr("dim"));

  array3d *myarray;

  int i, j, p, n, s, x, y;
  p = mydim[0];
  n = p * mydim[2];
 
  /* Initialize the matrix */
  NumericMatrix LL(n, n);
  for (i=0; i<n; i++){
  	for (j=i+1; j<n; j++){
  	  LL(i,j) = 0.0;
  	  LL(j,i) = 0.0;
  	}
  }
  
  /* Read the array from R */
  myarray = readarray(REAL(mod), INTEGER(mydim), 3);

  /* Get the Laplacian from the array by layer */
  arrLap(&myarray);
  
  /* Fill the matrix with inter-layer Laplacian */
  for (s=0; s<myarray->dim[0]; s++){
  	for (i=0; i<myarray->dim[1]; i++){
  	  y = (p * s) + i;
  	  for (j=i; j<myarray->dim[2]; j++){
  		x = (p * s) + j;
  		LL(x,y) = myarray->A[s][i][j];
  		LL(y,x) = myarray->A[s][i][j];
  	  }
  	}
  }

  /* Free multi-slice array */
  for (s=0; s<myarray->dim[0]; s++){
  	for (i=0; i<myarray->dim[1]; i++){
  	  free(myarray->A[s][i]);
  	}
  	free(myarray->A[s]);
  }
  free(myarray->A);
  free(myarray->dim);
  free(myarray);
  
  /* Return the intra-Laplacian */
  return LL;
}

/* 
   This function computes the Kroenecker direct product between a
   matrix Li and the Identity Matrix of dimension n 
   
   Inputs: 
   Matrix Li -> Laplacian of interlayer network defined by Si - Wi where Si is
   the degree matrix and Wi is the adjacency matrix defining the
   connection between layers of the network.
   
   Int n     -> Number of nodes of a layer network
   
   Outputs:
   NumericMatrix -> a numerical matrix of n x #layers composed as:
   Li KroenekercProd I_(n,n)
*/

// [[Rcpp::export]]
NumericMatrix directProd(NumericMatrix Li, int n) {
  IntegerVector mydim;
  mydim = as<IntegerVector>(Li.attr("dim"));
  
  int i, j, p, x, y;
  double tmp;
  
  /* Number of layers */
  p = mydim[0];

  /* Alloc the space to store the inter-Laplacian */
  NumericMatrix LI(n*p,n*p);

  /* Initialize and fill the matrix with zeros */
  for (i=0; i<(n*p); i++){
  	for (j=i+1; j<(n*p); j++){
  	  LI(i,j) = 0.0;
  	  LI(j,i) = 0.0;
  	}
  }
  
  /* Fill the matrix with correct values */
  for (x=0; x<p; x++){
	for (y=0; y<p; y++){
	  tmp = Li(x,y);
	  for (i=0; i<n; i++){
		LI((n*x) + i,(n*y) + i) = tmp;
	  }
	}
  }
  
  /* Return inter-Laplacian matrix */
  return LI;
}

// [[Rcpp::export]]
NumericMatrix IntraAdj(NumericVector mod) {
  IntegerVector mydim;
  mydim = as<IntegerVector>(mod.attr("dim"));
  int i, j, l, n, s, x, y, offset;
  double tmp;
  
  l = mydim[2];
  n = mydim[0];
 
  /* Initialize the matrix nodes by layer */
  NumericMatrix LL(n*l, n*l);
  for (i=0; i<(n*l); i++){
  	for (j=i+1; j<(n*l); j++){
  	  LL(i,j) = 0.0;
  	  LL(j,i) = 0.0;
  	}
  }

  /* Fill the matrix with inter-layer Laplacian */
  for (s=0; s<l; s++){
  	offset = n * n * s;
	for (i=0; i<n; i++){
	  y = (n * s) + i;
	  for (j=i; j<n; j++){
  		x = (n * s) + j;
  		tmp = mod[(offset + (n * i) + j)];
  		LL(x,y) = tmp;
  		LL(y,x) = tmp;
  	  }
  	}
  }

  /* Return the intra-Laplacian */
  return LL;
}

