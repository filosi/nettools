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


/* 3D array structure to store multilayer network */
typedef struct array3d
{
  int *dim;       // Array dimension of length 3
  double ***A;    // data Array[dim[2]][dim[1]][dim[0]]
} array3d;
