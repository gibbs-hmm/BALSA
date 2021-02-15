/****************************************************************/
/* BALSA - Bayesian Algorithm for Local Sequence Alignment      */
/*                                                              */
/* Please acknowledge the program authors on any publication of */
/* scientific results based in part on use of the program and   */
/* cite the following articles in which the program was         */
/* described.                                                   */
/*                                                              */
/* Webb-Robertson, B. J., L. A. McCue, et al. (2008). "Measuring*/
/* Global Credibility with Application to Local Sequence        */
/* Alignment." PLoS CompBiol. 4(5): e1000077.                   */
/* Webb, B. J., J. S. Liu, et al. (2002). "BALSA: Bayesian      */
/* algorithm for local sequence alignment." Nucleic Acids Res   */
/* 30(5): 1268-77.                                              */
/*                                                              */
/*                                                              */
/* Copyright (C) 2006   Health Research Inc.                    */
/* HEALTH RESEARCH INCORPORATED (HRI),                          */
/* ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.             */
/* Email:  gibbsamp@wadsworth.org                              */
/*                                                              */
/*                                                              */
/* Copyright (C) 2009   Brown University                        */
/* Brown University                                             */
/* Providence, RI 02912                                         */
/* Email:  gibbs@brown.edu                                     */
/****************************************************************/
/*                                                              */
/* This program is free software; you can redistribute it       */
/* and/or modify it under the terms of the GNU General Public   */
/* License as published by the Free Software Foundation;        */
/* either version 2 of the License, or (at your option)         */
/* any later version.                                           */
/*                                                              */
/* This program is distributed in the hope that it will be      */
/* useful, but WITHOUT ANY WARRANTY; without even the implied   */
/* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR      */
/* PURPOSE. See the GNU General Public License for more         */
/* details.                                                     */
/*                                                              */
/* You should have received a copy of the GNU General Public    */
/* License along with this program; if not, write to the        */
/* Free Software Foundation, Inc., 51 Franklin Street,          */
/* Fifth Floor, Boston, MA  02110-1301, USA.                    */
/****************************************************************/

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include "matrix.h"
#include "bay_matrix_all.h"
#include "BALSA_file_score_sums.h"

using namespace std;

//Time constructor
Score_Sums::Score_Sums() {

}

// Finds mode of gap opening and gap extension parameters
void Score_Sums::ScoreMatrix(int MAX1, int MAX2, int *seq1index, int *seq2index, int matrixnum)
{
  Matrix p(MAX1+1,MAX2+1);
  for(i = 1; i <= MAX1; ++i)
    p(i,0) = 0;
  for(j = 0; j <= MAX2; ++j)
    p(0,j) = 0;

  // Protien Matrix 
  for(i = 1; i <= MAX1; ++i) {   
    mat1 = seq1index[i];
    for(j = 1; j <= MAX2; ++j) {     
      mat2 = seq2index[j];
      p(i,j) = p_related[mat[matrixnum]][mat1][mat2];    
    }  
  }  
 
  gapo = gap_open[matrixnum];
  gape = gap_exte[matrixnum];
  numer = 0;
  denom = 0;
  Matrix MatM(MAX1+1, MAX2+1);
  Matrix MatI(MAX1+1, MAX2+1);
  Matrix MatD(MAX1+1, MAX2+1);
  Matrix MatN(MAX1+1, MAX2+1);
  Matrix MatE(MAX1+1, MAX2+1);
  
  for(row = 0; row <= MAX1; ++row) {
    MatM(row,0) = 0;
    MatI(row,0) = 0;
    MatD(row,0) = 0;   
    MatN(row,0) = 0;
    MatE(row,0) = 0;
  }
  for(col = 1; col <= MAX2; ++col) {
    MatM(0,col) = 0;
    MatI(0,col) = 0;
    MatD(0,col) = 0;
    MatN(0,col) = 0;
    MatE(0,col) = 0;
  }
  
  for(row = 1; row <= MAX1; ++row) {
    for(col = 1; col <= MAX2; ++col) {
      MatM(row,col) = (MatM(row-1,col-1) + MatI(row-1,col-1) + MatD(row-1,col-1) +
		       MatN(row-1,col-1)) * p(row,col);
      MatI(row,col) = gape*MatI(row-1,col) + gapo*(MatM(row-1,col) + MatN(row-1,col));
      MatD(row,col) = gape*MatD(row,col-1) + gapo*(MatM(row,col-1) + MatN(row,col-1));
      MatN(row,col) = p(row,col);   
      MatE(row,col) = MatM(row,col) + MatI(row,col) + MatD(row,col) + MatN(row,col);
      numer += MatE(row,col);     
    }
  }
  
  
  for(row = 0; row <= MAX1; ++row) {
    MatM(row,0) = 0;
    MatI(row,0) = 0;
    MatD(row,0) = 0;   
    MatN(row,0) = 0;
    MatE(row,0) = 0;
  }
  for(col = 1; col <= MAX2; ++col) {
    MatM(0,col) = 0;
    MatI(0,col) = 0;
    MatD(0,col) = 0;
    MatN(0,col) = 0;
    MatE(0,col) = 0;
  }
  
  for(row = 1; row <= MAX1; ++row) {
    for(col = 1; col <= MAX2; ++col) {
      MatM(row,col) = MatM(row-1,col-1) + MatI(row-1,col-1) + MatD(row-1,col-1) +
	MatN(row-1,col-1);
      MatI(row,col) = gape*MatI(row-1,col) + gapo*(MatM(row-1,col) + MatN(row-1,col));
      MatD(row,col) = gape*MatD(row,col-1) + gapo*(MatM(row,col-1) + MatN(row,col-1));
      MatN(row,col) = 1;   
      MatE(row,col) = MatM(row,col) + MatI(row,col) + MatD(row,col) + MatN(row,col);
      denom += MatE(row,col);
    }
  }
  
  PR1R2_theta = numer/denom;
}


