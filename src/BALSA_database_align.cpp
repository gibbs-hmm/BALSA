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
#include <string.h>
#include "matrix.h"
#include "bay_matrix_all.h"
#include "BALSA_file_align.h"
#include "BALSA_file_score_sums.h"

using namespace std;

//Time constructor
Align::Align() {

}

// Finds mode of gap opening and gap extension parameters
int* Align::Sample_Alignment(int MAX1, int MAX2, int *seq1index, int *seq2index, double *matscore)
{
  Matrix p(MAX1+1,MAX2+1);
  Matrix MatM(MAX1+1, MAX2+1);
  Matrix MatI(MAX1+1, MAX2+1);
  Matrix MatD(MAX1+1, MAX2+1);
  Matrix MatN(MAX1+1, MAX2+1);
  Matrix MatE(MAX1+1, MAX2+1);  
  Matrix Sample(MAX1+1, MAX2+1);
  sci = new int[MAX1+1];
  for(i = 1; i <= MAX1; ++i)
    sci[i] = 0;

  for(row = 1; row <= MAX1; ++row)   
    for(col = 1; col <= MAX2; ++col) 
      Sample(row,col) = 0;

  for(i = 0; i < num_mat; i++)
    denom += matscore[i];

  for(i = 0; i < num_mat; i++) 
    matscore[i] = (matscore[i]/denom)*1000;    

  for(int k = 0; k < num_mat; k++) {
    samMAX = (int) matscore[k];
    for(i = 1; i <= MAX1; ++i)
      p(i,0) = 0;
    for(j = 0; j <= MAX2; ++j)
      p(0,j) = 0;
    
    // Protien Matrix 
    for(i = 1; i <= MAX1; ++i) {   
      mat1 = seq1index[i];
      for(j = 1; j <= MAX2; ++j) {     
	mat2 = seq2index[j];
	p(i,j) = p_related[mat[k]][mat1][mat2];
      }  
    }   
    
    gapo = gap_open[k];
    gape = gap_exte[k];
    
    numer = 0;
    
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
    
    for(row = 1; row <= MAX1; ++row) 
      for(col = 1; col <= MAX2; ++col) 
	MatE(row,col) = MatE(row,col)/numer;
    
    srand(time(NULL));
    
    for(sam = 1; sam <= samMAX; ++sam) {  
      end_total = 0;
      randnum = rand();
      randnum = randnum/RAND_MAX;
      for(row = MAX1; row > 0; --row) {
	for(col = MAX2; col > 0; --col) {
	  end_total += MatE(row,col);
	  if(randnum < end_total) {
	    end_row = row;
	    end_col = col;
	    row = 0;
	    col = 0;
	  }
	}
      }
      
      row = end_row;      
      col = end_col;      
      denom = MatM(row,col) + MatI(row,col) + MatD(row,col) + MatN(row,col);
      MoveM = MatM(row,col)/denom;
      MoveI = MatI(row,col)/denom;
      MoveD = MatD(row,col)/denom;
      randnum = rand();
      randnum = randnum/RAND_MAX;
      if(randnum < MoveM) {
	++Sample(row,col);
	--row;
	--col;
	match = 1;
      }
      else if(randnum < (MoveM+MoveI)) {
	--row;
	insert_gap = 1;
      }
      else if(randnum < (MoveM+MoveI+MoveD)) {
	--col;     
	delete_gap = 1;
      }
      else {
	++Sample(row,col);
	++sci[row];
	row = 0;
	col = 0;
	match = 0;
      }
    
      while(row > 0 && col > 0) {
	if(match == 1) {
	  denom = MatM(row,col) + MatI(row,col) + MatD(row,col) + MatN(row,col);
	  MoveM = MatM(row,col)/denom;
	  MoveI = MatI(row,col)/denom;
	  MoveD = MatD(row,col)/denom;
	  randnum = rand();
	  randnum = randnum/RAND_MAX;
	  if(randnum < MoveM) {
	    ++Sample(row,col);
	    ++sci[row];
	    --row;
	    --col;	   
	  }
	  else if(randnum < (MoveM+MoveI)) {
	    --row;
	    match = 0;
	    insert_gap = 1;
	  }
	  else if(randnum < (MoveM+MoveI+MoveD)) {
	    --col;    
	    match = 0;
	    delete_gap = 1;
	  }
	  else {
	    ++Sample(row,col);
	    ++sci[row];
	    row = 0;
	    col = 0;
	    match = 0;
	  }
	}
	if(insert_gap == 1) {
	  denom = MatM(row,col) + MatI(row,col)+ MatN(row,col);
	  MoveM = MatM(row,col)/denom;
	  MoveI = MatI(row,col)/denom;	
	  randnum = rand();
	  randnum = randnum/RAND_MAX;
	  if(randnum < MoveM) {
	    ++Sample(row,col);
	    ++sci[row];
	    --row;
	    --col;	
	    insert_gap = 0;
	    match = 1;
	  }
	  else if(randnum < (MoveM+MoveI)) {
	    --row;	 
	  }       
	  else {
	    ++Sample(row,col);
	    ++sci[row];
	    insert_gap = 0;
	    row = 0;
	    col = 0;
	  }
	}	  
	else if(delete_gap = 1) {
	  denom = MatM(row,col) + MatD(row,col) + MatN(row,col);
	  MoveM = MatM(row,col)/denom;
	  MoveD = MatD(row,col)/denom;	  
	  randnum = rand();
	  randnum = randnum/RAND_MAX;
	  if(randnum < MoveM) {
	    ++Sample(row,col);
	    ++sci[row];
	    --row;
	    --col;	
	    delete_gap = 0;
	    match = 1;
	  }      
	  else if(randnum < (MoveM+MoveD)) {
	    --col;     	  
	  }
	  else {
	    ++Sample(row,col);
	    ++sci[row];
	    row = 0;
	    col = 0;
	    delete_gap = 0;
	  }
	}
      }
    }
  }

  ofstream histfile(histname);
  for(row = 1; row <= MAX1; ++row) {  
    for(col = 1; col <= MAX2; ++col) {
//      histfile << row << "     " << col << "     " << Sample(row,col)/1000 << endl;
      float tmpfloat = Sample(row,col)/1000;
      if (tmpfloat != 0)
        histfile << row << " " << col << " " << tmpfloat << endl;
    }
  }
  histfile.close();

  return NULL;
}
      
void Align::SetFileNames(char *name)
{
  strcpy( histname, name );
  strcat( histname, "_histogram" );
  strcpy( centname, name );
  strcat( centname, "_centroid" );
  strcpy( samplename, name );
  strcat( samplename, "_samples" );
  strcpy( credibilityname, name );
  strcat( credibilityname, "_credibility" );
}
