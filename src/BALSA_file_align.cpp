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
#include "msdefs.h"
#include "bay_matrix_all.h"
#include "BALSA_file_align.h"
#include "BALSA_file_score_sums.h"


using namespace std;


//Time constructor
Align::Align() {
  samples = 0;
}


// Finds mode of gap opening and gap extension parameters
int * Align::Sample_Alignment(int MAX1, int MAX2, int *seq1index, int *seq2index, double *matscore)
{ 
  Matrix p(MAX1+1,MAX2+1);
  Matrix MatM(MAX1+1, MAX2+1);
  Matrix MatI(MAX1+1, MAX2+1);
  Matrix MatD(MAX1+1, MAX2+1);
  Matrix MatN(MAX1+1, MAX2+1);
  Matrix MatE(MAX1+1, MAX2+1);  
  Matrix Sample(MAX1+1, MAX2+1);
  Matrix SampleTotal(MAX1+1, MAX2+1);
  Matrix Tmp_Sample(MAX1+1,MAX2+1);
  int *Samples_Matrix = new int[Sample_Size * (MAX1 + MAX2 + 2)]; // declare the toal sample matrix with one dimentaional array
  
  for(int i = 0; i<Sample_Size * (MAX1 + MAX2) ; ++i) // initialization of Samples_Matrix
    *(Samples_Matrix + i) = 0;

  sci = new int[MAX1+1];

  //for(row = 1; row <= MAX1; ++row)   

  //for(col = 1; col <= MAX2; ++col) 

  //  Sample(row,col) = 0;

  for(int i = 1; i <= MAX1; ++i)
    sci[i] = 0;

  for(int i = 0; i < num_mat; i++)
    denom += matscore[i];
   
  intscore[num_mat-1] = Sample_Size;
  for(int i = 0; i < num_mat-1; i++) {
    matscore[i] = (matscore[i]/denom)*100;    
    intscore[i] = (int) matscore[i];
    if((matscore[i] - intscore[i]) > 0.5)
      ++intscore[i];
    intscore[i] = intscore[i]*10;   
    intscore[num_mat-1] -= intscore[i];
  }

  /*for(mat = 0; mat <= 3; ++mat) {
    cout << intscore[mat] << endl;
    }*/
  for(int k=0; k<num_mat; k++) {

    for(row = 1; row <= MAX1; ++row)   
      for(col = 1; col <= MAX2; ++col) 
	Sample(row,col) = 0;

    samMAX = intscore[k];
    for(int i = 1; i <= MAX1; ++i)
      p(i,0) = 0;
    for(int j = 0; j <= MAX2; ++j)
      p(0,j) = 0;
    
    // Protien Matrix 
    for(int i = 1; i <= MAX1; ++i) {   
      mat1 = seq1index[i];
      for(int j = 1; j <= MAX2; ++j) {     
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
	//	cout << row << " " << col << " " << MatE(row,col) << " " << numer << endl;
      }
    }
    
    for(row = 1; row <= MAX1; ++row) 
      for(col = 1; col <= MAX2; ++col) 
	MatE(row,col) = MatE(row,col)/numer;
    
    srand(time(NULL));    

    int Indx = 0; // This is index of Samples
    
    for(row=1;row<=MAX1;++row)  // initialization for storing each accumulated samples      
      for(col=1;col<=MAX2;++col)	
	Tmp_Sample(row,col)=0;
    
    for(sam = 1; sam <= samMAX; ++sam) { // start sampling align.Sample_Size times         	 	
      end_total = 0;
      randnum = rand();
      randnum = randnum/RAND_MAX;
      int findlocation=0;
      for(row = MAX1; row > 0; --row) {
	for(col = MAX2; col > 0; --col) {
	  end_total += MatE(row,col);
	  if(randnum < end_total) {
	    end_row = row;
	    end_col = col;
	    row = 0;
	    col = 0;
	    findlocation = 1;
	    break;
	  }
	}	   
	if(findlocation == 1)
	  break;
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
	else if(delete_gap == 1) {
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
	  
      for(row = 1; row <= MAX1; ++row){  // store each sample's index which is including 1 into samples matrix
	for(col = 1; col <= MAX2; ++col){
	  if (Sample(row,col) - Tmp_Sample(row,col) != 0){
	    *(Samples_Matrix + Indx) = row;
	    *(Samples_Matrix + Indx + 1) = col;
	    Indx = Indx + 2;
	  }
	  Tmp_Sample(row,col) = Sample(row,col);  // store previous accumulated sample into temporary matrix
	}
      }
  
      *(Samples_Matrix + Indx) = -1; // indicator of each sample

      Indx++;  // index which start next sample to store its index which is including 1 into samples matrix
	
      // ********************************* samples matrix completion ******************************************
   

    }    // sample end bracket

    for(row = 1; row <= MAX1; ++row){
      for(col = 1; col <= MAX2; ++col){
	SampleTotal(row,col) += Sample(row,col);
      }
    }

    // cout<<*(Samples_Matrix + Indx-1)<<"   HHHHHHHHHH"<<endl;
  }     // num_mat end bracket 

 
  ofstream histfile(histname);
  for(row = 1; row <= MAX1; ++row) {  
    for(col = 1; col <= MAX2; ++col) {
      //      histfile << row << "     " << col << "     " << Sample(row,col)/align.Sample_Size << endl;
      double tmpfloat = SampleTotal(row,col)/Sample_Size;
      if (tmpfloat != 0)
        histfile << row << " " << col << " " << tmpfloat << endl;
    }
  }
  histfile.close();

  return Samples_Matrix;
}

int * Align::Centroid_Alignment(int MAX1,int MAX2)
{
	ifstream samplefile(histname);
	ofstream centroidfile(centname);
	Matrix Sample(MAX1+1,MAX2+1);
	int *Centroid = new int[MAX1 + MAX2]; // declare centroid matrix to store indice which is including 1. this can not exceed the sum of two sequence size
	
	for(int i =0; i < MAX1 + MAX2; ++i)
		*(Centroid + i) = 0;            // initialization of centroid matrix

	for(int row = 1; row <= MAX1; ++row)
		for(int col = 1; col <= MAX2; ++col)
			Sample(row,col)=0;		
	
	while(! samplefile.eof())  // store probabilities into Sample matrix from text file
	  {
	    int a = 0;
	    int b = 0;
	    double c;
	    samplefile>>a>>b>>c;
	    Sample(a,b)=c;
	  }

	int Indx=0;
	for(int row = 1;row <= MAX1; ++row){   // store indice which are greater than 0.5 into Centroid matrix 
		for(int col = 1; col <= MAX2; ++col){
			if(Sample(row,col)>0.5){
				*(Centroid + Indx) = row ;
				*(Centroid + (Indx+1)) = col;
				Indx = Indx + 2;
			}
		}
	}
	
    for(int row = 1; row <= MAX1; ++row) {    // write  row , column indice and probabilities into centroid.dat 
        for(int col = 1; col <= MAX2; ++col) {
           float tmpfloat = Sample(row,col);
           if (tmpfloat > 0.5)
           centroidfile << row << " " << col << " " << tmpfloat << endl;
        }
    }
    centroidfile.close();
	samplefile.close();
	return Centroid;
}

void Align::Credibility_Limit( int MAX1,int MAX2,int *Centroid_Indx , int *Samples_Indx )
{
  int * Centroid_Matrix = Centroid_Indx;   
  int * Samples_Matrix = Samples_Indx;	   
  int Distance_Array[Sample_Size]; // declare array of hamming distance 

  for (int i = 0 ; i < Sample_Size ; ++i)	// initialization of array of hamming distance
    Distance_Array[i] = 0;
	
  int score1 = 0; // score of centroid alignment     
  int Indx=0;
  double Element =*Centroid_Matrix;

  while ( Element != 0 ){ // evaluate the frequency of 1'value from the centroid matrix
    score1++;
    Indx=Indx+2;
    Element=*(Centroid_Matrix + Indx);
  }
 
  int Distance_Indx = 0;
  int Samples_Matrix_Indx = 0;
  int max_sample_length = -1;; 
	
  while (*(Samples_Matrix + Samples_Matrix_Indx) != 0 ) // evaluate the hamming distance from centroid alignment 
    {                                                     // to each sample alignment
      int score2 = 0;
      int common_count = 0;
      int i = 0;
      int j;
	    
      for ( i = Samples_Matrix_Indx ; *(Samples_Matrix + i ) != -1 ; i = i + 2) 
	{
	  int row1 = 0; 
	  int col1 = 0; 
	  row1 = *(Samples_Matrix + i);
	  col1 = *(Samples_Matrix + i + 1);
		
	  for ( j = 0; *(Centroid_Matrix + j) != 0 ; j = j + 2)
	    {
	      int row2 = 0;
	      int col2 = 0;
	      row2 = *(Centroid_Matrix + j);
	      col2 = *(Centroid_Matrix + j + 1);

	      if ( row1 == row2 && col1 == col2)
		common_count++;   // count for common pair from the two sequcne alignments
	    }
	  score2++;
	}
	    
      Samples_Matrix_Indx = i + 1; // starting point for next sample matrix
      Distance_Array[Distance_Indx++] = ( score1 - common_count) + (score2 - common_count);
      max_sample_length = max(max_sample_length, score2); // BT 06/30/08

      // cout<<score1<<","<<score2<<","<<common_count<<","<<( score1 + score2 ) - common_count << " " << max_sample_length <<endl;
    }   

 
  qsort (Distance_Array,Sample_Size, sizeof(int), &Align::compare); // increase sorting for array of hamming distance using  	
  
  ofstream credibilityfile(credibilityname);
  
  /* for(int i=0; i<20;++i)
     for(int j=0;j<50;++j)
     cout<<Distance_Array[i*50+j]<<" , ";
     cout<<endl;*/
  

//   credibilityfile << "               ************ Rank of Hamming Distance ************\n\n"; 
//   j = 0;
//   for (int i = 0; i < 20; i++)	{     // quick sort which is already implemented from the library      
	
//     for (int j = 0; j < 50; j++)	{
//       credibilityfile<< Distance_Array[i * 50 + j] <<","; //print the array of hamming distance
//     }
//     credibilityfile<< endl;
//   }

//   credibilityfile<<"\n\n\n";
  credibilityfile << "               ************ Credibility Limit(85,90,95) ************\n\n";
  int percent [] = {85,90,95};
  for (int i=0;i<3;++i)
    {
      int loc = int(Sample_Size *  (double(percent[i]) / 100));  // BT 06/30/08
      double norm_credibility = double(Distance_Array[loc]) / (score1 + max_sample_length);
      credibilityfile<< "The credibility limit of  " << percent[i] << " %  :  " << Distance_Array[loc] << " " << norm_credibility << endl;
    }
  
  credibilityfile.close();  

  if(samples)
    DumpSamples( Samples_Indx );
}


void Align::DumpSamples( int *Samples )
{
  ofstream samplefile(samplename);
  
  int index = 0;
  int index2 = 0;
  while( index < Sample_Size )
    {
      while( *(Samples + index2) != -1 )
	{
	  samplefile << index << " " << *(Samples + index2) << " " << *(Samples + index2 + 1) << endl;
	  index2 += 2;
	}
      index2++;
      index++;
    }
  samplefile.close();
}



int Align:: compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
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



