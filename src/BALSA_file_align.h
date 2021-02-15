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

// Declaration of the ALIGN class.
// Member functions are defined in the align.C

#ifndef ALIGN_H
#define ALIGN_H


using namespace std;
// static float gap_open[20] = {0.015,0.008,0.015,0.016};
// static float gap_exte[20] = {0.7,0.5,0.7,0.5};
extern float gap_open[];
extern float gap_exte[];
extern int num_mat;
extern int mat[];

class Align {
public:
  Align();
  int * Sample_Alignment(int, int, int *, int *, double *); // make and return samples matrix which store indice to include 1
  int * Centroid_Alignment(int MAX1,int MAX2);// make and return centroid matrix
  void Credibility_Limit( int MAX1,int MAX2,int *Centroid_Indx , int *Samples_Indx ); // make and print hamming distance
  static int compare (const void * a, const void * b); // used in quick sort

  void SetFileNames(char *name);
  void SetSampling() {samples = 1;}
  

private: 
  void DumpSamples( int *Samples );

  static const int Sample_Size = 1000; // initialization of sample size	 
  int i;
  int j;
  int matrixnum;
  int mat1;
  int mat2;
  float gapo;
  float gape;
  int sam;
  int row;
  int col;
  long double numer;
  long double end_total;
  int end_row;
  int end_col;
  float randnum;
  long double MoveM;
  long double MoveI;
  long double MoveD;
  int samMAX;
  long double denom;
  int match;
  int insert_gap;
  int delete_gap;
  int *sci;
//  int mat;
  int intscore[5];

  char histname[FILENAME_MAX];
  char centname[FILENAME_MAX];
  char samplename[FILENAME_MAX];
  char credibilityname[FILENAME_MAX];
  int  samples;
};

#endif
