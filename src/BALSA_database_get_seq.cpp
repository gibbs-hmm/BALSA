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

// get_seq.C
// Member function definitions for Get_Seq class
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "letter.h"
#include "BALSA_database_get_seq.h"

// Time constructor initializes each data member to zero.
// Ensures all Time objects start in a consistent state.
Get_Seq::Get_Seq() { 
}


void Get_Seq::readSequence1()
{
  char fileName[80];
  cout << "File name: ";
  cin >> fileName;
  
  ifstream sequenceFile(fileName);
  if (!sequenceFile) {
    cout << "File could not be opened" << endl;
    exit(1);
  }
  i = 1;
  sequenceFile.get(SeqFile[i]);
  while(SeqFile[i] != 10)
    sequenceFile.get(SeqFile[i]);
  sequenceFile.get(SeqFile[i]);
  while(SeqFile[i] != 62) {      
    if(SeqFile[i] != 10)	
      ++i;
    sequenceFile.get(SeqFile[i]);    
  }    
  MAX1 = i-1;
  seq1 = new char[i];    
  for(j = 1; j <= MAX1; ++j)
    seq1[j] = SeqFile[j];
  
  for(i = 1; i <= MAX1; ++i) {
    for(l = 0; l< 22; ++l) {
      if(seq1[i] == letter[l])
	seq1index[i] = l;
    }
  }
}


void Get_Seq::readSequence2(int num2)
{
  pos = 0;
  num2_temp = 0;

  ifstream sequenceFile("/home/bjwebb/SCOP/SCOP40");
    
  if (!sequenceFile) {
    cout << "File could not be opened" << endl;
    exit(1);
  }

  while(num2_temp <= num2) {
    sequenceFile.get(s[pos]);
    if(s[pos] == 62) 
      ++num2_temp;
  }
 
  while(s[pos]!= 10) 
    sequenceFile.get(s[pos]);

  sequenceFile.get(s[pos]);
  while(s[pos] != 62) {
    ++pos;
    sequenceFile.get(s[pos]);
    if(s[pos] == 10)      
      --pos;
  }   

  MAX2 = pos;
  seq2 = new char[MAX2];
  for(j = 0; j < MAX2; ++j) {
    seq2[j] = s[j];
    //   cout << seq2[j];
  }
  //  cout << endl;
  sequenceFile.close();

  for(j = 1; j <= MAX2; ++j) {     
    for(l = 0; l < 22; ++l) {
      if(seq2[j-1] == letter[l])
	seq2index[j] = l;
    }	   
  }
}



