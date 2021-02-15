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
#include <string.h>
#include "letter.h"
#include "msdefs.h"
#include "BALSA_file_get_seq.h"

using namespace std;

// Time constructor initializes each data member to zero.
// Ensures all Time objects start in a consistent state.
Get_Seq::Get_Seq() { 
}


void Get_Seq::readSequences(char *file1, char *file2)
{
  ifstream sequenceFile1(file1);
  if (!sequenceFile1) {
    cout << "File <" << file1 << "> could not be opened" << endl;
    exit(1);
  }
  i = 1;
  sequenceFile1.get(SeqFile[i]);
  while(SeqFile[i] != 10)
    sequenceFile1.get(SeqFile[i]);
  sequenceFile1.get(SeqFile[i]);
  while(!sequenceFile1.eof()) { 
    if(SeqFile[i] != 10)
      if (++i >= MAX_SEQ) {
        cout << "Error: the length of sequence 1 exceeds the current limit of " << MAX_SEQ << endl;
        exit(1);
      }
    sequenceFile1.get(SeqFile[i]);    
  }    
  MAX1 = i-1;
  seq1 = new char[i];    
  for(j = 1; j <= MAX1; ++j)
    seq1[j] = SeqFile[j];

  sequenceFile1.close();


  ifstream sequenceFile2(file2);
  if (!sequenceFile2) {
    cout << "File <" << file2 << "> could not be opened" << endl;
    exit(1);
  }
  i = 1;
  sequenceFile2.get(SeqFile[i]);
  while(SeqFile[i] != 10)
    sequenceFile2.get(SeqFile[i]);
  sequenceFile2.get(SeqFile[i]);
  while(!sequenceFile2.eof()) { 
    if(SeqFile[i] != 10)
      if (++i >= MAX_SEQ) {
        cout << "Error: the length of sequence 2 exceeds the current limit of " << MAX_SEQ << endl;
        exit(1);
      }
    sequenceFile2.get(SeqFile[i]);          
  }
  MAX2 = i-1;
  seq2 = new char[i];    
  for(j = 1; j <= MAX2; ++j)
    seq2[j] = SeqFile[j];

  sequenceFile2.close();

  char *tmpptr;

  for(i = 1; i <= MAX1; ++i) {
    if ((tmpptr=strchr(letter, seq1[i])) == NULL) {
      if (seq1[i] == 'x')
        seq1index[i] = sizeof(letter)-2;
      else {
        cout << "Error: unknown letter '" << seq1[i] << "' in sequence 1" << endl;
        exit(1);
      }
    } else seq1index[i] = tmpptr - letter;
  }

  for(j = 1; j <= MAX2; ++j) {
    if ((tmpptr=strchr(letter, seq2[j])) == NULL) {
      if (seq2[j] == 'x')
        seq2index[j] = sizeof(letter)-2;
      else {
        cout << "Error: unknown letter '" << seq2[j] << "' in sequence 2" << endl;
        exit(1);
      }
    } else seq2index[j] = tmpptr - letter;
  }
}
