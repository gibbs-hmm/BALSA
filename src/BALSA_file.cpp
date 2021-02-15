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

//BALSA_FILE.C
//  This is the version of BALSA which reads in two sequences from a file given by the user.

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cstdio>
#include "matrix.h"
#include "msdefs.h"
#include "BALSA_file_get_seq.h"
#include "BALSA_file_score_sums.h"
#include "BALSA_file_align.h"
#include "Blosum_name.h"

using namespace std;

#define GAPO_UPPER      -6
#define GAPO_LOWER      -50
#define GAPE_UPPER      -1
#define GAPE_LOWER      -10
#define OPTION_SIZE     50
#define VERSION         2.01

int mat[4];
double matscore[5], denom;
int intscore[5];
int num_mat;
float gap_open[4];
float gap_exte[4];

char amino_letter[] = { "ARNDCQEGHILKMFPSTWYVX" };

//  A   R   N   D  C
//  Q   E   G   H  I
//  L   K   M   F  P
//  S   T   W   Y  V

// char dna_letter[] = {"ABCDGHKMNRSTUVWXY"};
char dna_letter[] = {"ATCGN"};
char *letter = amino_letter;

char const *alignhist="./histodisp -FILENAME=histogram.dat -LABELX=SEQ1 -LABELY=SEQ2 -LABELZ=PROBABILITY_MATCH";

// local function declaration
void print_usage(char *);

int main(int argc, char* argv[])
{
  
  int i, j, k;
  Get_Seq seq;
  Score_Sums score;
  Align align;
  double tmpdouble, Score, posterior[4], gapo[4], gape[4];
  char option[OPTION_SIZE];
  char *tmpptr, *tmpptr2, *in1fname, *in2fname, *outfname;
  int  mat_type[4] = {0,0,0,0};

  num_mat = 0;
  for (i=1; i<argc; i++) {
    j = 0;
    while (argv[i][j] != '=' && argv[i][j] != '\0' && j < (OPTION_SIZE-1)) {
      option[j] = argv[i][j];
      j++;
    }
    option[j] = '\0';
    if (argv[i][j] != '\0')  
      j++;
    if (strcasecmp(option, "-INFILE1") == 0) {
      in1fname = &argv[i][j];
      if (strlen(in1fname) == 0) {
        cout << "Error: no input file 1 specified" << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
    } else if (strcasecmp(option, "-INFILE2") == 0) {
      in2fname = &argv[i][j];
      if (strlen(in2fname) == 0) {
        cout << "Error: no input file 2 specified" << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
    } else if (strcasecmp(option, "-OUTFILE") == 0) {
      outfname = &argv[i][j];
      if (strlen(outfname) == 0) {
        cout << "Error: no output filename specified" << endl << endl;
        print_usage(argv[0]);
        return -1;
      }      
    } else if (strcasecmp(option, "-SAMPLES") == 0) {
      align.SetSampling();
    } else if (strcasecmp(option, "-MATRIX") == 0) {
      if (num_mat >= 4) {
        // at most 4 scoring matrices can be specified
        cout << "Error: at most 4 scoring matrices can be specified" << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
      tmpptr = &argv[i][j];
      tmpptr2 = strtok(tmpptr, ",");
       
      if (tmpptr2 != NULL) {
        bool match = false;
        for (k=0; k<sizeof(name_related)/sizeof(*name_related); k++)
          if (strcasecmp(tmpptr2, name_related[k]) == 0) {
            match = true;
            mat[num_mat] = k;
	    if( strstr( tmpptr2, "DNA" ) )
	      {
		letter = dna_letter;
		mat_type[num_mat] = 1;
	      }
            break;
          }
        if (!match) {
          cout << "Error: unknown scoring matrix " << tmpptr2 << endl << endl;
          print_usage(argv[0]);
          return -1;
        }
      } else {
        cout << "Error: invalid matrix option " << tmpptr << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
      
      tmpptr2 = strtok(NULL, ",");
      if (tmpptr2 != NULL) {
        tmpdouble = atof(tmpptr2);
        if (tmpdouble >= 0 || tmpdouble > GAPO_UPPER || tmpdouble < GAPO_LOWER) {
          cout << "Error: invalid gap opening value " << tmpptr2 << endl << endl;
          print_usage(argv[0]);
          return -1;
        }
        gapo[num_mat] = tmpdouble;
        gap_open[num_mat] = pow(2.0, tmpdouble/2);
      } else {
        cout << "Error: invalid matrix option " << tmpptr << endl << endl;
        print_usage(argv[0]);
        return -1;
      }

      tmpptr2 = strtok(NULL, "");
      if (tmpptr2 != NULL) {
        tmpdouble = atof(tmpptr2);
        if (tmpdouble >= 0 || tmpdouble > GAPE_UPPER || tmpdouble < GAPE_LOWER) {
          cout << "Error: invalid gap extension value " << tmpptr2 << endl << endl;
          print_usage(argv[0]);
          return -1;
        }
        gape[num_mat] = tmpdouble;
        gap_exte[num_mat] = pow(2.0, tmpdouble/2);
      } else {
        cout << "Error: invalid matrix option " << tmpptr << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
      num_mat++;
    } else {
      cout << "Illegal option " << option << endl << endl;
      print_usage(argv[0]);
      return -1;
    }
  }
          
  for( i = 1; i < num_mat; i++ )
    {
      if( mat_type[i-1] != mat_type[i] )
	{
	  cout << "Error: Amino acid and dna matrix types can't be mixed." << endl << endl;
	  print_usage(argv[0]);
	  return -1;
	}
    }

  if (in1fname == NULL) {
    cout << "Error: no input file 1 specified" << endl << endl;
    print_usage(argv[0]);
    return -1;
  }
  if (in2fname == NULL) {
    cout << "Error: no input file 2 specified" << endl << endl;
    print_usage(argv[0]);
    return -1;
  }
  if (outfname == NULL) {
    cout << "Error: no output filename specified" << endl << endl;
    print_usage(argv[0]);
    return -1;
  }
  if (num_mat == 0) {
    cout << "Error: at least 1 scoring matrix has to be specified" << endl << endl;
    print_usage(argv[0]);
    return -1;
  }


  seq.readSequences(in1fname, in2fname);
  
  denom = 0;
  Score = 0;
  for(i=0; i<num_mat; i++) {
    score.ScoreMatrix(seq.MAX1, seq.MAX2, seq.seq1index, seq.seq2index, i);
    matscore[i] = score.PR1R2_theta;
    Score += matscore[i]/num_mat;
  }

  ofstream outfile(outfname);
  if (!outfile) {
    cout << "Output file <" << outfname << "> could not be created" << endl;
    exit(1);
  }

  for (j = 0; j < num_mat; j++) {
    posterior[j] = matscore[j] / num_mat / Score;
    outfile << "P(" << name_related[mat[j]] << ", Gap Opening Penalty=";
    outfile << gapo[j] << ", Gap Extension Penalty=" << gape[j];
    outfile << " | R1,R2) = " << posterior[j] << endl;
  }

  outfile.close();

  align.SetFileNames(outfname);

  int *Samples_Indx = align.Sample_Alignment(seq.MAX1, seq.MAX2, seq.seq1index, seq.seq2index, matscore); 
 
//  Matrix Centroid_Matrix;
  int *Centroid_Indx = align.Centroid_Alignment(seq.MAX1, seq.MAX2);

// Hamming Distance

  align.Credibility_Limit(seq.MAX1, seq.MAX2, Centroid_Indx , Samples_Indx );
  //  system(alignhist);
delete [] Samples_Indx;  // memory deletion of sample matrix
delete [] Centroid_Indx; // memory deletion of centroid matrix

   return 0;

}


void print_usage(char *prg_name) {
  cout << "BALSAfile " << VERSION << " " << __DATE__ << endl << endl;
  cout << "Usage: " << prg_name << " -SAMPLES -INFILE1=in1 -INFILE2=in2 -OUTFILE=outfile" << endl;
  cout << "        [-MATRIX=scoring_matrix,gapo,gape]" << endl << endl;
  cout << " in1 is the file containing sequence 1 in FASTA format;" << endl;
  cout << " in2 is the file containing sequence 2 in FASTA format;" << endl;
  cout << " outfile is the output file;" << endl;
  cout << " scoring_matrix can be BLOSUM_30, BLOSUM_40, BLOSUM_45," << endl;
  cout << "                       BLOSUM_50, BLOSUM_55, BLOSUM_62," << endl;
  cout << "                       BLOSUM_70, BLOSUM_80, or PAM1-500_DNA or GAP_DNA;" << endl;
  cout << " gapo is the gap opening penalty ("<< GAPO_LOWER << " < gapo < " << GAPO_UPPER << ");" << endl;
  cout << " gape is the gap extension penalty ("<< GAPE_LOWER << " < gape < " << GAPE_UPPER << ");" << endl;
  cout << "Note: At least 1 and at most 4 <scoring_matrix,gapo,gape> values can be specified" << endl;
  cout << " -SAMPLES is optional; if included, a file of 1000 sampled alignment will be created." << endl;
}











