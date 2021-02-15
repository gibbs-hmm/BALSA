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

//BALSA_database.C
//  This is the version of BALSA which reads in one sequence from a file given by the user
//  and align it against the database

#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include "matrix.h"
#include "letter.h"
#include "BALSA_file_score_sums.h"
#include "BALSA_file_align.h"
#include "BALSA_file_get_seq.h"
#include "Blosum_name.h"

using namespace std;


#define FNAME_PREFIX	"histogram."
#define GAPO_UPPER	-6
#define GAPO_LOWER	-20
#define GAPE_UPPER	-1
#define GAPE_LOWER	-10
#define OPTION_SIZE	50

int mat[4], num_mat;
double matscore[4];
float gap_open[4];
float gap_exte[4];

char const *alignhist="histodisp -FILENAME=histogram.dat -LABELX=SEQ1 -LABELY=SEQ2 -LABELZ=PROBABILITY_MATCH";
char amino_letter[] = { "ARNDCQEGHILKMFPSTWYVX" };
char dna_letter[] = {"ATCGN"};
char *letter = amino_letter;

// local function declaration
void print_usage(char *);
int get_dbsize(char *);

int main(int argc, char* argv[])
{
//  int numberalign, num, numindex[1500];
  int numberalign, num;
  int i, j, k, MAX1, MAX2, headerlength;
  char Header[1000], SeqFile[MAX_SEQ], *seq1, *seq2;
  int seq1index[MAX_SEQ], seq2index[MAX_SEQ];
  double Score, BayesianPValue, posterior[4], gapo[4], gape[4];
  double tmpdouble;
  char fname_ext[25], fname[25];
  char tmpstr[MAX_SEQ+1];
  Score_Sums score;
  Align align;
  char option[OPTION_SIZE];
  char *tmpptr, *tmpptr2, *infname, *dbfname, *outfname;
  long db_homologs = 0, dbsize;
  double final_term = 0, pvalue = 0.01;

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
    if (strcasecmp(option, "-INFILE") == 0) {
      infname = &argv[i][j];
      if (strlen(infname) == 0) {
        cout << "Error: no input filename specified" << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
    } else if (strcasecmp(option, "-DBFILE") == 0) {
      dbfname = &argv[i][j];
      if (strlen(dbfname) == 0) {
        cout << "Error: no database filename specified" << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
    } else if (strcasecmp(option, "-DB_HOMOLOGS") == 0) {
      tmpptr = &argv[i][j];
      db_homologs = atoi(tmpptr);
      if (db_homologs <= 0) {
        cout << "Error: invalid number of homologs " << tmpptr << " in database" << endl << endl;
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
    } else if (strcasecmp(option, "-PVALUE") == 0) {
      tmpptr = &argv[i][j];
      pvalue = atof(tmpptr);
      if (pvalue <= 0 || pvalue >= 1) {
        cout << "Error: invalid level of significance for a return: " << tmpptr << endl << endl;
        print_usage(argv[0]);
        return -1;
      }
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

  if (infname == NULL) {
    cout << "Error: no input filename specified" << endl << endl;
    print_usage(argv[0]);
    return -1;
  }
  if (dbfname == NULL) {
    cout << "Error: no database filename specified" << endl << endl;
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

  // Get Sequence File
  ifstream sequenceFile(infname);
  if (!sequenceFile) {
    cout << "File <" << infname << "> could not be opened" << endl;
    exit(1);
  }
  i = 1;
  sequenceFile.get(SeqFile[i]);
  while(SeqFile[i] != 10)
    sequenceFile.get(SeqFile[i]);
  sequenceFile.get(SeqFile[i]);
  while(!sequenceFile.eof()) {
    if(SeqFile[i] != 10)
      if (++i >= MAX_SEQ) {
        cout << "Error: the length of input sequence exceeds the current limit of " << MAX_SEQ << endl;
        exit(1);
      }
    sequenceFile.get(SeqFile[i]);
  }
  MAX1 = i-1;
  seq1 = new char[i];
  for(j = 1; j <= MAX1; ++j) {
    seq1[j] = SeqFile[j];
  }
  sequenceFile.close();

//  char *tmpptr;
        
  for(i = 1; i <= MAX1; ++i) {
    if ((tmpptr=strchr(letter, seq1[i])) == NULL) {
      if (seq1[i] == 'x')
        seq1index[i] = sizeof(letter)-2;
      else {
        cout << "Error: unknown letter '" << seq1[i] << "' in input sequence" << endl;
        exit(1);
      }
    } else seq1index[i] = tmpptr - letter;
  }

  dbsize = get_dbsize(dbfname);
  if (db_homologs > 0)
    final_term = (double) db_homologs / ((dbsize*dbsize-dbsize)/2 - db_homologs);
  else
    final_term = 1.0 / dbsize;

  ifstream SCOPfile(dbfname);
  if (!SCOPfile) {
    cout << "Database file <" << dbfname << "> could not be opened" << endl;
    exit(1);
  }

  ofstream outfile(outfname);
  if (!outfile) {
    cout << "Output file could not be created" << endl;
    exit(1);
  }

  num = 0;
  for(i = 1; i <= dbsize; ++i) {

    while (!SCOPfile.eof()) {
      SCOPfile.getline(Header, MAX_SEQ);
      if (Header[0] == '>')
        break;
    }
    j = strlen(Header);
    if (j > 0) {
      if (Header[j-1] == '\r')
        j--;
    }
    headerlength = j;

    strcpy(SeqFile, "");
    char firstchar;
    while (!SCOPfile.eof()) {
      firstchar = SCOPfile.get();
      SCOPfile.unget();
      if (firstchar == '>')
        break;
      SCOPfile.getline(tmpstr, MAX_SEQ);

      // we should not count \r at the end of a string
      j = strlen(tmpstr);
      if (j > 0) {
        if (tmpstr[j-1] == '\r')
          j--;
      }

      if (j == 0)
        break;
      else
        strncat(SeqFile, tmpstr, j);
    }

    MAX2 = strlen(SeqFile);
    seq2 = new char[MAX2+1];


    for(j = 1; j <= MAX2; ++j) {
      seq2[j] = SeqFile[(j-1)]; 
    }

    for(j = 1; j <= MAX2; ++j) {
      if ((tmpptr=strchr(letter, seq2[j])) == NULL) {
        if (seq2[j] == 'x')
          seq2index[j] = sizeof(letter)-2;
        else {
          cout << "Error: unknown letter '" << seq2[j] << "' in database sequence " << i << endl;
          exit(1);
        }
      } else seq2index[j] = tmpptr - letter;
    }

    Score = 0;
    for(j = 0; j < num_mat; j++) {
      score.ScoreMatrix(MAX1, MAX2, seq1index, seq2index, j);
      matscore[j] = score.PR1R2_theta;
      Score += matscore[j]/num_mat;
    }

    BayesianPValue = (Score*final_term) + 1; 
    BayesianPValue = 1/BayesianPValue;
    if(BayesianPValue < pvalue) {
      ++num;

      // create output file name from num
      strcpy(fname, FNAME_PREFIX);
      sprintf(fname_ext, "%d", num);
      strcat(fname, fname_ext);
      outfile << fname << endl;

//      numindex[num] = i;
      for(j = 0; j < headerlength; ++j) 
	outfile << Header[j];
      outfile << endl;
      for(j = 1; j <= MAX2; ++j) {
	outfile << seq2[j];
        if (j % 71 == 0 && j != MAX2)
          outfile << endl;
      }
      outfile << endl;
      for (j = 0; j < num_mat; j++) {
        posterior[j] = matscore[j] / num_mat / Score;
        outfile << "P(" << name_related[mat[j]] << ", Gap Opening Penalty=";
        outfile << gapo[j] << ", Gap Extension Penalty=" << gape[j];
        outfile << " | R1,R2) = " << posterior[j] << endl;
      }
      outfile << "Bayesian P-Value = " << BayesianPValue << endl << endl;

      align.SetFileNames(outfname);

      align.Sample_Alignment(MAX1, MAX2, seq1index, seq2index, matscore);
//      system(alignhist);
    }
    delete[] seq2;
  }

  delete[] seq1;

  outfile.close();
  SCOPfile.close();

  if (num == 0) {
    cout << "There are no sequences with the Bayesian P-Value less than " << pvalue << endl;
  }

  return 0;
}


void print_usage(char *prg_name) {
  cout << "Usage: " << prg_name << " -INFILE=infile -DBFILE=dbfile -OUTFILE=outfile" << endl;
  cout << "        [-MATRIX=scoring_matrix,gapo,gape]" << endl;
  cout << "        [-DB_HOMOLOGS=number_of_homologs_in_database]" << endl;
  cout << "        [-PVALUE=level_of_significance_for_a_return (def=0.01)]" << endl << endl;
  cout << " infile is the file containing input sequence in FASTA format;" << endl;
  cout << " dbfile is the database file to be used;" << endl;
  cout << " outfile is the output file;" << endl;
  cout << " scoring_matrix can be BLOSUM_30, BLOSUM_40, BLOSUM_45," << endl;
  cout << "                       BLOSUM_50, BLOSUM_55, BLOSUM_62," << endl;
  cout << "                       BLOSUM_70, BLOSUM_80 or PAM250;" << endl;
  cout << " gapo is the gap opening penalty ("<< GAPO_LOWER << " < gapo < " << GAPO_UPPER << ");" << endl;
  cout << " gape is the gap extension penalty ("<< GAPE_LOWER << " < gape < " << GAPE_UPPER << ");" << endl;
  cout << "Note: At least 1 and at most 4 <scoring_matrix,gapo,gape> values can be specified" << endl;
}


// get the number of sequences in the database file
int get_dbsize(char *fname) {
  ifstream DBfile(fname);
  char tmpstr[MAX_SEQ+1];
  int dbsize = 0, len = 0, tmplen;

  if (!DBfile) {
    cout << "Database file <" << fname << "> could not be opened" << endl;
    exit(1);
  }

  while (!DBfile.eof()) {
    DBfile.getline(tmpstr, MAX_SEQ);

    // this check is to prevent a single line to be longer than MAX_SEQ
    if (strlen(tmpstr) == MAX_SEQ-1) {
      cout << "Error: the length of a database sequence exceeds the current limit of " << MAX_SEQ << endl;
      exit(1);
    }

    // we should not count \r at the end of a string
    tmplen = strlen(tmpstr);
    if (tmpstr[tmplen-1] == '\r')
      tmplen--;

    if (tmpstr[0] == '>') {
      if (len >= MAX_SEQ) {
        cout << "Error: the length of a database sequence exceeds the current limit of " << MAX_SEQ << endl;
        exit(1);
      }
      len = 0;
      dbsize++;
    } else len += tmplen;
  }
  if (len >= MAX_SEQ) {
    cout << "Error: the length of a database sequence exceeds the current limit of " << MAX_SEQ << endl;
    exit(1);
  }

  DBfile.close();

  return dbsize;
}
