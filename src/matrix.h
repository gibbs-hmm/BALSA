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

// Declaration of the MATRIX class.

#ifndef _MATRIX_H_
#define _MATRIX_H_

class Matrix
{
public:
  Matrix(int MaxRows = 1, int MaxCols = 1)
    {
      maxCols = MaxCols;
      maxRows = MaxRows;
      pData = new long double[maxRows * maxCols];
    }
  ~Matrix()
    { delete [] pData; }
  void ReSize(int NewRows, int NewCols)
    {
      delete [] pData;
      maxCols = NewCols;
      maxRows = NewRows;
      pData = new long double[maxRows * maxCols];
    }
  long double& operator ()(int row, int col)
    { return pData[row + maxRows * col]; }
  int checkRowCol(int row, int col)
    { return (row >= 0 && row < maxRows && col >= 0 && col < maxCols) ? 1: 0; }
  int getRows()
    { return maxRows; }
  int getCols()
    { return maxCols; }
  
protected:
  int maxRows;
  int maxCols;
  long double *pData;
};


#endif
