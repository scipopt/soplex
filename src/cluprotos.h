/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2001 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cluprotos.h,v 1.5 2001/11/29 14:00:25 bzfkocht Exp $"


#ifndef _CLUPROTOS_H_
#define _CLUPROTOS_H_


#define WITH_L_ROWS             1

#define CLU_OK                  0
#define CLU_INSTABLE            1
#define CLU_SINGULAR            2
#define CLU_OUT_OF_MEMORY       8

#include "svector.h"
#include "clutypes.h"

namespace soplex
{

int factor(CLUFactor*,
            SVector**,       /* Array of column vector pointers   */
            double,          /* pivoting threshold                */
            double           /* epsilon for zero detection        */
         );

void solveRight (CLUFactor*, double*, double*);
int solveRight4update (CLUFactor*, double*, int*, double,
                        double*, double*, int*, int*);

void solveRight2(CLUFactor* fac,
                  double* vec1,
                  double* vec2,
                  double* rhs1,
                  double* rhs2);

int solveRight2update(CLUFactor* fac,
                       double* vec1,
                       double* vec2,
                       double* rhs1,
                       double* rhs2,
                       int* nonz,
                       double eps,
                       double* tmp,
                       int* forestNum,
                       int* forestIdx);

void solveLeft (double*, CLUFactor*, double*);
int solveLeftEps (double*, CLUFactor*, double*, int*, double);

int solveLeft2(CLUFactor* fac,
                double* vec1,
                int* nonz,
                double* vec2,
                double eps,
                double* rhs1,
                double* rhs2);

int updateCLUFactor (CLUFactor*, int, double*, const int*, int);
int updateCLUFactorNoClear(CLUFactor*, int, const double*, const int*, int);
int forestUpdateCLUFactor (CLUFactor* fac, int col, double* work, int n, int *nonz);
int CLUFactorIsConsistent (const CLUFactor*);

void solveLright(CLUFactor* fac, double* vec);

/********************************************************************************
        very sparse solution methods
*/
int vSolveRight4update(CLUFactor* fac, double eps,
                        double* vec, int* idx,               /* result       */
                        double* rhs, int* ridx, int rn,      /* rhs & Forest */
                        double* forest, int* forestNum, int* forestIdx);
int vSolveRight4update2(CLUFactor* fac, double eps,
                         double* vec, int* idx,                      /* result1 */
                         double* rhs, int* ridx, int rn,             /* rhs1    */
                         double* vec2, double eps2,                  /* result2 */
                         double* rhs2, int* ridx2, int rn2,          /* rhs2    */
                         double* forest, int* forestNum, int* forestIdx);
void vSolveRightNoNZ(CLUFactor* fac,
                      double* vec2, double eps2,              /* result2 */
                      double* rhs2, int* ridx2, int rn2);   /* rhs2    */

int vSolveLeft(CLUFactor* fac, double eps,
                double* vec, int* idx,                               /* result */
                double* rhs, int* ridx, int rn);                   /* rhs    */

void vSolveLeftNoNZ(CLUFactor* fac, double eps,
                     double* vec,                             /* result */
                     double* rhs, int* ridx, int rn);       /* rhs    */

int vSolveLeft2(CLUFactor* fac, double eps,
                 double* vec, int* idx,                              /* result */
                 double* rhs, int* ridx, int rn,                     /* rhs    */
                 double* vec2,                                       /* result2 */
                 double* rhs2, int* ridx2, int rn2);               /* rhs2    */

/********************************************************************************/
void remaxRow(CLUFactor*, int, int);
void remaxCol(CLUFactor*, int, int);
int makeLvec(CLUFactor*, int, int);
void packRows(CLUFactor* fac);

void dumpCLUFactor(const CLUFactor* fac);

} // namespace soplex
#endif // _CLUPROTOS_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
