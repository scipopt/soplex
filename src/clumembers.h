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
#pragma ident "@(#) $Id: clumembers.h,v 1.5 2001/11/29 22:52:54 bzfkocht Exp $"

#ifndef _CLUMEMBERS_H_
#define _CLUMEMBERS_H_

#include "clutypes.h"

namespace soplex
{

class CLUFactor
{
public:
   int thedim;                 /* dimension of factorized matrix   */
   int stat;                   /* Status indicator. */
   int nzCnt;                  /* number of nonzeros in U      */
   double initMaxabs;     /* maximum abs number in initail Matrix */
   double maxabs;         /* maximum abs number in L and U        */

   double rowMemMult;     /* factor of minimum Memory * #of nonzeros */
   double colMemMult;     /* factor of minimum Memory * #of nonzeros */
   double lMemMult;       /* factor of minimum Memory * #of nonzeros */

   Perm row;            /* row permutation matrices */
   Perm col;            /* column permutation matrices */

   L l;              /* L matrix */
   double *diag;          /* Array of pivot elements          */
   U u;              /* U matrix */

   double *work;          /* Working array: must always be left as 0! */
   double *work2;         /* Working array: must always be left as 0! */

private:
   void solveUright(double* wrk, double* vec);
   int  solveUrightEps(double* vec, int* nonz, double eps, double* rhs);
   void solveUright2(double* work1, double* vec1, double* work2, double* vec2);
   int  solveUright2eps(double* work1, double* vec1, double* work2, 
      double* vec2, int* nonz, double eps);
   void solveLright2(double* vec1, double* vec2);
   void solveUpdateRight(double* vec);
   void solveUpdateRight2(double* vec1, double* vec2);
   void solveUleft(double* work, double* vec);
   void solveUleft2(
      double* work1, double* vec1, double* work2, double* vec2);
   int solveLleft2forest(double* vec1, int* /* nonz */,
      double* vec2, double /* eps */);
   void solveLleft2(double* vec1, int* /* nonz */,
      double* vec2, double /* eps */);
   int solveLleftForest(double* vec, int* /* nonz */, double /* eps */);
   void solveLleft(double* vec);
   int solveLleftEps(double* vec, int* nonz, double eps);
   void solveUpdateLeft(double* vec);
   void solveUpdateLeft2(double* vec1, double* vec2);

public:
   void solveLright(double* vec);
   int  solveRight4update(double* vec, int* nonz, double eps, double* rhs,
      double* forest, int* forestNum, int* forestIdx);
   void solveRight(double* vec, double* rhs);
   int  solveRight2update(double* vec1, double* vec2, double* rhs1,
      double* rhs2, int* nonz, double eps, double* forest,
      int* forestNum, int* forestIdx);
   void solveRight2(double* vec1, double* vec2, double* rhs1, double* rhs2);
   void solveLeft(double* vec, double* rhs);
   int solveLeftEps(double* vec, double* rhs, int* nonz, double eps);
   int solveLeft2(double* vec1, int* nonz, double* vec2, 
      double eps, double* rhs1, double* rhs2);


};

} // namespace soplex
#endif // _CLUMEMBERS_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
