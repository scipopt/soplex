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
#pragma ident "@(#) $Id: clumembers.h,v 1.7 2001/11/30 22:14:59 bzfkocht Exp $"

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

   // From solve.cpp
   ///
   void solveUright(double* wrk, double* vec);
   ///
   int  solveUrightEps(double* vec, int* nonz, double eps, double* rhs);
   ///
   void solveUright2(double* work1, double* vec1, double* work2, double* vec2);
   ///
   int  solveUright2eps(double* work1, double* vec1, double* work2, 
      double* vec2, int* nonz, double eps);
   ///
   void solveLright2(double* vec1, double* vec2);
   ///
   void solveUpdateRight(double* vec);
   ///
   void solveUpdateRight2(double* vec1, double* vec2);
   ///
   void solveUleft(double* work, double* vec);
   ///
   void solveUleft2(
      double* work1, double* vec1, double* work2, double* vec2);
   ///
   int solveLleft2forest(double* vec1, int* /* nonz */,
      double* vec2, double /* eps */);
   ///
   void solveLleft2(double* vec1, int* /* nonz */,
      double* vec2, double /* eps */);
   ///
   int solveLleftForest(double* vec, int* /* nonz */, double /* eps */);
   ///
   void solveLleft(double* vec);
   ///
   int solveLleftEps(double* vec, int* nonz, double eps);
   ///
   void solveUpdateLeft(double* vec);
   ///
   void solveUpdateLeft2(double* vec1, double* vec2);

   // From vsolve.cpp 
   ///
   int vSolveLright(double* vec, int* ridx, int rn, double eps);
      void vSolveLright2(double* vec, int* ridx, int* rnptr, double eps,
         double* vec2, int* ridx2, int* rn2ptr, double eps2);
   ///
   int vSolveUright(double* vec, int* vidx,
      double* rhs, int* ridx, int rn, double eps);
   ///
   void vSolveUrightNoNZ(double* vec,
      double* rhs, int* ridx, int rn, double eps);
   ///
   int vSolveUright2(double* vec, int* vidx, double* rhs, 
      int* ridx, int rn, double eps,
      double* vec2, double* rhs2, int* ridx2, int rn2, double eps2);
   ///
   int vSolveUpdateRight(double* vec, int* ridx, int n, double eps);
   ///
   void vSolveUpdateRightNoNZ(double* vec, double /*eps*/);
   ///
   int solveUleft(double eps, double* vec, int* vecidx,
      double* rhs, int* rhsidx, int rhsn);
   ///
   void solveUleftNoNZ(double eps, double* vec,
      double* rhs, int* rhsidx, int rhsn);
   ///
   int solveLleftForest(double eps, double* vec, int* nonz, int n);
   ///
   void solveLleftForestNoNZ(double* vec);
   ///
   int solveLleft(double eps, double* vec, int* nonz, int rn);
   ///
   void solveLleftNoNZ(double* vec);
   ///
   int solveUpdateLeft(double eps, double* vec, int* nonz, int n);

public:
   // From solve.cpp 
   ///
   void solveLright(double* vec);
   ///
   int  solveRight4update(double* vec, int* nonz, double eps, double* rhs,
      double* forest, int* forestNum, int* forestIdx);
   ///
   void solveRight(double* vec, double* rhs);
   ///
   int  solveRight2update(double* vec1, double* vec2, double* rhs1,
      double* rhs2, int* nonz, double eps, double* forest,
      int* forestNum, int* forestIdx);
   ///
   void solveRight2(double* vec1, double* vec2, double* rhs1, double* rhs2);
   ///
   void solveLeft(double* vec, double* rhs);
   ///
   int solveLeftEps(double* vec, double* rhs, int* nonz, double eps);
   ///
   int solveLeft2(double* vec1, int* nonz, double* vec2, 
      double eps, double* rhs1, double* rhs2);

   // From vsolve.cpp: Very sparse solution methods.
   ///
   int vSolveRight4update(double eps, double* vec, int* idx, /* result       */
                        double* rhs, int* ridx, int rn,      /* rhs & Forest */
                        double* forest, int* forestNum, int* forestIdx);
   ///
   int vSolveRight4update2(double eps,
      double* vec, int* idx,                      /* result1 */
      double* rhs, int* ridx, int rn,             /* rhs1    */
      double* vec2, double eps2,                  /* result2 */
      double* rhs2, int* ridx2, int rn2,          /* rhs2    */
      double* forest, int* forestNum, int* forestIdx);
   ///
   void vSolveRightNoNZ(double* vec2, double eps2,              /* result2 */
      double* rhs2, int* ridx2, int rn2);   /* rhs2    */
   ///
   int vSolveLeft(double eps,
      double* vec, int* idx,                               /* result */
      double* rhs, int* ridx, int rn);                   /* rhs    */
   ///
   void vSolveLeftNoNZ(double eps,
      double* vec,                             /* result */
      double* rhs, int* ridx, int rn);       /* rhs    */
   ///
   int vSolveLeft2(double eps,
      double* vec, int* idx,                              /* result */
      double* rhs, int* ridx, int rn,                     /* rhs    */
      double* vec2,                                       /* result2 */
      double* rhs2, int* ridx2, int rn2);               /* rhs2    */

#warning MAKE Pring and Temp private when factor.cpp is finished
private:
public: 
   class Pring   /* Pivot ring */
   {
   public:
      Pring(): next(NULL), prev(NULL){}      
      Pring *next;
      Pring *prev;
      int idx;            /* index of pivot row */
      int pos;            /* position of pivot column in row */
      int mkwtz;          /* markowitz number of pivot */
   private:
      Pring(const Pring&);
      Pring& operator= (const Pring&);
   };
   class Temp  /* Temporary data structures. */
   {
   public: 
      Temp(int p_dim);
      ~Temp();
      int*    s_mark;
      double* s_max;           /* maximum absolute value per row (or -1) */
      int*    s_cact;          /* lengths of columns of active submatrix */
   private:
      Temp( const Temp& );
      Temp& operator= ( const Temp& );
   };
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
