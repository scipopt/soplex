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
#pragma ident "@(#) $Id: clufactor.h,v 1.8 2001/12/25 14:32:37 bzfkocht Exp $"

#ifndef _CLUFACTOR_H_
#define _CLUFACTOR_H_

#include "slinsolver.h"
#include "svector.h"

#define WITH_L_ROWS 1

namespace soplex
{

/**todo I think that CLUFactor should be a descendent of SLinSolver.
        This has to be muddled out with SLUFactor also.
        Probably CLUFactor should be inlined in SLUFactor.
 */
class CLUFactor
{
public:
   /** Double linked ring structure for garbage collection of column or
    *  row file in working matrix.
    */
   struct Dring
   {  
      Dring* next; 
      Dring* prev;
      int    idx;
   };

   /// Pivot Ring
   class Pring
   {
   public:
      Pring* next;
      Pring* prev;
      int   idx;            ///< index of pivot row
      int   pos;            ///< position of pivot column in row 
      int   mkwtz;          ///< markowitz number of pivot 

      Pring() : next(NULL), prev(NULL) 
      {}      

   private:
      Pring(const Pring&);
      Pring& operator= (const Pring&);
   };
protected:
   /// Temporary data structures.
   class Temp 
   {
   public: 
      int*    s_mark;
      double* s_max;        ///< maximum absolute value per row (or -1) 
      int*    s_cact;       ///< lengths of columns of active submatrix 
      int     stage;           
      Pring   pivots;       ///< ring of selected pivot rows 
      Pring*  pivot_col;    ///< column index handlers for double linked list 
      Pring*  pivot_colNZ;  ///< lists for columns to number of nonzeros      
      Pring*  pivot_row;    ///< row index handlers for double linked list 
      Pring*  pivot_rowNZ;  ///< lists for rows to number of nonzeros

      Temp();
      ~Temp();
      void init(int p_dim);
      void clear();

   private:
      Temp( const Temp& );
      Temp& operator= ( const Temp& );
   };

   /// Data structures for saving the row and column permutations.
   struct Perm
   {
      int* orig;          ///< orig[p] original index from p 
      int* perm;          ///< perm[i] permuted index from i 
   };

   /// Data structures for saving the working matrix and U factor.
   struct U
   {
      struct Row
      {
         Dring list;         /* Double linked ringlist of vector
                                indices in the order they appear
                                in the row file */
         Dring *elem;        /* Array of ring elements.            */
         int size;           /* size of arrays val and idx         */
         int used;           /* used entries of arrays idx and val */
         double *val;        /* hold nonzero values                */
         int *idx;           /* hold nonzero indices               */
         int *start;         /* starting positions in val and idx  */
         int *len;           /* used nonzeros per row vectors      */
         int *max;           /* maximum available nonzeros per row:
                                start[i] + max[i] == start[elem[i].next->idx] 
                                len[i] <= max[i]. */
      } row;
      struct Col
      {
         Dring list;         /* Double linked ringlist of vector
                                indices in the order they appear
                                in the column file */
         Dring *elem;        /* Array of ring elements.            */
         int size;           /* size of array idx                  */
         int used;           /* used entries of array idx          */
         int *idx;           /* hold nonzero indices               */
         double *val;        /* hold nonzero values: this is only initialized
                                in the end of the factorization with DEFAULT
                                updates. */
         int *start;         /* starting positions in val and idx  */
         int *len;           /* used nonzeros per column vector    */
         int *max;           /* maximum available nonzeros per colunn:
                                start[i] + max[i] == start[elem[i].next->idx] 
                                len[i] <= max[i]. */
      } col;
      int lastColSing;       /* stage of last eliminated column singleton */
      int lastRowSing;       /* stage of last eliminated row singleton */
   };


   /// Data structures for saving the working matrix and U factor.
   struct L
   {
      int size;           /* size of arrays val and idx        */
      double *val;        /* values of L vectors               */
      int *idx;           /* indices of L vectors              */
      int startSize;      /* size of array start               */
      int firstUpdate;    /* number of first update L vector   */
      int firstUnused;    /* number of first unused L vector   */
      int *start;         /* starting positions in val and idx */
      int *row;           /* column indices of L vectors       */
      int updateType;     /* type of updates to be used.       */

      /* The following arrays have length |firstUpdate|, since they keep
       * rows of the L-vectors occuring during the factorization (without
       * updates), only:
       */
      double *rval;       /* values of rows of L               */
      int *ridx;          /* indices of rows of L              */
      int *rbeg;          /* start of rows in rval and ridx    */
      int *rorig;         /* original row permutation          */
      int *rperm;         /* original row permutation          */
   };

   SLinSolver::Status stat;   ///< Status indicator.

   int     thedim;            ///< dimension of factorized matrix   
   int     nzCnt;             ///< number of nonzeros in U      
   double  initMaxabs;        ///< maximum abs number in initail Matrix 
   double  maxabs;            ///< maximum abs number in L and U        

   double  rowMemMult;        ///< factor of minimum Memory * #of nonzeros 
   double  colMemMult;        ///< factor of minimum Memory * #of nonzeros 
   double  lMemMult;          ///< factor of minimum Memory * #of nonzeros 

   Perm    row;               ///< row permutation matrices 
   Perm    col;               ///< column permutation matrices 

   L       l;                 ///< L matrix 
   double* diag;              ///< Array of pivot elements          
   U       u;                 ///< U matrix 

   double* work;              ///< Working array: must always be left as 0! 
   double* work2;             ///< Working array: must always be left as 0! 

private:
   Temp    temp;              ///< Temporary storage

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

   // from forest.cpp
   ///
   void forestPackColumns();
   ///
   void forestMinColMem(int size);
   ///
   void forestReMaxCol(int col, int len);

   // from factor.cpp
   ///
   void initPerm();
   ///
   void initFactorMatrix(SVector** vec, const double eps );
   ///
   void minLMem(int size);
   ///
   void setPivot(const int p_stage, const int p_col, 
                 const int p_row, const double val);
   ///
   void colSingletons();
   ///
   void rowSingletons();

   ///
   void initFactorRings();
   ///
   void freeFactorRings();
      
   ///
   int setupColVals();
   ///
   void setupRowVals();

   ///
   void eliminateRowSingletons();
   ///
   void eliminateColSingletons();
   ///
   void selectPivots(double threshold);
   ///
   int updateRow(int r, int lv, int prow, int pcol, double pval, double eps);

   ///
   void eliminatePivot(int prow, int pos, double eps);
   ///
   void eliminateNucleus(const double eps, const double threshold);
   ///
   void minRowMem(int size);
   ///
   void minColMem(int size);
   ///
   void remaxCol(int p_col, int len);
   ///
   void packRows();
   ///
   void packColumns();
   ///
   void remaxRow(int p_row, int len);
   ///
   int makeLvec(int p_len, int p_row);

   /// no copy construtor.
   CLUFactor(const CLUFactor&);
   /// no assignment operator.
   CLUFactor& operator=(const CLUFactor&);

protected:
   /// default construtor. 
   /** Since there is no sense in constructing a CLUFactor object
    *  per se, this is protected.
    */
   CLUFactor() {}

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

   // from forest.cpp
   void forestUpdate(int col, double* work, int num, int *nonz);

   // from update.cpp
   void update(int p_col, double* p_work, const int* p_idx, int num);
   void updateNoClear(
      int p_col, const double* p_work, const int* p_idx, int num);


   // from factor.cpp
   ///
   void factor(SVector** vec,     ///< Array of column vector pointers  
      double threshold,           ///< pivoting threshold                
      double eps);                ///< epsilon for zero detection        

   ///
   void dump() const;
   ///
   bool isConsistent() const;
};

} // namespace soplex
#endif // _CLUFACTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
