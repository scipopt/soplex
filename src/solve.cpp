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
#pragma ident "@(#) $Id: solve.cpp,v 1.4 2001/11/20 16:43:27 bzfpfend Exp $"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "clutypes.h"
#include "clumembers.h"
#include "cluprotos.h"
#include "cring.h"

namespace soplex
{


/*****************************************************************************/

static void solveUright(double* wrk, CLUFactor* fac, double* vec)
{
   int i, j, r, c;
   int *rorig, *corig;
   int *cidx, *clen, *cbeg;
   double *cval, *diag;
   double x;

   int *idx;
   double *val;
   double *work;
   work = wrk;

   rorig = fac->row.orig;
   corig = fac->col.orig;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;

   diag = fac->diag;

   for (i = fac->thedim - 1; i >= 0; --i)
   {
      r = rorig[i];
      c = corig[i];
      work[c] = x = diag[r] * vec[r];
      vec[r] = 0;
      if (x != 0.0)
      {
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
            vec[*idx++] -= x * (*val++);
      }
   }
}

static int solveUrightEps(CLUFactor* fac, double* vec,
                           int* nonz, double eps, double* rhs)
{
   int i, j, r, c, n;
   int *rorig, *corig;
   int *cidx, *clen, *cbeg;
   double *cval, *diag;
   double x, meps;

   int *idx;
   double *val;
   double *work;
   work = vec;

   rorig = fac->row.orig;
   corig = fac->col.orig;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;

   diag = fac->diag;
   meps = -eps;
   n = 0;

   for (i = fac->thedim - 1; i >= 0; --i)
   {
      r = rorig[i];
      x = diag[r] * rhs[r];
      if (x > eps || x < meps)
      {
         c = corig[i];
         work[c] = x;
         nonz[n++] = c;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
            rhs[*idx++] -= x * (*val++);
      }
   }

   return n;
}

static void solveUright2
(
   CLUFactor* fac,
   double* work1,
   double* vec1,
   double* work2,
   double* vec2
)
{
   int i, j, r, c;
   int *rorig, *corig;
   int *cidx, *clen, *cbeg;
   double *cval, *diag;
   double x1, x2;

   int* idx;
   double* val;

   rorig = fac->row.orig;
   corig = fac->col.orig;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;

   diag = fac->diag;

   for (i = fac->thedim - 1; i >= 0; --i)
   {
      r = rorig[i];
      c = corig[i];
      work1[c] = x1 = diag[r] * vec1[r];
      work2[c] = x2 = diag[r] * vec2[r];
      vec1[r] = vec2[r] = 0;
      if (x1 != 0.0 && x2 != 0.0)
      {
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
         {
            vec1[*idx] -= x1 * (*val);
            vec2[*idx++] -= x2 * (*val++);
         }
      }
      else if (x1 != 0.0)
      {
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
            vec1[*idx++] -= x1 * (*val++);
      }
      else if (x2 != 0.0)
      {
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
            vec2[*idx++] -= x2 * (*val++);
      }
   }
}

static int solveUright2eps
(
   CLUFactor* fac,
   double* work1,
   double* vec1,
   double* work2,
   double* vec2,
   int* nonz,
   double eps
)
{
   int i, j, r, c, n;
   int *rorig, *corig;
   int *cidx, *clen, *cbeg;
   int notzero1, notzero2;
   double *cval, *diag;
   double x1, x2;
   double meps;

   int* idx;
   double* val;

   rorig = fac->row.orig;
   corig = fac->col.orig;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;

   diag = fac->diag;
   meps = -eps;
   n = 0;

   for (i = fac->thedim - 1; i >= 0; --i)
   {
      c = corig[i];
      r = rorig[i];
      work1[c] = x1 = diag[r] * vec1[r];
      work2[c] = x2 = diag[r] * vec2[r];
      vec1[r] = vec2[r] = 0;
      notzero1 = (x1 < meps || x1 > eps);
      notzero2 = (x2 < meps || x2 > eps);
      if (notzero1 && notzero2)
      {
         *nonz++ = c;
         n++;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
         {
            vec1[*idx] -= x1 * (*val);
            vec2[*idx++] -= x2 * (*val++);
         }
      }
      else if (notzero1)
      {
         work2[c] = 0;
         *nonz++ = c;
         n++;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
            vec1[*idx++] -= x1 * (*val++);
      }
      else if (notzero2)
      {
         work1[c] = 0;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
            vec2[*idx++] -= x2 * (*val++);
      }
      else
      {
         work1[c] = 0;
         work2[c] = 0;
      }
   }

   return n;
}

void solveLright(CLUFactor* fac, double* vec)
{
   int i, j, k;
   int end;
   double x;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   end = fac->l.firstUpdate;
   for (i = 0; i < end; ++i)
   {
      if ((x = vec[lrow[i]]) != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            vec[*idx++] -= x * (*val++);
      }
   }

   if (fac->l.updateType)                     /* Forest-Tomlin Updates */
   {
      end = fac->l.firstUnused;
      for (; i < end; ++i)
      {
         x = 0;
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            x += vec[*idx++] * (*val++);
         vec[lrow[i]] -= x;
      }
   }
}

static void solveLright2(CLUFactor* fac, double* vec1, double* vec2)
{
   int i, j, k;
   int end;
   double x2;
   double x1;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   end = fac->l.firstUpdate;
   for (i = 0; i < end; ++i)
   {
      x1 = vec1[lrow[i]];
      x2 = vec2[lrow[i]];
      if (x1 != 0 && x2 != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            vec1[*idx] -= x1 * (*val);
            vec2[*idx++] -= x2 * (*val++);
         }
      }
      else if (x1 != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            vec1[*idx++] -= x1 * (*val++);
      }
      else if (x2 != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            vec2[*idx++] -= x2 * (*val++);
      }
   }

   if (fac->l.updateType)                     /* Forest-Tomlin Updates */
   {
      end = fac->l.firstUnused;
      for (; i < end; ++i)
      {
         x1 = 0;
         x2 = 0;
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            x1 += vec1[*idx] * (*val);
            x2 += vec2[*idx++] * (*val++);
         }
         vec1[lrow[i]] -= x1;
         vec2[lrow[i]] -= x2;
      }
   }
}

static void solveUpdateRight(CLUFactor* fac, double* vec)
{
   int i, j, k;
   int end;
   double x;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   assert(!fac->l.updateType);               /* no Forest-Tomlin Updates */

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   end = fac->l.firstUnused;
   for (i = fac->l.firstUpdate; i < end; ++i)
   {
      if ((x = vec[lrow[i]]) != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            vec[*idx++] -= x * (*val++);
      }
   }
}

static void solveUpdateRight2(CLUFactor* fac, double* vec1, double* vec2)
{
   int i, j, k;
   int end;
   double x1, x2;
   double *lval;
   int *lrow, *lidx;
   int *lbeg;

   int* idx;
   double* val;

   assert(!fac->l.updateType);               /* no Forest-Tomlin Updates */

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   end = fac->l.firstUnused;
   for (i = fac->l.firstUpdate; i < end; ++i)
   {
      x1 = vec1[lrow[i]];
      x2 = vec2[lrow[i]];
      if (x1 != 0.0 && x2 != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            vec1[*idx] -= x1 * (*val);
            vec2[*idx++] -= x2 * (*val++);
         }
      }
      else if (x1 != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            vec1[*idx++] -= x1 * (*val++);
      }
      else if (x2 != 0.0)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
            vec2[*idx++] -= x2 * (*val++);
      }
   }
}

int solveRight4update(CLUFactor* fac, double* vec, int* nonz, double eps,
                       double* rhs,
                       double* forest, int* forestNum, int* forestIdx)
{
   solveLright(fac, rhs);

   if (forest)
   {
      double* r = rhs;
      int n = 0;
      int i = 0;
      int e = fac->thedim;
      int* idx = forestIdx;
      for (; i < e;)
      {
         idx[n] = i++;
         n += ((*forest++ = *r++) != 0);
      }
      *forestNum = n;
   }

   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      solveUright(vec, fac, rhs);
      solveUpdateRight(fac, vec);
      return 0;
   }
   else
      return solveUrightEps(fac, vec, nonz, eps, rhs);
}

void solveRight(CLUFactor* fac, double* vec, double* rhs)
{
   solveLright(fac, rhs);
   solveUright(vec, fac, rhs);
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
      solveUpdateRight(fac, vec);
}


int solveRight2update(CLUFactor* fac,
                       double* vec1,
                       double* vec2,
                       double* rhs1,
                       double* rhs2,
                       int* nonz,
                       double eps,
                       double* forest,
                       int* forestNum,
                       int* forestIdx)
{
   solveLright2(fac, rhs1, rhs2);

   if (forest)
   {
      double* r = rhs1;
      int n = 0;
      int i = 0;
      int e = fac->thedim;
      int* idx = forestIdx;
      for (; i < e;)
      {
         idx[n] = i++;
         n += ((*forest++ = *r++) != 0);
      }
      *forestNum = n;
   }

   if (! fac->l.updateType)           /* no Forest-Tomlin Updates */
   {
      solveUright2(fac, vec1, rhs1, vec2, rhs2);
      solveUpdateRight2(fac, vec1, vec2);
      return 0;
   }
   else
      return solveUright2eps(fac, vec1, rhs1, vec2, rhs2, nonz, eps);
}

void solveRight2
(
   CLUFactor* fac,
   double* vec1,
   double* vec2,
   double* rhs1,
   double* rhs2
)
{
   solveLright2(fac, rhs1, rhs2);
   if (fac->l.updateType)             /* Forest-Tomlin Updates */
      solveUright2(fac, vec1, rhs1, vec2, rhs2);
   else
   {
      solveUright2(fac, vec1, rhs1, vec2, rhs2);
      solveUpdateRight2(fac, vec1, vec2);
   }
}

/*****************************************************************************/

static void solveUleft(double* work, double* vec, CLUFactor* fac)
{
   double x;
   int i, k, l, r, c;
   int end;
   int *rorig, *corig;
   int *ridx, *rlen, *rbeg, *idx;
   double *rval, *diag, *val;

   rorig = fac->row.orig;
   corig = fac->col.orig;

   ridx = fac->u.row.idx;
   rval = fac->u.row.val;
   rlen = fac->u.row.len;
   rbeg = fac->u.row.start;

   diag = fac->diag;

   end = fac->thedim;
   for (i = 0; i < end; ++i)
   {
      c = corig[i];
      r = rorig[i];
      x = vec[c];
      vec[c] = 0;
      if (x)
      {
         x *= diag[r];
         work[r] = x;
         k = rbeg[r];
         idx = &ridx[k];
         val = &rval[k];
         for (l = rlen[r]; l; --l)
            vec[*idx++] -= x * (*val++);
      }
   }
}


static void solveUleft2
(
   double* work1,
   double* vec1,
   double* work2,
   double* vec2,
   CLUFactor* fac
)
{
   double x1;
   double x2;
   int i, k, l, r, c;
   int end;
   int *rorig, *corig;
   int *ridx, *rlen, *rbeg, *idx;
   double *rval, *diag, *val;

   rorig = fac->row.orig;
   corig = fac->col.orig;

   ridx = fac->u.row.idx;
   rval = fac->u.row.val;
   rlen = fac->u.row.len;
   rbeg = fac->u.row.start;

   diag = fac->diag;

   end = fac->thedim;
   for (i = 0; i < end; ++i)
   {
      c = corig[i];
      r = rorig[i];
      x1 = vec1[c];
      x2 = vec2[c];
      if (x1 && x2)
      {
         x1 *= diag[r];
         x2 *= diag[r];
         work1[r] = x1;
         work2[r] = x2;
         k = rbeg[r];
         idx = &ridx[k];
         val = &rval[k];
         for (l = rlen[r]; l; --l)
         {
            vec1[*idx] -= x1 * (*val);
            vec2[*idx++] -= x2 * (*val++);
         }
      }
      else if (x1)
      {
         x1 *= diag[r];
         work1[r] = x1;
         k = rbeg[r];
         idx = &ridx[k];
         val = &rval[k];
         for (l = rlen[r]; l; --l)
            vec1[*idx++] -= x1 * (*val++);
      }
      else if (x2)
      {
         x2 *= diag[r];
         work2[r] = x2;
         k = rbeg[r];
         idx = &ridx[k];
         val = &rval[k];
         for (l = rlen[r]; l; --l)
            vec2[*idx++] -= x2 * (*val++);
      }
   }
}

static int solveLleft2forest
(
   double* vec1,
   int* /* nonz */,
   double* vec2,
   double /* eps */,
   CLUFactor* fac
)
{
   int i;
   int j;
   int k;
   int end;
   double x1, x2;
   double *lval, *val;
   int *lidx, *idx, *lrow;
   int *lbeg;
   int *rorig;

   rorig = fac->row.orig;
   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   end = fac->l.firstUpdate;
   for (i = fac->l.firstUnused - 1; i >= end; --i)
   {
      j = lrow[i];
      x1 = vec1[j];
      x2 = vec2[j];
      if (x1 != 0.0)
      {
         if (x2 != 0.0)
         {
            k = lbeg[i];
            val = &lval[k];
            idx = &lidx[k];
            for (j = lbeg[i + 1]; j > k; --j)
            {
               vec1[*idx] -= x1 * (*val);
               vec2[*idx++] -= x2 * (*val++);
            }
         }
         else
         {
            k = lbeg[i];
            val = &lval[k];
            idx = &lidx[k];
            for (j = lbeg[i + 1]; j > k; --j)
               vec1[*idx++] -= x1 * (*val++);
         }
      }
      else if (x2 != 0.0)
      {
         k = lbeg[i];
         val = &lval[k];
         idx = &lidx[k];
         for (j = lbeg[i + 1]; j > k; --j)
            vec2[*idx++] -= x2 * (*val++);
      }
   }
   return 0;
}

static void solveLleft2
(
   double* vec1,
   int* /* nonz */,
   double* vec2,
   double /* eps */,
   CLUFactor* fac
)
{
   int i, j, k, r;
   int x1not0, x2not0;
   double x1, x2;

   double *rval, *lval, *val;
   int *ridx, *lidx, *idx, *lrow;
   int *rbeg, *lbeg;
   int *rorig;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   ridx = fac->l.ridx;
   rval = fac->l.rval;
   rbeg = fac->l.rbeg;
   rorig = fac->l.rorig;

#ifndef WITH_L_ROWS
   i = fac->l.firstUpdate - 1;
   for (; i >= 0; --i)
   {
      k = lbeg[i];
      val = &lval[k];
      idx = &lidx[k];
      x1 = 0;
      x2 = 0;
      for (j = lbeg[i + 1]; j > k; --j)
      {
         x1 += vec1[*idx] * (*val);
         x2 += vec2[*idx++] * (*val++);
      }
      vec1[lrow[i]] -= x1;
      vec2[lrow[i]] -= x2;
   }
#else
for (i = fac->thedim; i--;)
   {
      r = rorig[i];
      x1 = vec1[r];
      x2 = vec2[r];
      x1not0 = (x1 != 0);
      x2not0 = (x2 != 0);

      if (x1not0 && x2not0)
      {
         k = rbeg[r];
         j = rbeg[r + 1] - k;
         val = &rval[k];
         idx = &ridx[k];
         while (j-- > 0)
         {
            assert(fac->row.perm[*idx] < i);
            vec1[*idx] -= x1 * *val;
            vec2[*idx++] -= x2 * *val++;
         }
      }
      else if (x1not0)
      {
         k = rbeg[r];
         j = rbeg[r + 1] - k;
         val = &rval[k];
         idx = &ridx[k];
         while (j-- > 0)
         {
            assert(fac->row.perm[*idx] < i);
            vec1[*idx++] -= x1 * *val++;
         }
      }
      else if (x2not0)
      {
         k = rbeg[r];
         j = rbeg[r + 1] - k;
         val = &rval[k];
         idx = &ridx[k];
         while (j-- > 0)
         {
            assert(fac->row.perm[*idx] < i);
            vec2[*idx++] -= x2 * *val++;
         }
      }
   }
#endif
}

static int solveLleftForest(double* vec, int* /* nonz */, double /* eps */, CLUFactor* fac)
{
   int i, j, k, end;
   double x;
   double *val, *lval;
   int *idx, *lidx, *lrow, *lbeg;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   end = fac->l.firstUpdate;
   for (i = fac->l.firstUnused - 1; i >= end; --i)
   {
      if ((x = vec[lrow[i]]) != 0.0)
      {
         k = lbeg[i];
         val = &lval[k];
         idx = &lidx[k];
         for (j = lbeg[i + 1]; j > k; --j)
            vec[*idx++] -= x * (*val++);
      }
   }

   return 0;
}

static void solveLleft(double* vec, CLUFactor* fac)
{
   int i, j, k;
   int r;
   double x;
   double *rval, *lval, *val;
   int *ridx, *lidx, *idx, *lrow;
   int *rbeg, *lbeg;
   int* rorig;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   ridx = fac->l.ridx;
   rval = fac->l.rval;
   rbeg = fac->l.rbeg;
   rorig = fac->l.rorig;

#ifndef WITH_L_ROWS
   i = fac->l.firstUpdate - 1;
   for (; i >= 0; --i)
   {
      k = lbeg[i];
      val = &lval[k];
      idx = &lidx[k];
      x = 0;
      for (j = lbeg[i + 1]; j > k; --j)
         x += vec[*idx++] * (*val++);
      vec[lrow[i]] -= x;
   }
#else
   for (i = fac->thedim; i--;)
   {
      r = rorig[i];
      x = vec[r];
      if (x != 0.0)
      {
         k = rbeg[r];
         j = rbeg[r + 1] - k;
         val = &rval[k];
         idx = &ridx[k];
         while (j-- > 0)
         {
            assert(fac->l.rperm[*idx] < i);
            vec[*idx++] -= x * *val++;
         }
      }
   }
#endif
}

static int solveLleftEps(double* vec, CLUFactor* fac, int* nonz, double eps)
{
   int i, j, k, n;
   int r;
   double x, meps;
   double *rval, *lval, *val;
   int *ridx, *lidx, *idx, *lrow;
   int *rbeg, *lbeg;
   int* rorig;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   ridx = fac->l.ridx;
   rval = fac->l.rval;
   rbeg = fac->l.rbeg;
   rorig = fac->l.rorig;
   meps = -eps;
   n = 0;

   i = fac->l.firstUpdate - 1;
#ifndef WITH_L_ROWS
   for (; i >= 0; --i)
   {
      k = lbeg[i];
      val = &lval[k];
      idx = &lidx[k];
      x = 0;
      for (j = lbeg[i + 1]; j > k; --j)
         x += vec[*idx++] * (*val++);
      vec[lrow[i]] -= x;
   }
#else
   for (i = fac->thedim; i--;)
   {
      r = rorig[i];
      x = vec[r];
      if (x < meps || x > eps)
      {
         *nonz++ = r;
         n++;
         k = rbeg[r];
         j = rbeg[r + 1] - k;
         val = &rval[k];
         idx = &ridx[k];
         while (j-- > 0)
         {
            assert(fac->row.perm[*idx] < i);
            vec[*idx++] -= x * *val++;
         }
      }
      else
         vec[r] = 0;
   }
#endif

   return n;
}

static void solveUpdateLeft(double* vec, CLUFactor* fac)
{
   int i, j, k, end;
   double x;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   assert(!fac->l.updateType);               /* Forest-Tomlin Updates */

   end = fac->l.firstUpdate;
   for (i = fac->l.firstUnused - 1; i >= end; --i)
   {
      k = lbeg[i];
      val = &lval[k];
      idx = &lidx[k];
      x = 0;
      for (j = lbeg[i + 1]; j > k; --j)
         x += vec[*idx++] * (*val++);
      vec[lrow[i]] -= x;
   }
}

static void solveUpdateLeft2(double* vec1, double* vec2, CLUFactor* fac)
{
   int i, j, k, end;
   double x1, x2;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;

   assert(!fac->l.updateType);               /* Forest-Tomlin Updates */

   end = fac->l.firstUpdate;
   for (i = fac->l.firstUnused - 1; i >= end; --i)
   {
      k = lbeg[i];
      val = &lval[k];
      idx = &lidx[k];
      x1 = 0;
      x2 = 0;
      for (j = lbeg[i + 1]; j > k; --j)
      {
         x1 += vec1[*idx] * (*val);
         x2 += vec2[*idx++] * (*val++);
      }
      vec1[lrow[i]] -= x1;
      vec2[lrow[i]] -= x2;
   }
}

void solveLeft(double* vec, CLUFactor* fac, double* rhs)
{
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      solveUpdateLeft(rhs, fac);
      solveUleft(vec, rhs, fac);
      solveLleft(vec, fac);
   }
   else
   {
      solveUleft(vec, rhs, fac);
      solveLleftForest(vec, 0, 0, fac);
      solveLleft(vec, fac);
   }
}

int solveLeftEps(double* vec, CLUFactor* fac, double* rhs, int* nonz, double eps)
{
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      solveUpdateLeft(rhs, fac);
      solveUleft(vec, rhs, fac);
      return solveLleftEps(vec, fac, nonz, eps);
   }
   else
   {
      solveUleft(vec, rhs, fac);
      solveLleftForest(vec, nonz, eps, fac);
      return solveLleftEps(vec, fac, nonz, eps);
   }
}

int solveLeft2
(
   CLUFactor* fac,
   double* vec1,
   int* nonz,
   double* vec2,
   double eps,
   double* rhs1,
   double* rhs2
)
{
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      solveUpdateLeft2(rhs1, rhs2, fac);
      solveUleft2(vec1, rhs1, vec2, rhs2, fac);
      solveLleft2(vec1, nonz, vec2, eps, fac);
      return 0;
   }
   else
   {
      solveUleft2(vec1, rhs1, vec2, rhs2, fac);
      solveLleft2forest(vec1, nonz, vec2, eps, fac);
      solveLleft2(vec1, nonz, vec2, eps, fac);
      return 0;
   }
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
