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
#pragma ident "@(#) $Id: vsolve.cpp,v 1.3 2001/11/08 14:27:29 bzfkocht Exp $"


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

#define ZERO    1e-100


static const double verySparseFactor4right = 0.2;
static const double verySparseFactor4left  = 0.1;

/*****************************************************************************/
static
void enQueueMax(int* heap, int* size, int elem)
{
   int i, j;

   j = (*size)++;
   while (j > 0)
   {
      i = (j - 1) / 2;
      if (elem > heap[i])
      {
         heap[j] = heap[i];
         j = i;
      }
      else
         break;
   }
   heap[j] = elem;

#ifdef  DEBUG
   for (i = 1; i < *size; ++i)
      for (j = 0; j < i; ++j)
         if (heap[i] == heap[j])
            printf("ERROR\n");
#endif  /* DEBUG */
}

static
int deQueueMax(int* heap, int* size)
{
   int e, elem;
   int i, j, s;
   int e1, e2;

   elem = *heap;
   e = heap[s = --(*size)];
   --s;
   for (j = 0, i = 1; i < s; i = 2 * j + 1)
   {
      e1 = heap[i];
      e2 = heap[i + 1];
      if (e1 > e2)
      {
         if (e < e1)
         {
            heap[j] = e1;
            j = i;
         }
         else
         {
            heap[j] = e;
            return elem;
         }
      }
      else
      {
         if (e < e2)
         {
            heap[j] = e2;
            j = i + 1;
         }
         else
         {
            heap[j] = e;
            return elem;
         }
      }
   }

   if (i < *size && e < heap[i])
   {
      heap[j] = heap[i];
      j = i;
   }

   heap[j] = e;
   return elem;
}

static
void enQueueMin(int* heap, int* size, int elem)
{
   int i, j;

   j = (*size)++;
   while (j > 0)
   {
      i = (j - 1) / 2;
      if (elem < heap[i])
      {
         heap[j] = heap[i];
         j = i;
      }
      else
         break;
   }
   heap[j] = elem;

#ifdef  DEBUG
   for (i = 1; i < *size; ++i)
      for (j = 0; j < i; ++j)
         if (heap[i] == heap[j])
            printf("ERROR\n");
#endif  /* DEBUG */
}

static
int deQueueMin(int* heap, int* size)
{
   int e, elem;
   int i, j, s;
   int e1, e2;

   elem = *heap;
   e = heap[s = --(*size)];
   --s;
   for (j = 0, i = 1; i < s; i = 2 * j + 1)
   {
      e1 = heap[i];
      e2 = heap[i + 1];
      if (e1 < e2)
      {
         if (e > e1)
         {
            heap[j] = e1;
            j = i;
         }
         else
         {
            heap[j] = e;
            return elem;
         }
      }
      else
      {
         if (e > e2)
         {
            heap[j] = e2;
            j = i + 1;
         }
         else
         {
            heap[j] = e;
            return elem;
         }
      }
   }

   if (i < *size && e > heap[i])
   {
      heap[j] = heap[i];
      j = i;
   }

   heap[j] = e;
   return elem;
}

/*****************************************************************************/
static
int vSolveLright(CLUFactor* fac, double* vec, int* ridx, int rn, double eps)
{
   int i, j, k, n;
   int end;
   double x;
   double meps;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   meps = -eps;

   end = fac->l.firstUpdate;
   for (i = 0; i < end; ++i)
   {
      x = vec[lrow[i]];
      if (x > eps || x < meps)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            ridx[rn] = n = *idx++;
            rn += (vec[n] == 0);
            vec[n] -= x * (*val++);
            vec[n] += ZERO * (vec[n] == 0);
         }
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
         ridx[rn] = j = lrow[i];
         rn += (vec[j] == 0);
         vec[j] -= x;
         vec[j] += ZERO * (vec[j] == 0);
      }
   }

   return rn;
}

static
void vSolveLright2(CLUFactor* fac,
                    double* vec, int* ridx, int* rnptr, double eps,
                    double* vec2, int* ridx2, int* rn2ptr, double eps2)
{
   int i, j, k, n;
   int end;
   double x, y;
   double x2, y2;
   double meps2, meps;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   int rn = *rnptr;
   int rn2 = *rn2ptr;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   meps = -eps;
   meps2 = -eps2;

   end = fac->l.firstUpdate;
   for (i = 0; i < end; ++i)
   {
      j = lrow[i];
      x2 = vec2[j];
      x = vec[j];
      if (x > eps || x < meps)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         if (x2 > eps2 || x2 < meps2)
         {
            for (j = lbeg[i + 1]; j > k; --j)
            {
               ridx[rn] = ridx2[rn2] = n = *idx++;
               y = vec[n];
               y2 = vec2[n];
               rn += (y == 0);
               rn2 += (y2 == 0);
               y -= x * (*val);
               y2 -= x2 * (*val++);
               vec[n] = y + ZERO * (y == 0);
               vec2[n] = y2 + ZERO * (y2 == 0);
            }
         }
         else
         {
            for (j = lbeg[i + 1]; j > k; --j)
            {
               ridx[rn] = n = *idx++;
               y = vec[n];
               rn += (y == 0);
               y -= x * (*val++);
               vec[n] = y + ZERO * (y == 0);
            }
         }
      }
      else if (x2 > eps2 || x2 < meps2)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            ridx2[rn2] = n = *idx++;
            y2 = vec2[n];
            rn2 += (y2 == 0);
            y2 -= x2 * (*val++);
            vec2[n] = y2 + ZERO * (y2 == 0);
         }
      }
   }

   if (fac->l.updateType)                     /* Forest-Tomlin Updates */
   {
      end = fac->l.firstUnused;
      for (; i < end; ++i)
      {
         x = x2 = 0;
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            x += vec[*idx] * (*val);
            x2 += vec2[*idx++] * (*val++);
         }
         ridx[rn] = ridx2[rn2] = j = lrow[i];
         rn += (vec[j] == 0);
         rn2 += (vec2[j] == 0);
         vec[j] -= x;
         vec2[j] -= x2;
         vec[j] += ZERO * (vec[j] == 0);
         vec2[j] += ZERO * (vec2[j] == 0);
      }
   }

   *rnptr = rn;
   *rn2ptr = rn2;
}

static
int vSolveUright(CLUFactor* fac, double* vec, int* vidx,
                  double* rhs, int* ridx, int rn, double eps)
{
   int i, j, k, r, c, n;
   int dim;
   int *rorig, *corig;
   int *rperm, *cperm;
   int *cidx, *clen, *cbeg;
   double *cval, *diag;
   double x, y, meps;
   int lastRowSing, lastColSing;

   int *idx;
   double *val;
   double *work;
   work = vec;

   rorig = fac->row.orig;
   corig = fac->col.orig;
   rperm = fac->row.perm;
   cperm = fac->col.perm;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;
   diag = fac->diag;
   dim = fac->thedim;

   lastRowSing = fac->u.lastRowSing;
   lastColSing = fac->u.lastColSing;

   meps = -eps;
   n = 0;

   while (rn > 0)
   {
      /*      Find nonzero with highest permuted row index and setup i and r
       */
      i = deQueueMax(ridx, &rn);
      r = rorig[i];

      x = diag[r] * rhs[r];
      rhs[r] = 0;
      if (x > eps || x < meps)
      {
         c = corig[i];
         vidx[n++] = c;
         work[c] = x;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
         {
            k = *idx++;
            y = rhs[k];
            if (y == 0)
            {
               y = -x * (*val++);
               if (y < meps || y > eps)
               {
                  rhs[k] = y;
                  enQueueMax(ridx, &rn, rperm[k]);
               }
            }
            else
            {
               y -= x * (*val++);
               y += ZERO * (y == 0);
               rhs[k] = y;
            }
         }

         if (rn > i*verySparseFactor4right)
         {                                   /* continue with dense case */
            for (i = *ridx; i >= 0; --i)
            {
               r = rorig[i];
               x = diag[r] * rhs[r];
               rhs[r] = 0;
               if (x > eps || x < meps)
               {
                  c = corig[i];
                  vidx[n++] = c;
                  work[c] = x;
                  val = &cval[cbeg[c]];
                  idx = &cidx[cbeg[c]];
                  j = clen[c];
                  while (j-- > 0)
                     rhs[*idx++] -= x * (*val++);
               }
            }
            break;
         }
      }
   }

   return n;
}

static
void vSolveUrightNoNZ(CLUFactor* fac, double* vec,
                       double* rhs, int* ridx, int rn, double eps)
{
   int i, j, k, r, c;
   int dim;
   int *rorig, *corig;
   int *rperm, *cperm;
   int *cidx, *clen, *cbeg;
   double *cval, *diag;
   double x, y, meps;
   int lastRowSing, lastColSing;

   int *idx;
   double *val;
   double *work;
   work = vec;

   rorig = fac->row.orig;
   corig = fac->col.orig;
   rperm = fac->row.perm;
   cperm = fac->col.perm;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;
   diag = fac->diag;
   dim = fac->thedim;

   lastRowSing = fac->u.lastRowSing;
   lastColSing = fac->u.lastColSing;

   meps = -eps;
   while (rn > 0)
   {
      if (rn > *ridx * verySparseFactor4right)
      {                                       /* continue with dense case */
         for (i = *ridx; i >= 0; --i)
         {
            r = rorig[i];
            x = diag[r] * rhs[r];
            rhs[r] = 0;
            if (x > eps || x < meps)
            {
               c = corig[i];
               work[c] = x;
               val = &cval[cbeg[c]];
               idx = &cidx[cbeg[c]];
               j = clen[c];
               while (j-- > 0)
                  rhs[*idx++] -= x * (*val++);
            }
         }
         break;
      }

      /*      Find nonzero with highest permuted row index and setup i and r
       */
      i = deQueueMax(ridx, &rn);
      r = rorig[i];

      x = diag[r] * rhs[r];
      rhs[r] = 0;
      if (x > eps || x < meps)
      {
         c = corig[i];
         work[c] = x;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
         {
            k = *idx++;
            y = rhs[k];
            if (y == 0)
            {
               y = -x * (*val++);
               if (y < meps || y > eps)
               {
                  rhs[k] = y;
                  enQueueMax(ridx, &rn, rperm[k]);
               }
            }
            else
            {
               y -= x * (*val++);
               y += ZERO * (y == 0);
               rhs[k] = y;
            }
         }
      }
   }
}

static
int vSolveUright2(CLUFactor* fac,
                   double* vec, int* vidx, double* rhs, int* ridx, int rn, double eps,
                   double* vec2, double* rhs2, int* ridx2, int rn2, double eps2)
{
   int i, j, k, r, c, n;
   int dim;
   int *rorig, *corig;
   int *rperm, *cperm;
   int *cidx, *clen, *cbeg;
   double *cval, *diag;
   double x, y, meps;
   double x2, y2, meps2;
   int lastRowSing, lastColSing;

   int *idx;
   double *val;
   double *work;
   double *work2;

   work = vec;
   work2 = vec2;
   rorig = fac->row.orig;
   corig = fac->col.orig;
   rperm = fac->row.perm;
   cperm = fac->col.perm;

   cidx = fac->u.col.idx;
   cval = fac->u.col.val;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;
   diag = fac->diag;
   dim = fac->thedim;

   lastRowSing = fac->u.lastRowSing;
   lastColSing = fac->u.lastColSing;

   meps = -eps;
   meps2 = -eps2;
   n = 0;

   while (rn + rn2 > 0)
   {
      /*      Find nonzero with highest permuted row index and setup i and r
       */
      if (rn <= 0)
         i = deQueueMax(ridx2, &rn2);
      else if (rn2 <= 0)
         i = deQueueMax(ridx, &rn);
      else if (*ridx2 > *ridx)
         i = deQueueMax(ridx2, &rn2);
      else if (*ridx2 < *ridx)
         i = deQueueMax(ridx, &rn);
      else
      {
         i = deQueueMax(ridx, &rn);
         i = deQueueMax(ridx2, &rn2);
      }
      r = rorig[i];

      x = diag[r] * rhs[r];
      x2 = diag[r] * rhs2[r];
      rhs[r] = 0;
      rhs2[r] = 0;
      if (x > eps || x < meps)
      {
         c = corig[i];
         vidx[n++] = c;
         work[c] = x;
         work2[c] = x2;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         if (x2 > eps2 || x2 < meps2)
         {
            while (j-- > 0)
            {
               k = *idx++;
               y2 = rhs2[k];
               if (y2 == 0)
               {
                  y2 = -x2 * (*val);
                  if (y2 < meps2 || y2 > eps2)
                  {
                     rhs2[k] = y2;
                     enQueueMax(ridx2, &rn2, rperm[k]);
                  }
               }
               else
               {
                  y2 -= x2 * (*val);
                  rhs2[k] = y2 + ZERO * (y2 == 0);
               }
               y = rhs[k];
               if (y == 0)
               {
                  y = -x * (*val++);
                  if (y < meps || y > eps)
                  {
                     rhs[k] = y;
                     enQueueMax(ridx, &rn, rperm[k]);
                  }
               }
               else
               {
                  y -= x * (*val++);
                  y += ZERO * (y == 0);
                  rhs[k] = y;
               }
            }
         }
         else
         {
            while (j-- > 0)
            {
               k = *idx++;
               y = rhs[k];
               if (y == 0)
               {
                  y = -x * (*val++);
                  if (y < meps || y > eps)
                  {
                     rhs[k] = y;
                     enQueueMax(ridx, &rn, rperm[k]);
                  }
               }
               else
               {
                  y -= x * (*val++);
                  y += ZERO * (y == 0);
                  rhs[k] = y;
               }
            }
         }
      }
      else if (x2 > eps2 || x2 < meps2)
      {
         c = corig[i];
         work2[c] = x2;
         val = &cval[cbeg[c]];
         idx = &cidx[cbeg[c]];
         j = clen[c];
         while (j-- > 0)
         {
            k = *idx++;
            y2 = rhs2[k];
            if (y2 == 0)
            {
               y2 = -x2 * (*val++);
               if (y2 < meps2 || y2 > eps2)
               {
                  rhs2[k] = y2;
                  enQueueMax(ridx2, &rn2, rperm[k]);
               }
            }
            else
            {
               y2 -= x2 * (*val++);
               rhs2[k] = y2 + ZERO * (y2 == 0);
            }
         }
      }

      if (rn + rn2 > i*verySparseFactor4right)
      {                                       /* continue with dense case */
         if (*ridx > *ridx2)
            i = *ridx;
         else
            i = *ridx2;
         for (; i >= 0; --i)
         {
            r = rorig[i];
            x = diag[r] * rhs[r];
            x2 = diag[r] * rhs2[r];
            rhs[r] = 0;
            rhs2[r] = 0;
            if (x2 > eps2 || x2 < meps2)
            {
               c = corig[i];
               work2[c] = x2;
               val = &cval[cbeg[c]];
               idx = &cidx[cbeg[c]];
               j = clen[c];
               if (x > eps || x < meps)
               {
                  vidx[n++] = c;
                  work[c] = x;
                  while (j-- > 0)
                  {
                     rhs[*idx] -= x * (*val);
                     rhs2[*idx++] -= x2 * (*val++);
                  }
               }
               else
               {
                  while (j-- > 0)
                     rhs2[*idx++] -= x2 * (*val++);
               }
            }
            else if (x > eps || x < meps)
            {
               c = corig[i];
               vidx[n++] = c;
               work[c] = x;
               val = &cval[cbeg[c]];
               idx = &cidx[cbeg[c]];
               j = clen[c];
               while (j-- > 0)
                  rhs[*idx++] -= x * (*val++);
            }
         }
         break;
      }
   }

   return n;
}


static
int vSolveUpdateRight(CLUFactor* fac, double* vec, int* ridx, int n, double eps)
{
   int i, j, k, l;
   int end;
   double meps, x, y;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   assert(!fac->l.updateType);               /* no Forest-Tomlin Updates */

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   meps = -eps;
   end = fac->l.firstUnused;

   for (i = fac->l.firstUpdate; i < end; ++i)
   {
      x = vec[lrow[i]];
      if (x > eps || x < meps)
      {
         k = lbeg[i];
         idx = &(lidx[k]);
         val = &(lval[k]);
         for (j = lbeg[i + 1]; j > k; --j)
         {
            ridx[n] = l = *idx++;
            y = vec[l];
            n += (y == 0);
            y = y - x * (*val++);
            vec[l] = y + (y == 0) * 1e-100;
         }
      }
   }

   return n;
}


static
void vSolveUpdateRightNoNZ(CLUFactor* fac, double* vec, double eps)
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


int vSolveRight4update(CLUFactor* fac, double eps,
                        double* vec, int* idx,                       /* result */
                        double* rhs, int* ridx, int rn,              /* rhs    */
                        double* forest, int* forestNum, int* forestIdx)
{
   rn = vSolveLright(fac, rhs, ridx, rn, eps);

   /*  turn index list into a heap
    */
   if (forest)
   {
      double meps, x;
      int i, j, k;
      int* rperm;
      int* it = forestIdx;

      rperm = fac->row.perm;
      meps = -eps;
      for (i = j = 0; i < rn; ++i)
      {
         k = ridx[i];
         x = rhs[k];
         if (x < meps || x > eps)
         {
            enQueueMax(ridx, &j, rperm[*it++ = k]);
            forest[k] = x;
         }
         else
            rhs[k] = 0;
      }
      *forestNum = rn = j;
   }
   else
   {
      double meps, x;
      int i, j, k;
      int* rperm;

      rperm = fac->row.perm;
      meps = -eps;
      for (i = j = 0; i < rn; ++i)
      {
         k = ridx[i];
         x = rhs[k];
         if (x < meps || x > eps)
            enQueueMax(ridx, &j, rperm[k]);
         else
            rhs[k] = 0;
      }
      rn = j;
   }

   rn = vSolveUright(fac, vec, idx, rhs, ridx, rn, eps);
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
      rn = vSolveUpdateRight(fac, vec, idx, rn, eps);
   return rn;
}

int vSolveRight4update2(CLUFactor* fac, double eps,
                         double* vec, int* idx,                  /* result1 */
                         double* rhs, int* ridx, int rn,         /* rhs1    */
                         double* vec2, double eps2,              /* result2 */
                         double* rhs2, int* ridx2, int rn2,      /* rhs2    */
                         double* forest, int* forestNum, int* forestIdx)
{
   /*
   rn  = vSolveLright(fac, rhs,  ridx,  rn,  eps);
   rn2 = vSolveLright(fac, rhs2, ridx2, rn2, eps2);
    */
   vSolveLright2(fac, rhs, ridx, &rn, eps, rhs2, ridx2, &rn2, eps2);

   /*  turn index list into a heap
    */
   if (forest)
   {
      double meps = -eps;
      double x;
      int i, j, k;
      int* rperm;
      int* it = forestIdx;

      rperm = fac->row.perm;
      for (i = j = 0; i < rn; ++i)
      {
         k = ridx[i];
         x = rhs[k];
         if (x < meps || x > eps)
         {
            enQueueMax(ridx, &j, rperm[*it++ = k]);
            forest[k] = x;
         }
         else
            rhs[k] = 0;
      }
      *forestNum = rn = j;
   }
   else
   {
      double meps = -eps;
      double x;
      int i, j, k;
      int* rperm;

      rperm = fac->row.perm;
      for (i = j = 0; i < rn; ++i)
      {
         k = ridx[i];
         x = rhs[k];
         if (x < meps || x > eps)
            enQueueMax(ridx, &j, rperm[k]);
         else
            rhs[k] = 0;
      }
      rn = j;
   }
   if (rn2 > fac->thedim*verySparseFactor4right)
   {
      ridx2[0] = fac->thedim - 1;
      /* ridx2[1] = fac->thedim - 2; */
   }
   else
   {
      double x;
      /*      double  maxabs; */
      int i, j, k;
      int* rperm;
      double meps = -eps2;

      /*      maxabs = 1;    */
      rperm = fac->row.perm;
      for (i = j = 0; i < rn2; ++i)
      {
         k = ridx2[i];
         x = rhs2[k];
         if (x < meps)
         {
            /*              maxabs = (maxabs < -x) ? -x : maxabs;  */
            enQueueMax(ridx2, &j, rperm[k]);
         }
         else if (x > eps2)
         {
            /*              maxabs = (maxabs < x) ? x : maxabs;    */
            enQueueMax(ridx2, &j, rperm[k]);
         }
         else
            rhs2[k] = 0;
      }
      rn2 = j;
      /*      eps2 = maxabs * eps2;  */
   }

   rn = vSolveUright(fac, vec, idx, rhs, ridx, rn, eps);
   vSolveUrightNoNZ(fac, vec2, rhs2, ridx2, rn2, eps2);

   /*
   rn = vSolveUright2(fac, vec, idx, rhs, ridx, rn, eps, vec2, rhs2, ridx2, rn2, eps2);
   */

   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      rn = vSolveUpdateRight(fac, vec, idx, rn, eps);
      vSolveUpdateRightNoNZ(fac, vec2, eps2);
   }

   return rn;
}

void vSolveRightNoNZ(CLUFactor* fac,
                      double* vec2, double eps2,              /* result2 */
                      double* rhs2, int* ridx2, int rn2)    /* rhs2    */
{
   rn2 = vSolveLright(fac, rhs2, ridx2, rn2, eps2);

   if (rn2 > fac->thedim*verySparseFactor4right)
   {
      *ridx2 = fac->thedim - 1;
   }
   else
   {
      double x;
      /*      double  maxabs; */
      int i, j, k;
      int* rperm;
      double meps = -eps2;

      /*      maxabs = 1;    */
      rperm = fac->row.perm;
      for (i = j = 0; i < rn2; ++i)
      {
         k = ridx2[i];
         x = rhs2[k];
         if (x < meps)
         {
            /*              maxabs = (maxabs < -x) ? -x : maxabs;  */
            enQueueMax(ridx2, &j, rperm[k]);
         }
         else if (x > eps2)
         {
            /*              maxabs = (maxabs < x) ? x : maxabs;    */
            enQueueMax(ridx2, &j, rperm[k]);
         }
         else
            rhs2[k] = 0;
      }
      rn2 = j;
      /*      eps2 = maxabs * eps2;  */
   }

   vSolveUrightNoNZ(fac, vec2, rhs2, ridx2, rn2, eps2);

   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
      vSolveUpdateRightNoNZ(fac, vec2, eps2);
}


/*****************************************************************************/
static
int solveUleft(CLUFactor* fac, double eps,
                double* vec, int* vecidx,
                double* rhs, int* rhsidx, int rhsn)
{
   double x, y, meps;
   int i, j, k, l, n, r, c;
   int *rorig, *corig, *cperm;
   int *ridx, *rlen, *rbeg, *idx;
   double *rval, *diag, *val;

   rorig = fac->row.orig;
   corig = fac->col.orig;
   cperm = fac->col.perm;

   /*  move rhsidx to a heap
    */
   for (i = 0; i < rhsn;)
      enQueueMin(rhsidx, &i, cperm[rhsidx[i]]);

   ridx = fac->u.row.idx;
   rval = fac->u.row.val;
   rlen = fac->u.row.len;
   rbeg = fac->u.row.start;

   diag = fac->diag;
   meps = -eps;
   n = 0;

   while (rhsn > 0)
   {
      i = deQueueMin(rhsidx, &rhsn);
      c = corig[i];
      x = rhs[c];
      rhs[c] = 0;
      if (x < meps || x > eps)
      {
         r = rorig[i];
         vecidx[n++] = r;
         x *= diag[r];
         vec[r] = x;
         k = rbeg[r];
         idx = &ridx[k];
         val = &rval[k];
         for (l = rlen[r]; l; --l)
         {
            j = *idx++;
            y = rhs[j];
            if (y == 0)
            {
               y = -x * (*val++);
               if (y < meps || y > eps)
               {
                  rhs[j] = y;
                  enQueueMin(rhsidx, &rhsn, cperm[j]);
               }
            }
            else
            {
               y -= x * (*val++);
               rhs[j] = y + ZERO * (y == 0);
            }
         }
      }
   }

   return n;
}

static
void solveUleftNoNZ(CLUFactor* fac, double eps, double* vec,
                     double* rhs, int* rhsidx, int rhsn)
{
   double x, y, meps;
   int i, j, k, l, r, c;
   int *rorig, *corig, *cperm;
   int *ridx, *rlen, *rbeg, *idx;
   double *rval, *diag, *val;

   rorig = fac->row.orig;
   corig = fac->col.orig;
   cperm = fac->col.perm;

   /*  move rhsidx to a heap
    */
   for (i = 0; i < rhsn;)
      enQueueMin(rhsidx, &i, cperm[rhsidx[i]]);

   ridx = fac->u.row.idx;
   rval = fac->u.row.val;
   rlen = fac->u.row.len;
   rbeg = fac->u.row.start;

   diag = fac->diag;
   meps = -eps;

   while (rhsn > 0)
   {
      i = deQueueMin(rhsidx, &rhsn);
      c = corig[i];
      x = rhs[c];
      rhs[c] = 0;
      if (x < meps || x > eps)
      {
         r = rorig[i];
         x *= diag[r];
         vec[r] = x;
         k = rbeg[r];
         idx = &ridx[k];
         val = &rval[k];
         for (l = rlen[r]; l; --l)
         {
            j = *idx++;
            y = rhs[j];
            if (y == 0)
            {
               y = -x * (*val++);
               if (y < meps || y > eps)
               {
                  rhs[j] = y;
                  enQueueMin(rhsidx, &rhsn, cperm[j]);
               }
            }
            else
            {
               y -= x * (*val++);
               rhs[j] = y + ZERO * (y == 0);
            }
         }
      }
   }
}

static
int solveLleftForest(CLUFactor* fac, double eps, double* vec, int* nonz, int n)
{
   int i, j, k, l, end;
   double x, y, meps;
   double *val, *lval;
   int *idx, *lidx, *lrow, *lbeg;

   meps = -eps;
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
         {
            l = *idx++;
            y = vec[l];
            if (y == 0)
            {
               y = -x * (*val++);
               if (y < meps || y > eps)
               {
                  vec[l] = y;
                  nonz[n++] = l;
               }
            }
            else
            {
               y -= x * (*val++);
               vec[l] = y + ZERO * (y == 0);
            }
         }
      }
   }

   return n;
}

static
void solveLleftForestNoNZ(CLUFactor* fac, double* vec)
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
}

static
int solveLleft(CLUFactor* fac, double eps, double* vec, int* nonz, int rn)
{
   int i, j, k, l, n;
   int r;
   double x, y, meps;
   double *rval, *lval, *val;
   int *ridx, *lidx, *idx, *lrow;
   int *rbeg, *lbeg;
   int *rorig, *rperm;
   int *last;

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   ridx = fac->l.ridx;
   rval = fac->l.rval;
   rbeg = fac->l.rbeg;
   rorig = fac->l.rorig;
   rperm = fac->l.rperm;
   meps = -eps;
   n = 0;

   i = fac->l.firstUpdate - 1;
#ifndef WITH_L_ROWS
   ERROR NOT YET IMPLEMENTED
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
   /*  move rhsidx to a heap
   */
   for (i = 0; i < rn;)
      enQueueMax(nonz, &i, rperm[nonz[i]]);
   last = nonz + fac->thedim;

   while (rn > 0)
   {
      i = deQueueMax(nonz, &rn);
      r = rorig[i];
      x = vec[r];
      if (x < meps || x > eps)
      {
         *(--last) = r;
         n++;
         k = rbeg[r];
         j = rbeg[r + 1] - k;
         val = &rval[k];
         idx = &ridx[k];
         while (j-- > 0)
         {
            assert(fac->l.rperm[*idx] < i);
            l = *idx++;
            y = vec[l];
            if (y == 0)
            {
               y = -x * *val++;
               if (y < meps || y > eps)
               {
                  vec[l] = y;
                  enQueueMax(nonz, &rn, rperm[l]);
               }
            }
            else
            {
               y -= x * *val++;
               vec[l] = y + ZERO * (y == 0);
            }
         }
      }
      else
         vec[r] = 0;
   }

   for (i = 0; i < n; ++i)
      *nonz++ = *last++;
#endif

   return n;
}

static
void solveLleftNoNZ(CLUFactor* fac, double* vec)
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

static
int solveUpdateLeft(CLUFactor* fac, double eps, double* vec, int* nonz, int n)
{
   int i, j, k, end;
   double x, y, meps;
   double *lval, *val;
   int *lrow, *lidx, *idx;
   int *lbeg;

   assert(!fac->l.updateType);               /* no Forest-Tomlin Updates! */

   lval = fac->l.val;
   lidx = fac->l.idx;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   meps = -eps;

   end = fac->l.firstUpdate;
   for (i = fac->l.firstUnused - 1; i >= end; --i)
   {
      k = lbeg[i];
      val = &lval[k];
      idx = &lidx[k];
      x = 0;
      for (j = lbeg[i + 1]; j > k; --j)
         x += vec[*idx++] * (*val++);
      k = lrow[i];
      y = vec[k];
      if (y == 0)
      {
         y = -x;
         if (y < meps || y > eps)
         {
            nonz[n++] = k;
            vec[k] = y;
         }
      }
      else
      {
         y -= x;
         vec[k] = y + ZERO * (y == 0);
      }
   }

   return n;
}

int vSolveLeft(CLUFactor* fac, double eps,
                double* vec, int* idx,                       /* result */
                double* rhs, int* ridx, int rn)            /* rhs    */
{
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      rn = solveUpdateLeft(fac, eps, rhs, ridx, rn);
      rn = solveUleft(fac, eps, vec, idx, rhs, ridx, rn);
   }
   else
   {
      rn = solveUleft(fac, eps, vec, idx, rhs, ridx, rn);
      rn = solveLleftForest(fac, eps, vec, idx, rn);
   }
   if (rn + fac->l.firstUpdate > verySparseFactor4left * fac->thedim)
   {
      solveLleftNoNZ(fac, vec);
      return 0;
   }
   else
      return solveLleft(fac, eps, vec, idx, rn);
}

int vSolveLeft2(CLUFactor* fac, double eps,
                 double* vec, int* idx,                      /* result */
                 double* rhs, int* ridx, int rn,             /* rhs    */
                 double* vec2,                               /* result2 */
                 double* rhs2, int* ridx2, int rn2)        /* rhs2    */
{
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      rn = solveUpdateLeft(fac, eps, rhs, ridx, rn);
      rn = solveUleft(fac, eps, vec, idx, rhs, ridx, rn);
      rn2 = solveUpdateLeft(fac, eps, rhs2, ridx2, rn2);
      solveUleftNoNZ (fac, eps, vec2, rhs2, ridx2, rn2);
   }
   else
   {
      rn = solveUleft(fac, eps, vec, idx, rhs, ridx, rn);
      rn = solveLleftForest(fac, eps, vec, idx, rn);
      solveUleftNoNZ(fac, eps, vec2, rhs2, ridx2, rn2);
      solveLleftForestNoNZ(fac, vec2);
   }
   rn = solveLleft(fac, eps, vec, idx, rn);
   solveLleftNoNZ (fac, vec2);

   return rn;
}

void vSolveLeftNoNZ(CLUFactor* fac, double eps,
                     double* vec2,                            /* result2 */
                     double* rhs2, int* ridx2, int rn2)     /* rhs2    */
{
   if (!fac->l.updateType)            /* no Forest-Tomlin Updates */
   {
      rn2 = solveUpdateLeft(fac, eps, rhs2, ridx2, rn2);
      solveUleftNoNZ (fac, eps, vec2, rhs2, ridx2, rn2);
   }
   else
   {
      solveUleftNoNZ(fac, eps, vec2, rhs2, ridx2, rn2);
      solveLleftForestNoNZ(fac, vec2);
   }
   solveLleftNoNZ (fac, vec2);
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
