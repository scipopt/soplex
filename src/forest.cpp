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
#pragma ident "@(#) $Id: forest.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "clutypes.h"
#include "clumembers.h"
#include "cluprotos.h"
#include "cring.h"

namespace soplex
{

/* TK20011102 
 * Nicht aendern, die Konstante wird benutzt um zwischen
 * explizit auf 0 gesetzt und ist 0 zu unterscheiden.
 * Sehr mehrkwuerdige Konstruktion.
 */
#define ZERO    1e-100

static const double verySparseFactor = 0.001;

/*****************************************************************************/
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
/*
 *      Perform garbage collection on column file
 */
static void packColumns(CLUFactor* fac)
{
   int n, i, j, col;
   Dring *ring, *list;

   double *cval = fac->u.col.val;
   int *cidx = fac->u.col.idx;
   int *clen = fac->u.col.len;
   int *cmax = fac->u.col.max;
   int *cbeg = fac->u.col.start;

   n = 0;
   list = &(fac->u.col.list);
   for (ring = list->next; ring != list; ring = ring->next)
   {
      col = ring->idx;
      if (cbeg[col] != n)
      {
         do
         {
            col = ring->idx;
            i = cbeg[col];
            cbeg[col] = n;
            cmax[col] = clen[col];
            j = i + clen[col];
            for (; i < j; ++i)
            {
               cval[n] = cval[i];
               cidx[n++] = cidx[i];
            }
            ring = ring->next;
         }
         while (ring != list);
         goto terminatePackColumns;
      }
      n += clen[col];
      cmax[col] = clen[col];
   }

terminatePackColumns :
   fac->u.col.used = n;
   fac->u.col.max[fac->thedim] = 0;
}

/*
 *      Ensure that column memory is at least size.
 */
static void minColMem(CLUFactor* fac, int size)
{
   if (fac->u.col.size < size)
   {
      fac->u.col.size = size;
      fac->u.col.idx = (int *)Realloc(fac->u.col.idx, size * sizeof(int));
      fac->u.col.val = (double*)Realloc(fac->u.col.val, size * sizeof(double));
      assert(fac->u.col.idx);
      assert(fac->u.col.val);
   }
}


/*
 *      Make column col of fac large enough to hold len nonzeros.
 */
static void reMaxCol(CLUFactor* fac, int col, int len)
{
   assert(fac->u.col.max[col] < len);

   if (fac->u.col.elem[col].next == &(fac->u.col.list)) /* last in column file */
   {
      int delta = len - fac->u.col.max[col];

      if (delta > fac->u.col.size - fac->u.col.used)
      {
         packColumns(fac);
         delta = len - fac->u.col.max[col];
         if (fac->u.col.size < fac->colMemMult * fac->u.col.used + len)
            minColMem(fac, int(fac->colMemMult * fac->u.col.used + len));
      }
      assert(delta <= fac->u.col.size - fac->u.col.used
             && "ERROR: could not allocate memory for column file");

      fac->u.col.used += delta;
      fac->u.col.max[col] = len;
   }

   else                        /* column must be moved to end of column file */
   {
      int i, j, k;
      int *idx;
      double *val;
      Dring *ring;

      if (len > fac->u.col.size - fac->u.col.used)
      {
         packColumns(fac);
         if (fac->u.col.size < fac->colMemMult * fac->u.col.used + len)
            minColMem(fac, int(fac->colMemMult * fac->u.col.used + len));
      }
      assert(len <= fac->u.col.size - fac->u.col.used
             && "ERROR: could not allocate memory for column file");

      j = fac->u.col.used;
      i = fac->u.col.start[col];
      k = fac->u.col.len[col] + i;
      fac->u.col.start[col] = j;
      fac->u.col.used += len;

      fac->u.col.max[fac->u.col.elem[col].prev->idx] += fac->u.col.max[col];
      fac->u.col.max[col] = len;
      removeDR(fac->u.col.elem[col]);
      ring = fac->u.col.list.prev;
      init2DR (fac->u.col.elem[col], *ring);

      idx = fac->u.col.idx;
      val = fac->u.col.val;
      for (; i < k; ++i)
      {
         val[j] = val[i];
         idx[j++] = idx[i];
      }
   }
}

/*****************************************************************************/


int forestUpdateCLUFactor(CLUFactor* fac, int col, double* work, int num, int *nonz)
{
   int i, j, k, l, m, n;
   int ll, c, r, row;
   double x;

   double *lval;
   int *lidx;
   int *lbeg = fac->l.start;

   double *cval = fac->u.col.val;
   int *cidx = fac->u.col.idx;
   int *cmax = fac->u.col.max;
   int *clen = fac->u.col.len;
   int *cbeg = fac->u.col.start;

   double *rval = fac->u.row.val;
   int *ridx = fac->u.row.idx;
   int *rmax = fac->u.row.max;
   int *rlen = fac->u.row.len;
   int *rbeg = fac->u.row.start;

   int *rperm = fac->row.perm;
   int *rorig = fac->row.orig;
   int *cperm = fac->col.perm;
   int *corig = fac->col.orig;

   double *diag = fac->diag;
   double maxabs = fac->maxabs;
   int dim = fac->thedim;
   int stat = CLU_OK;


   /*  Remove column col form U
    */
   j = cbeg[col];
   i = clen[col];
   fac->nzCnt -= i;
   for (i += j - 1; i >= j; --i)
   {
      m = cidx[i];
      k = rbeg[m];
      l = --(rlen[m]) + k;
      while (ridx[k] != col)
         ++k;
      ridx[k] = ridx[l];
      rval[k] = rval[l];
   }


   /*  Insert new vector column col
    *  thereby determining the highest permuted row index r.
    */
   if (num)
   {
      clen[col] = 0;
      if (num > cmax[col])
         reMaxCol(fac, col, num);
      cidx = fac->u.col.idx;
      cval = fac->u.col.val;
      k = cbeg[col];
      r = 0;
      for (j = 0; j < num; ++j)
      {
         i = *nonz++;
         x = work[i];
         work[i] = 0;
         if (x > 1.e-12 || x < -1.e-12)
         {
            if (x > maxabs)
               maxabs = x;
            else if (-x > maxabs)
               maxabs = -x;

            /* insert to column file */
            assert(k - cbeg[col] < cmax[col]);
            cval[k] = x;
            cidx[k++] = i;

            /* insert to row file */
            if (rmax[i] <= rlen[i])
            {
               remaxRow(fac, i, rlen[i] + 1);
               rval = fac->u.row.val;
               ridx = fac->u.row.idx;
            }
            l = rbeg[i] + (rlen[i])++;
            rval[l] = x;
            ridx[l] = col;

            /* check permuted row index */
            if (rperm[i] > r)
               r = rperm[i];
         }
      }
      fac->nzCnt += (clen[col] = k - cbeg[col]);
   }

   else
   {
      /*
      clen[col] = 0;
      reMaxCol(fac, col, dim);
       */
      cidx = fac->u.col.idx;
      cval = fac->u.col.val;
      k = cbeg[col];
      j = k + cmax[col];
      r = 0;
      for (i = 0; i < dim; ++i)
      {
         x = work[i];
         work[i] = 0;
         if (x > 1.e-12 || x < -1.e-12)
         {
            if (x > maxabs)
               maxabs = x;
            else if (-x > maxabs)
               maxabs = -x;

            /* insert to column file */
            if (k >= j)
            {
               clen[col] = k - cbeg[col];
               reMaxCol(fac, col, dim - i);
               cidx = fac->u.col.idx;
               cval = fac->u.col.val;
               k = cbeg[col];
               j = k + cmax[col];
               k += clen[col];
            }
            assert(k - cbeg[col] < cmax[col]);
            cval[k] = x;
            cidx[k++] = i;

            /* insert to row file */
            if (rmax[i] <= rlen[i])
            {
               remaxRow(fac, i, rlen[i] + 1);
               rval = fac->u.row.val;
               ridx = fac->u.row.idx;
            }
            l = rbeg[i] + (rlen[i])++;
            rval[l] = x;
            ridx[l] = col;

            /* check permuted row index */
            if (rperm[i] > r)
               r = rperm[i];
         }
      }
      fac->nzCnt += (clen[col] = k - cbeg[col]);
      if (cbeg[col] + cmax[col] == fac->u.col.used)
      {
         fac->u.col.used -= cmax[col];
         cmax[col] = clen[col];
         fac->u.col.used += cmax[col];
      }
   }


   /*  Adjust stages of column and row singletons in U.
    */
   fac->u.lastRowSing = fac->u.lastColSing;

   c = cperm[col];
   if (r > c)                         /* Forest Tomlin update */
   {
      /*      update permutations
       */
      j = rorig[c];
      for (i = c; i < r; ++i)
         rorig[i] = rorig[i + 1];
      rorig[r] = j;
      for (i = c; i <= r; ++i)
         rperm[rorig[i]] = i;

      j = corig[c];
      for (i = c; i < r; ++i)
         corig[i] = corig[i + 1];
      corig[r] = j;
      for (i = c; i <= r; ++i)
         cperm[corig[i]] = i;


      row = rorig[r];
      j = rbeg[row];
      i = rlen[row];
      fac->nzCnt -= i;

      if (i < verySparseFactor*(dim - c))      /* few nonzeros to be eliminated        */
      {
         /*  move row r from U to work
         fprintf(stderr, ".");
          */
         num = 0;
         for (i += j - 1; i >= j; --i)
         {
            k = ridx[i];
            work[k] = rval[i];
            enQueueMin(nonz, &num, cperm[k]);
            m = --(clen[k]) + cbeg[k];
            for (l = m; cidx[l] != row; --l)
              ;
            cidx[l] = cidx[m];
            cval[l] = cval[m];
         }


         /*  Eliminate row r from U to L file
          */
         ll = makeLvec(fac, r - c, row);
         lval = fac->l.val;
         lidx = fac->l.idx;
         /* for(i = c; i < r; ++i)       */
         while (num)
         {
#ifndef NDEBUG
            for (i = 0; i < num; ++i)
            {
               if (work[corig[nonz[i]]] == 0)
                  assert(0);
            }
#endif  // NDEBUG
            i = deQueueMin(nonz, &num);
            if (i == r)
               break;
            k = corig[i];
            assert(work[k]);
            {
               n = rorig[i];
               x = work[k] * diag[n];
               lidx[ll] = n;
               lval[ll] = x;
               work[k] = 0;
               ll++;

               if (x > maxabs)
                  maxabs = x;
               else if (-x > maxabs)
                  maxabs = -x;

               j = rbeg[n];
               m = rlen[n] + j;
               for (; j < m; ++j)
               {
                  int jj = ridx[j];
                  double y = work[jj];
                  if (y == 0)
                     enQueueMin(nonz, &num, cperm[jj]);
                  y -= x * rval[j];
                  work[jj] = y + (y == 0) * ZERO;
               }
            }
         }
         if (lbeg[fac->l.firstUnused - 1] == ll)
            (fac->l.firstUnused)--;
         else
            lbeg[fac->l.firstUnused] = ll;


         /*  Set diagonal value
          */
         if (i != r)
            return fac->stat = CLU_SINGULAR;
         k = corig[r];
         x = work[k];
         diag[row] = 1 / x;
         work[k] = 0;


         /*  make row large enough to fit all nonzeros.
          */
         if (rmax[row] < num)
         {
            rlen[row] = 0;
            remaxRow(fac, row, num);
            rval = fac->u.row.val;
            ridx = fac->u.row.idx;
         }
         fac->nzCnt += num;

         /*  Insert work to updated row thereby clearing work;
          */
         n = rbeg[row];
         for (i = 0; i < num; ++i)
         {
            j = corig[nonz[i]];
            x = work[j];
            assert(x != 0.0);
            {
               if (x > maxabs)
                  maxabs = x;
               else if (-x > maxabs)
                  maxabs = -x;

               ridx[n] = j;
               rval[n] = x;
               work[j] = 0;
               ++n;

               if (clen[j] >= cmax[j])
               {
                  reMaxCol(fac, j, clen[j] + 1);
                  cidx = fac->u.col.idx;
                  cval = fac->u.col.val;
               }
               cval[cbeg[j] + clen[j]] = x;
               cidx[cbeg[j] + clen[j]++] = row;
            }
         }
         rlen[row] = n - rbeg[row];
      }

      else            /* few nonzeros to be eliminated        */
      {
         /*  move row r from U to work
         fprintf(stderr, " ");
          */
         for (i += j - 1; i >= j; --i)
         {
            k = ridx[i];
            work[k] = rval[i];
            m = --(clen[k]) + cbeg[k];
            for (l = m; cidx[l] != row; --l)
              ;
            cidx[l] = cidx[m];
            cval[l] = cval[m];
         }


         /*  Eliminate row r from U to L file
          */
         ll = makeLvec(fac, r - c, row);
         lval = fac->l.val;
         lidx = fac->l.idx;
         for (i = c; i < r; ++i)
         {
            k = corig[i];
            if (work[k])
            {
               n = rorig[i];
               x = work[k] * diag[n];
               lidx[ll] = n;
               lval[ll] = x;
               work[k] = 0;
               ll++;

               if (x > maxabs)
                  maxabs = x;
               else if (-x > maxabs)
                  maxabs = -x;

               j = rbeg[n];
               m = rlen[n] + j;
               for (; j < m; ++j)
                  work[ridx[j]] -= x * rval[j];
            }
         }
         if (lbeg[fac->l.firstUnused - 1] == ll)
            (fac->l.firstUnused)--;
         else
            lbeg[fac->l.firstUnused] = ll;


         /*  Set diagonal value
          */
         k = corig[r];
         x = work[k];
         if (x == 0.0)
            return fac->stat = CLU_SINGULAR;
         diag[row] = 1 / x;
         work[k] = 0;


         /*  count remaining nonzeros in work and make row large enough
          *  to fit them all.
          */
         n = 0;
         for (i = r + 1; i < dim; ++i)
            n += (work[corig[i]] != 0.0);
         if (rmax[row] < n)
         {
            rlen[row] = 0;
            remaxRow(fac, row, n);
            rval = fac->u.row.val;
            ridx = fac->u.row.idx;
         }
         fac->nzCnt += n;

         /*  Insert work to updated row thereby clearing work;
          */
         n = rbeg[row];
         for (i = r + 1; i < dim; ++i)
         {
            j = corig[i];
            x = work[j];
            if (x != 0.0)
            {
               if (x > maxabs)
                  maxabs = x;
               else if (-x > maxabs)
                  maxabs = -x;

               ridx[n] = j;
               rval[n] = x;
               work[j] = 0;
               ++n;

               if (clen[j] >= cmax[j])
               {
                  reMaxCol(fac, j, clen[j] + 1);
                  cidx = fac->u.col.idx;
                  cval = fac->u.col.val;
               }
               cval[cbeg[j] + clen[j]] = x;
               cidx[cbeg[j] + clen[j]++] = row;
            }
         }
         rlen[row] = n - rbeg[row];
      }
   }

   else if (r == c)
   {
      /*  Move diagonal element to diag.  Note, that it must be the last
       *  element, since it has just been inserted above.
       */
      row = rorig[r];
      i = rbeg[row] + --(rlen[row]);
      diag[row] = 1 / rval[i];
      for (j = i = --(clen[col]) + cbeg[col]; cidx[i] != row; --i)
        ;
      cidx[i] = cidx[j];
      cval[i] = cval[j];
   }

   else /* r < c */
      return fac->stat = CLU_SINGULAR;

   fac->maxabs = maxabs;

   assert(CLUFactorIsConsistent(fac));
   return stat;
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
