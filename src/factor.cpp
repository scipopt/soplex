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
#pragma ident "@(#) $Id: factor.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"


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
/*
 *      Global data.
 */
static CLUFactor* fac;            /* currently factorized matrix */
static double epsilon;        /* epsilon to use */
static double threshold;      /* epsilon to use */
static double initMaxabs;     /* maximum abs value in initial matrix */
static int stage;
static int dim;

/*      Woring matrix U
 */
static int *cidx;
static int *cmax;
static int *clen;
static int *cbeg;

static double *rval;
static int *ridx;
static int *rmax;
static int *rlen;
static int *rbeg;

static double *diag;
static double *work;

/*      matrix L
 */
static double *lval;
static int *lidx;
static int *lbeg;
static int *lrow;


/*****************************************************************************/
/*
 *      Temporary data structures.
 */
static int*    s_mark;
static double* s_max;           /* maximum absolute value per row (or -1) */
static int*    s_cact;          /* lengths of columns of active submatrix */
/*
        For the i=th column the situation might look like this:
 
\begin{verbatim}
        idx     = ....................iiiIIIIII-----..............
        cbeg[i] =                     ^
        cact[i] =                        +----+
        clen[i] =                     +-------+
        cmax[i] =                     +------------+
 
        Indices clen[i]-cbeg[i]:      ^^^
\end{verbatim}
        belong to column i, but have allready been pivotal and don't belong to
        the active submatrix.
 */

typedef struct _pring   /* Pivot ring */
{
   struct _pring *next;
   struct _pring *prev;
   int idx;            /* index of pivot row */
   int pos;            /* position of pivot column in row */
   int mkwtz;          /* markowitz number of pivot */
}
Pring;

static Pring pivots;         /* ring of selected pivot rows */
static Pring *col,            /* column index handlers for double linked list */
*colNZ,          /* lists for columns to number of nonzeros      */
*row,            /* same for rows */
*rowNZ;         /* same for rows */


static void newTmp(int dim)
{
   s_max = (double*)Malloc(dim * sizeof(double));
   assert(s_max);

   s_cact = (int *)Malloc(dim * sizeof(int));
   assert(s_cact);

   s_mark = (int *)Malloc(dim * sizeof(int));
   assert(s_mark);
}

static void deleteTmp()
{
   Free(s_max);
   s_max = 0;

   Free(s_mark);
   s_mark = 0;

   Free(s_cact);
   s_cact = 0;
}

/*****************************************************************************/
/*
 *      Ensure that row memory is at least size.
 */
static void minRowMem(int size)
{
   if (fac->u.row.size < size)
   {
      fac->u.row.size = size;
      fac->u.row.val = (double*)Realloc(fac->u.row.val, size * sizeof(double));
      fac->u.row.idx = (int *)Realloc(fac->u.row.idx, size * sizeof(int));
      assert(fac->u.row.idx);
      assert(fac->u.row.val);
      ridx = fac->u.row.idx;
      rval = fac->u.row.val;
   }
}

/*****************************************************************************/
/*
 *      Ensure that column memory is at least size.
 */
static void minColMem(int size)
{
   if (fac->u.col.size < size)
   {
      fac->u.col.size = size;
      fac->u.col.idx = (int *)Realloc(fac->u.col.idx, size * sizeof(int));
      assert(fac->u.col.idx);
      cidx = fac->u.col.idx;
   }
}


/*****************************************************************************/
/*
 *      Make new Lvector to fit len nonzeros.
 *      return index to the 
 */
static void minLMem(CLUFactor *fac, int size)
{
   if (size > fac->l.size)
   {
      fac->l.size = int(0.2 * fac->l.size + size);
      fac->l.val = (double*)Realloc(fac->l.val, fac->l.size * sizeof(double));
      fac->l.idx = (int *)Realloc(fac->l.idx, fac->l.size * sizeof(int));
      assert(fac->l.idx);
      assert(fac->l.val);
      lidx = fac->l.idx;
      lval = fac->l.val;
   }
}

int makeLvec(CLUFactor* fac, int len, int row)
{
   int* lrow = fac->l.row;
   int* lbeg = fac->l.start;
   int first = lbeg[fac->l.firstUnused];

   if (fac->l.firstUnused >= fac->l.startSize)
   {
      fac->l.startSize += 100;
      fac->l.start = (int *)Realloc(fac->l.start, fac->l.startSize * sizeof(int));
      lbeg = fac->l.start;
   }

   assert(len > 0 && "ERROR: no empty columns allowed in L vectors");

   minLMem(fac, first + len);
   lrow[fac->l.firstUnused] = row;
   lbeg[++(fac->l.firstUnused)] = first + len;

   assert(lbeg[fac->l.firstUnused] <= fac->l.size);
   assert(fac->l.firstUnused <= fac->l.startSize);
   return first;
}

/*****************************************************************************/

static void initPerm()
{
   int* p_or = fac->row.orig;
   int* p_pr = fac->row.perm;
   int* p_oc = fac->col.orig;
   int* p_pc = fac->col.perm;
   int  i;

   for (i = fac->thedim - 1; i >= 0; --i)
      p_or[i] = p_pr[i] = p_oc[i] = p_pc[i] = -1;
}

static void setPivot(int stage, int col, int row, double val)
{
   /*@
   assert(printf("%5d:   (%3d, %3d) = %g\n", stage, row, col, val) || 1);
   */
   assert(fac->row.perm[row] < 0);
   assert(fac->col.perm[col] < 0);
   fac->row.orig[stage] = row;
   fac->col.orig[stage] = col;
   fac->row.perm[row] = stage;
   fac->col.perm[col] = stage;
   diag[row] = 1 / val;
   if (diag[row] > fac->maxabs)
      fac->maxabs = diag[row];
   else if (-diag[row] > fac->maxabs)
      fac->maxabs = -diag[row];
}


/*****************************************************************************/
/*
 *      Initialize row and column file of working matrix and
 *      mark column singletons.
 */
static void initMatrix(SVector** vec)
{
   double x;
   int i, j, l, k, m;
   int tot;
   Dring *rring, *lastrring;
   Dring *cring, *lastcring;
   SVector *psv;
   int *sing = s_mark;

   /*  Initialize:
    *  - column file thereby remembering column singletons in |sing|.
    *  - nonzeros counts per row
    *  - total number of nonzeros
    */
   for (i = 0; i < dim; ++i)
      rmax[i] = rlen[i] = 0;

   tot = 0;

   for (i = 0; i < dim; ++i)
   {
      psv = vec[i];
      l = psv->size();
      if (l > 1)
      {
         tot += l;
         for (j = 0; j < l; ++j)
            ++(rmax[psv->index(j)]);
      }
      else if (l == 0)
      {
         fac->stat = CLU_SINGULAR;
         return;
      }
   }



   /*  Resize nonzero memory if necessary
    */
   minRowMem(int(fac->rowMemMult * tot));
   minColMem(int(fac->colMemMult * tot));
   minLMem(fac, int(fac->lMemMult * tot));


   /*  Initialize:
    *  - row ring lists
    *  - row vectors in file
    *  - column ring lists
    */
   rbeg[0] = 0;

   rring = fac->u.row.elem;
   lastrring = &(fac->u.row.list);
   lastrring->idx = dim;
   lastrring->next = rring;

   cring = fac->u.col.elem;
   lastcring = &(fac->u.col.list);
   lastcring->idx = dim;
   lastcring->next = cring;

   m = 0;
   for (i = 0; i < dim; ++i)
   {
      rbeg[i] = m;
      m += rmax[i];

      rring->idx = i;
      rring->prev = lastrring;
      lastrring->next = rring;
      lastrring = rring;
      ++rring;

      cring->idx = i;
      cring->prev = lastcring;
      lastcring->next = cring;
      lastcring = cring;
      ++cring;
   }
   fac->u.row.used = m;

   lastrring->next = &(fac->u.row.list);
   lastrring->next->prev = lastrring;

   lastcring->next = &(fac->u.col.list);
   lastcring->next->prev = lastcring;



   /*  Copy matrix to row and column file
    *  excluding and marking column singletons!
    */
   m = 0;
   stage = 0;

   initMaxabs = 0;
   for (i = 0; i < dim; ++i)
   {
      psv = vec[i];
      l = psv->size();
      cbeg[i] = m;
      if (l > 1)                               /* exclude column singletons */
      {
         int kk, ll;
         for (j = ll = 0; j < l; ++j)
         {
            x = psv->value(j);
            if (isNonZero(x, epsilon))
            {
               k = psv->index(j);
               kk = rbeg[k] + (rlen[k]++);
               cidx[m++] = k;
               ridx[kk] = i;
               rval[kk] = x;
               ++ll;
               if (x > initMaxabs)
                  initMaxabs = x;
               else if (-x > initMaxabs)
                  initMaxabs = -x;
            }
         }
         l = ll;
         --m;
      }
      if (l > 1)
      {
         ++m;
         s_cact[i] = clen[i] = cmax[i] = l;
      }
      else if (l <= 0)       /* singular */
      {
         fac->stat = CLU_SINGULAR;
         return;
      }
      else                    /* singleton */
      {
         clen[i] = cmax[i] = 0;
         for (j = 0; isZero(psv->value(j), epsilon); ++j)
           ;
         if (fac->row.perm[psv->index(j)] >= 0)
         {
            fac->stat = CLU_SINGULAR;
            return;
         }
         x = psv->value(j);
         if (x > initMaxabs)
            initMaxabs = x;
         else if (-x > initMaxabs)
            initMaxabs = -x;
         setPivot(stage, i, psv->index(j), x);
         sing[stage++] = i;
      }
   }

   fac->u.col.used = m;
}



/*****************************************************************************/
/*
 *      Remove column singletons
 */

static void colSingletons(void)
{
   int i, j, k, n;
   int len;
   int col, row, newrow;
   int *idx;
   int *rorig = fac->row.orig;
   int *rperm = fac->row.perm;
   int *sing = s_mark;


   /*  Iteratively update column counts due to removed column singletons
    *  thereby removing new arising columns singletons
    *  and computing the index of the first row singleton (-1)
    *  until no more can be found.
    */
   fac->u.lastColSing = -1;
   for (i = 0; i < stage; ++i)
   {
      row = rorig[i];
      assert(row >= 0);
      idx = &(ridx[rbeg[row]]);
      len = rlen[row];

      if (len)
         fac->u.lastColSing = i;

      for (j = 0; j < len; ++j)
      {
         /*  Move pivotal nonzeros to front of column.
          */
         col = idx[j];
         n = cbeg[col] + clen[col] - s_cact[col];
         for (k = n; cidx[k] != row; ++k)
           ;
         assert(k < cbeg[col] + clen[col]);
         cidx[k] = cidx[n];
         cidx[n] = row;

         n = --(s_cact[col]);          /* column nonzeros of ACTIVE matrix */

         if (n == 1)                  /* Here is another singleton */
         {
            newrow = cidx[--clen[col] + cbeg[col]];

            /*      Ensure, matrix not singular
             */
            if (rperm[newrow] >= 0)
            {
               fac->stat = CLU_SINGULAR;
               return;
            }

            /*      Find singleton in row.
             */
            n = rbeg[newrow] + (--(rlen[newrow]));
            for (k = n; ridx[k] != col; --k)
              ;

            /*      Remove singleton from column.
             */
            setPivot(stage, col, newrow, rval[k]);
            sing[stage++] = col;

            /*      Move pivot element to diag.
             */
            rval[k] = rval[n];
            ridx[k] = ridx[n];
         }
         else if (n == 0)
         {
            fac->stat = CLU_SINGULAR;
            return;
         }
      }
   }

   assert(stage <= fac->thedim);
}


/*****************************************************************************/
/*
 *      Remove row singletons
 */
static void rowSingletons(void)
{
   double pval;
   int i, j, k, l, r;
   int row, col, len, rs, lk;
   int *idx;
   int *rperm = fac->row.perm;
   int *sing = s_mark;

   /*  Mark row singletons
    */
   rs = stage;
   for (i = 0; i < dim; ++i)
   {
      if (rperm[i] < 0 && rlen[i] == 1)
         sing[stage++] = i;
   }

   /*  Eliminate row singletons
    *  thereby marking newly arising ones
    *  until no more can be found.
    */
   for (; rs < stage; ++rs)
   {
      /*      Move pivot element from row file to diag
       */
      row = sing[rs];
      j = rbeg[row];
      col = ridx[j];
      pval = rval[j];
      setPivot(rs, col, row, pval);
      rlen[row] = 0;

      /*      Remove pivot column form workingmatrix
       *      thereby building up L vector.
       */
      idx = &(cidx[cbeg[col]]);
      i = s_cact[col];                /* nr. nonzeros of new L vector */
      lk = makeLvec(fac, i - 1, row);
      len = clen[col];
      i = (clen[col] -= i);         /* remove pivot column from U */

      for (; i < len; ++i)
      {
         r = idx[i];
         if (r != row)
         {
            /*      Find pivot column in row.
             */
            l = --(rlen[r]);
            k = rbeg[r] + l;
            for (j = k; ridx[j] != col; --j)
              ;
            assert(k >= rbeg[r]);

            /*      Initialize L vector
             */
            lidx[lk] = r;
            lval[lk] = rval[j] / pval;
            ++lk;

            /*      Remove pivot column from row.
             */
            ridx[j] = ridx[k];
            rval[j] = rval[k];

            /*      Check new row length.
             */
            if (l == 1)
               sing[stage++] = r;
            else if (l == 0)
            {
               fac->stat = CLU_SINGULAR;
               return;
            }
         }
      }
   }

   fac->u.lastRowSing = stage - 1;
}


/*****************************************************************************/
/*
 *      Init nonzero number Ring lists
 *      and required entries of arrays max and mark
 */

static void initRings(void)
{
   int i;
   int *rperm = fac->row.perm;
   int *cperm = fac->col.perm;
   Pring *ring;

   col = (Pring*)Malloc((dim + 1) * sizeof(Pring));
   colNZ = (Pring*)Malloc((dim + 1) * sizeof(Pring));
   row = (Pring*)Malloc((dim + 1) * sizeof(Pring));
   rowNZ = (Pring*)Malloc((dim + 1) * sizeof(Pring));
   assert(col && colNZ && row && rowNZ);

   for (i = dim - stage; i >= 0; --i)
   {
      initDR(colNZ[i]);
      initDR(rowNZ[i]);
   }

   for (i = 0; i < dim; ++i)
   {
      if (rperm[i] < 0)
      {
         assert(rlen[i] > 1);
         ring = &(rowNZ[rlen[i]]);
         init2DR(row[i], *ring);
         row[i].idx = i;
         s_max[i] = -1;
      }
      if (cperm[i] < 0)
      {
         assert(s_cact[i] > 1);
         ring = &(colNZ[s_cact[i]]);
         init2DR(col[i], *ring);
         col[i].idx = i;
         s_mark[i] = 0;
      }
   }
}

static void freeRings(void)
{
   Free(col);
   Free(colNZ);
   Free(row);
   Free(rowNZ);
   col = colNZ = row = rowNZ = 0;
}


/*****************************************************************************/

/*
 *      Eliminate all row singletons from nucleus.
 *      A row singleton may well be column singleton at the same time!
 */
static void eliminateRowSingletons(void)
{
   int i, j, k, l, r;
   int len, lk;
   int pcol, prow;
   double pval;
   int *idx;
   Pring *sing;

   for (sing = rowNZ[1].prev; sing != &(rowNZ[1]); sing = sing->prev)
   {
      prow = sing->idx;
      i = rbeg[prow];
      pcol = ridx[i];
      pval = rval[i];
      setPivot(stage++, pcol, prow, pval);
      rlen[prow] = 0;
      removeDR(col[pcol]);

      /*      Eliminate pivot column and build L vector.
       */
      i = s_cact[pcol];
      if (i > 1)
      {
         idx = &(cidx[cbeg[pcol]]);
         len = clen[pcol];
         lk = makeLvec(fac, i - 1, prow);
         i = clen[pcol] -= i;

         for (; (r = idx[i]) != prow; ++i)
         {
            /*      Find pivot column in row.
             */
            l = --(rlen[r]);
            k = rbeg[r] + l;
            for (j = k; ridx[j] != pcol; --j)
              ;
            assert(j >= rbeg[r]);

            /*      Initialize L vector
             */
            lidx[lk] = r;
            lval[lk] = rval[j] / pval;
            ++lk;

            /*      Remove pivot column from row.
             */
            ridx[j] = ridx[k];
            rval[j] = rval[k];

            /*      Move column to appropriate nonzero ring.
             */
            removeDR(row[r]);
            init2DR (row[r], rowNZ[l]);
            assert(fac->row.perm[r] < 0);
            s_max[r] = -1;
         }

         /* skip pivot element */
         assert(i < len && "ERROR: pivot column does not contain pivot row");

         for (++i; i < len; ++i)
         {
            /*      Find pivot column in row.
             */
            r = idx[i];
            l = --(rlen[r]);
            k = rbeg[r] + l;
            for (j = k; ridx[j] != pcol; --j)
              ;
            assert(j >= rbeg[r]);

            /*      Initialize L vector
             */
            lidx[lk] = r;
            lval[lk] = rval[j] / pval;
            ++lk;

            /*      Remove pivot column from row.
             */
            ridx[j] = ridx[k];
            rval[j] = rval[k];

            /*      Move column to appropriate nonzero ring.
             */
            removeDR(row[r]);
            init2DR (row[r], rowNZ[l]);
            assert(fac->row.perm[r] < 0);
            s_max[r] = -1;
         }
      }
      else
         clen[pcol] -= i;
   }

   initDR(rowNZ[1]);           /* Remove all row singletons from list */
}



/*
 *      Eliminate all column singletons from nucleus.
 *      A column singleton must not be row singleton at the same time!
 */
static void eliminateColSingletons(void)
{
   int i, j, k, l, c;
   int pcol, prow;
   Pring *sing;

   for (sing = colNZ[1].prev; sing != &(colNZ[1]); sing = sing->prev)
   {
      /*      Find pivot value
       */
      pcol = sing->idx;
      j = --(clen[pcol]) + cbeg[pcol];     /* remove pivot column */
      prow = cidx[j];
      removeDR(row[prow]);

      j = --(rlen[prow]) + rbeg[prow];
      for (i = j; (c = ridx[i]) != pcol; --i)
      {
         l = clen[c] + cbeg[c] - (s_cact[c])--;
         for (k = l; cidx[k] != prow; ++k)
           ;
         cidx[k] = cidx[l];
         cidx[l] = prow;
         l = s_cact[c];
         removeDR(col[c]);
         init2DR(col[c], colNZ[l]);
         assert(fac->col.perm[c] < 0);
      }

      /*      remove pivot element from pivot row
       */
      setPivot(stage++, pcol, prow, rval[i]);
      ridx[i] = ridx[j];
      rval[i] = rval[j];

      j = rbeg[prow];
      for (--i; i >= j; --i)
      {
         c = ridx[i];
         l = clen[c] + cbeg[c] - (s_cact[c])--;
         for (k = l; cidx[k] != prow; ++k)
           ;
         cidx[k] = cidx[l];
         cidx[l] = prow;
         l = s_cact[c];
         removeDR(col[c]);
         init2DR(col[c], colNZ[l]);
         assert(fac->col.perm[c] < 0);
      }
   }

   initDR(colNZ[1]);           /* Remove all column singletons from list */
}


/*
 *      No singletons available:        Select pivot elements.
 */
static void selectPivots(void)
{
   int ii, i, j, k, ll, l, m;
   int count, num, rw, cl;
   int len, beg;
   double maxabs, x;
   long mkwtz;
   int candidates;

   candidates = dim - stage - 1;
   if (candidates > 4)
      candidates = 4;

   num = 0;
   count = 2;

   for (;;)
   {
      ii = -1;

      if (rowNZ[count].next != &(rowNZ[count]))
      {
         rw = rowNZ[count].next->idx;
         beg = rbeg[rw];
         len = rlen[rw] + beg - 1;

         /*  set maxabs to maximum absolute value in row
          *  (compute it if necessary).
          */
         if ((maxabs = s_max[rw]) < 0)
         {
            maxabs = rval[len];
            if (maxabs < 0)
               maxabs = -maxabs;
            for (i = len - 1; i >= beg; --i)
               if (maxabs < rval[i])
                  maxabs = rval[i];
               else if (maxabs < -rval[i])
                  maxabs = -rval[i];
            s_max[rw] = maxabs;               /* ##### */
         }
         maxabs *= threshold;

         /*  select pivot element with lowest markowitz number in row
          */
         mkwtz = dim + 1;
         for (i = len; i >= beg; --i)
         {
            k = ridx[i];
            j = s_cact[k];
            x = rval[i];
            if (j < mkwtz && (x > maxabs || -x > maxabs))
            {
               mkwtz = j;
               cl = k;
               ii = i;
               if (j <= count)              /* ##### */
                  break;
            }
         }
      }

      else if (colNZ[count].next != &(colNZ[count]))
      {
         cl = colNZ[count].next->idx;
         beg = cbeg[cl];
         len = clen[cl] + beg - 1;
         beg = len - s_cact[cl] + 1;
         assert(count == s_cact[cl]);

         /*  select pivot element with lowest markowitz number in column
          */
         mkwtz = dim + 1;
         for (i = len; i >= beg; --i)
         {
            k = cidx[i];
            j = rlen[k];
            if (j < mkwtz)
            {
               /*  ensure that element (cl,k) is stable.
                */
               if (s_max[k] > 0)
               {
                  /*  case 1: maxabs is known
                   */
                  for (m = rbeg[k], l = m + rlen[k] - 1; l >= m; --l)
                  {
                     if (ridx[l] == cl)
                     {
                        x = rval[l];
                        ll = l;
                        break;
                     }
                  }
                  maxabs = s_max[k];
               }
               else
               {
                  /*  case 2: maxabs needs to be computed
                   */
                  m = rbeg[k];
                  maxabs = rval[m];
                  maxabs = (maxabs > 0) ? maxabs : -maxabs;
                  for (l = m + rlen[k] - 1; l >= m; --l)
                  {
                     if (maxabs < rval[l])
                        maxabs = rval[l];
                     else if (maxabs < -rval[l])
                        maxabs = -rval[l];
                     if (ridx[l] == cl)
                     {
                        x = rval[l];
                        ll = l;
                        break;
                     }
                  }
                  for (--l; l > m; --l)
                  {
                     if (maxabs < rval[l])
                        maxabs = rval[l];
                     else if (maxabs < -rval[l])
                        maxabs = -rval[l];
                  }
                  s_max[k] = maxabs;
               }
               maxabs *= threshold;

               if (x > maxabs || -x > maxabs)
               {
                  mkwtz = j;
                  rw = k;
                  ii = ll;
                  if (j <= count + 1)
                     break;
               }
            }
         }
      }

      else
      {
         ++count;
         continue;
      }

      removeDR(col[cl]);
      initDR(col[cl]);

      if (ii >= 0)
      {
         /*  Initialize selected pivot element
          */
         Pring *pr;
         row[rw].pos = ii - rbeg[rw];
         row[rw].mkwtz = mkwtz = (mkwtz - 1) * (count - 1);
         for (pr = pivots.next; pr->idx >= 0; pr = pr->next)
         {
            if (pr->idx == rw || pr->mkwtz >= mkwtz)
               break;
         }
         pr = pr->prev;
         if (pr->idx != rw)
         {
            removeDR(row[rw]);
            init2DR (row[rw], *pr);
         }
         num++;
         if (num >= candidates)
            break;
      }
   }

   /*
       while(pivots.next->mkwtz < pivots.prev->mkwtz)
       {
           Pring   *pr;
           pr = pivots.prev;
           removeDR(*pr);
           init2DR (*pr, rowNZ[rlen[pr->idx]]);
       }
   */

   assert(fac->row.perm[rw] < 0);
   assert(fac->col.perm[cl] < 0);
}


/*
 *      Perform L and update loop for row r
 */
static int updateRow
(
   int r,
   int lv,
   int prow,
   int pcol,
   double pval
)
{
   int fill;
   double x, lx;
   int c, i, j, k, l, m, n;

   n = rbeg[r];
   m = --(rlen[r]) + n;

   /*  compute L vector entry and
    *  and remove pivot column form row file
    */
   for (j = m; ridx[j] != pcol; --j)
     ;
   lx = rval[j] / pval;
   lval[lv] = lx;
   lidx[lv] = r;
   ++lv;

   ridx[j] = ridx[m];
   rval[j] = rval[m];


   /*  update loop (I) and
    *  computing expected fill
    */
   fill = rlen[prow];
   for (j = m - 1; j >= n; --j)
   {
      c = ridx[j];
      if (s_mark[c])
      {
         /*  count fill elements.
          */
         s_mark[c] = 0;
         --fill;

         /*  update row values
          */
         x = rval[j] -= work[c] * lx;
         if (isZero(x, epsilon))
         {
            /* Eliminate zero from row r
             */
            --rlen[r];
            --m;
            rval[j] = rval[m];
            ridx[j] = ridx[m];

            /* Eliminate zero from column c
             */
            --(s_cact[c]);
            k = --(clen[c]) + cbeg[c];
            for (i = k; cidx[i] != r; --i)
              ;
            cidx[i] = cidx[k];
         }
      }
   }


   /*  create space for fill in row file
    */
   l = rlen[r];
   if (l + fill > rmax[r])
      remaxRow(fac, r, l + fill);
   l += rbeg[r];

   /*  fill creating update loop (II)
    */
   for (j = rbeg[prow], m = j + rlen[prow]; j < m; ++j)
   {
      c = ridx[j];
      if (s_mark[c])
      {
         x = - work[c] * lx;
         if (isNonZero(x, epsilon))
         {
            /* produce fill element in row r
             */
            rval[l] = x;
            ridx[l] = c;
            l++;
            rlen[r]++;

            /* produce fill element in column c
             */
            if (clen[c] >= cmax[c])
               remaxCol(fac, c, clen[c] + 1);
            cidx[cbeg[c] + (clen[c])++] = r;
            s_cact[c]++;
         }
      }
      else
         s_mark[c] = 1;
   }

   /*  move row to appropriate list.
    */
   removeDR(row[r]);
   init2DR(row[r], rowNZ[rlen[r]]);
   assert(fac->row.perm[r] < 0);
   s_max[r] = -1;

   return lv;
}

/*
 *      Eliminate pivot element
 */
static void eliminatePivot(int prow, int pos)
{
   int i, j, k, l;
   int lv;
   int pcol;
   double pval;
   int pbeg = rbeg[prow];
   int plen = --(rlen[prow]);
   int pend = pbeg + plen;


   /*  extract pivot element   */
   i = pbeg + pos;
   pcol = ridx[i];
   pval = rval[i];
   removeDR(col[pcol]);
   initDR(col[pcol]);

   /*  remove pivot from pivot row     */
   ridx[i] = ridx[pend];
   rval[i] = rval[pend];

   /*  set pivot element and construct L vector */
   setPivot(stage++, pcol, prow, pval);
   if (s_cact[pcol] - 1 > 0)
      lv = makeLvec(fac, s_cact[pcol] - 1, prow);


   /*  init working vector,
    *  remove pivot row from working matrix
    *  and remove columns from list.
    */
   for (i = pbeg; i < pend; ++i)
   {
      j = ridx[i];
      s_mark[j] = 1;
      work[j] = rval[i];
      removeDR(col[j]);
      l = cbeg[j] + clen[j] - s_cact[j];
      for (k = l; cidx[k] != prow; ++k)
        ;
      cidx[k] = cidx[l];
      cidx[l] = prow;
      s_cact[j]--;
   }

   /*  perform L and update loop
    */
   for
   (
      i = clen[pcol] - s_cact[pcol];
      (l = cidx[cbeg[pcol] + i]) != prow;
      ++i
  )
   {
      assert(fac->row.perm[l] < 0);
      updateRow(l, lv++, prow, pcol, pval);
   }

   /*  skip pivot row  */

   l = clen[pcol];
   for (++i; i < l; ++i)
      updateRow(cidx[cbeg[pcol] + i], lv++, prow, pcol, pval);


   /*  remove pivot column from column file.
    */
   clen[pcol] -= s_cact[pcol];

   /*  clear working vector and reinsert columns to lists
    */
   for (i = rbeg[prow], pend = i + plen; i < pend; ++i)
   {
      j = ridx[i];
      work[j] = 0;
      s_mark[j] = 0;
      init2DR(col[j], colNZ[s_cact[j]]);
      assert(fac->col.perm[j] < 0);
   }
}


/*
 *      Factorize nucleus.
 */
static void eliminateNucleus(void)
{
   int r, c;
   Pring *pivot;

   pivots.mkwtz = -1;
   pivots.idx = -1;
   pivots.pos = -1;

   while (stage < dim - 1)
   {
#ifdef DEBUG
      int i;
      CLUFactorIsConsistent(fac);
      for (i = 0; i < dim; ++i)
         if (fac->col.perm[i] < 0)
         {
            assert(s_mark[i] == 0);
         }
#endif

      if (rowNZ[1].next != &(rowNZ[1]))        /* row singleton available */
         eliminateRowSingletons();
      else if (colNZ[1].next != &(colNZ[1]))   /* column singleton available */
         eliminateColSingletons();
      else
      {
         initDR(pivots);
         selectPivots();

         assert (\
                  pivots.next != &pivots && \
                  "ERROR: no pivot element selected"
               );

         for (pivot = pivots.next; pivot != &pivots; pivot = pivot->next)
         {
            eliminatePivot(pivot->idx, pivot->pos);
         }
      }

      if (rowNZ->next != rowNZ || colNZ->next != colNZ)
      {
         fac->stat = CLU_SINGULAR;
         return;
      }
   }

   if (stage < dim)
   {
      /*      Eliminate remaining element.
       *      Note, that this must be both, column and row singleton.
       */
      assert(rowNZ[1].next != &(rowNZ[1]) && "ERROR: one row must be left");
      assert(colNZ[1].next != &(colNZ[1]) && "ERROR: one col must be left");
      r = rowNZ[1].next->idx;
      c = colNZ[1].next->idx;
      rlen[r] = 0;
      clen[c]--;
      setPivot(stage, c, r, rval[rbeg[r]]);
   }
}

/*****************************************************************************/

static int setupColVals(CLUFactor* fc)
{
   int i, j, k, n;
   int* idx;
   double* val;
   double* cval;
   double maxabs;

   if (fc->u.col.val)
      free(fc->u.col.val);
   fc->u.col.val = (double*)Malloc(fc->u.col.size * sizeof(double));
   cval = fc->u.col.val;

   cval = fc->u.col.val;

   for (i = dim - 1; i >= 0; --i)
      clen[i] = 0;

   maxabs = 0;
   for (n = dim, i = dim - 1; i >= 0; --i)
   {
      j = rbeg[i];
      idx = &ridx[j];
      val = &rval[j];
      j = rlen[i];
      n += j;
      while (j--)
      {
         k = cbeg[*idx] + clen[*idx]++;
         cidx[k] = i;
         cval[k] = *val;
         if (*val > maxabs)
            maxabs = *val;
         else if (-*val > maxabs)
            maxabs = -*val;
         ++idx;
         ++val;
      }
   }
   fac->maxabs = maxabs;
   return n;
}

/*****************************************************************************/

#ifdef WITH_L_ROWS
static void setupRowVals(CLUFactor* fc)
{
   int i, j, k, l;
   int dim, vecs, mem;
   int* row;
   int* idx;
   double* val;
   int* beg;
   int* ridx;
   double* rval;
   int* rbeg;
   int *rorig, *rrorig;
   int *rperm, *rrperm;

   dim = fc->thedim;
   vecs = fc->l.firstUpdate;
   row = fc->l.row;
   idx = fc->l.idx;
   val = fc->l.val;
   beg = fc->l.start;
   mem = beg[vecs];

   if (fc->l.rval)
   {
      free(fc->l.rval);
      free(fc->l.ridx);
      free(fc->l.rbeg);
      free(fc->l.rorig);
      free(fc->l.rperm);
   }
   fc->l.rval = (double*)Malloc(mem * sizeof(double));
   fc->l.ridx = (int *)Malloc(mem * sizeof(int));
   fc->l.rbeg = (int *)Malloc((dim + 1) * sizeof(int));
   fc->l.rorig = (int *)Malloc(dim * sizeof(int));
   fc->l.rperm = (int *)Malloc(dim * sizeof(int));

   ridx = fc->l.ridx;
   rval = fc->l.rval;
   rbeg = fc->l.rbeg;
   rorig = fc->l.rorig;
   rrorig = fc->row.orig;
   rperm = fc->l.rperm;
   rrperm = fc->row.perm;

   for (i = dim; i--; *rbeg++ = 0)
   {
      *rorig++ = *rrorig++;
      *rperm++ = *rrperm++;
   }
   *rbeg = 0;

   rbeg = fc->l.rbeg + 1;
   for (i = mem; i--;)
      rbeg[*idx++]++;
   idx = fc->l.idx;

   for (l = 0, i = dim; i--; rbeg++)
   {
      j = *rbeg;
      *rbeg = l;
      l += j;
   }
   assert(l == mem);

   rbeg = fc->l.rbeg + 1;
   for (i = j = 0; i < vecs; ++i)
   {
      l = row[i];
      assert(idx == &fc->l.idx[fc->l.start[i]]);
      for (; j < beg[i + 1]; j++)
      {
         k = rbeg[*idx++]++;
         ridx[k] = l;
         rval[k] = *val++;
      }
   }
   assert(fc->l.rbeg[dim] == mem);
   assert(fc->l.rbeg[0] == 0);
}
#endif

/*****************************************************************************/

int factor(CLUFactor* fc,
            SVector** vec,          /* Array of column vector pointers   */
            double th,           /* pivoting threshold                */
            double eps           /* epsilon for zero detection        */
         )
{
   fac = fc;
   epsilon = eps;
   threshold = th;
   dim = fac->thedim;
   stage = 0;
   fac->stat = CLU_OK;

   lidx = fac->l.idx;
   lval = fac->l.val;
   lrow = fac->l.row;
   lbeg = fac->l.start;
   lbeg[0] = 0;
   fac->l.firstUpdate = 0;
   fac->l.firstUnused = 0;

   cidx = fac->u.col.idx;
   cmax = fac->u.col.max;
   clen = fac->u.col.len;
   cbeg = fac->u.col.start;

   ridx = fac->u.row.idx;
   rval = fac->u.row.val;
   rmax = fac->u.row.max;
   rlen = fac->u.row.len;
   rbeg = fac->u.row.start;

   diag = fac->diag;
   work = fac->work;

   newTmp(dim);
   initPerm();

   initMatrix(vec);
   if (fac->stat)
      goto TERMINATE;
   fac->initMaxabs = initMaxabs;

   colSingletons();
   if (fac->stat)
      goto TERMINATE;

   rowSingletons();
   if (fac->stat)
      goto TERMINATE;

   if (stage < dim)
   {
      initRings();
      eliminateNucleus();
      freeRings();
   }

TERMINATE:
   deleteTmp();
   fac->l.firstUpdate = fac->l.firstUnused;

   if (!fac->stat)
   {
#ifdef WITH_L_ROWS
      setupRowVals(fac);
#endif
      fac->nzCnt = setupColVals(fac);
   }

   /* assert(dumpCLUFactor(fac)); */

   return fac->stat;
}

void dumpCLUFactor(CLUFactor* fac)
{
   int i, j, k;

   /*  Dump U:
    */
   for (i = 0; i < fac->thedim; ++i)
   {
      if (fac->row.perm[i] >= 0)
         printf("diag[%2d]: [%2d] = %g\n",
                i, fac->col.orig[fac->row.perm[i]], fac->diag[i]);
      for (j = 0; j < fac->u.row.len[i]; ++j)
         printf
         (
            "   u[%2d]:  [%2d] = %g\n",
            i,
            fac->u.row.idx[fac->u.row.start[i] + j],
            fac->u.row.val[fac->u.row.start[i] + j]
        );
   }

   /*  Dump L:
    */
   for (i = 0; i < fac->thedim; ++i)
   {
      for (j = 0; j < fac->l.firstUnused; ++j)
         if (fac->col.orig[fac->row.perm[fac->l.row[j]]] == i)
         {
            printf("l[%2d]:\n", i);
            for (k = fac->l.start[j]; k < fac->l.start[j + 1]; ++k)
               printf
               (
                  "   l[%2d]:  [%2d] = %g\n",
                  k - fac->l.start[j],
                  fac->l.idx[k],
                  fac->l.val[k]
              );
            break;
         }
   }

   return;
}

/*****************************************************************************/
/*
 *      Perform garbage collection on row file
 */
void packRows(CLUFactor* fac)
{
   int n, i, j, row;
   Dring *ring, *list;

   int *ridx = fac->u.row.idx;
   double *rval = fac->u.row.val;
   int *rlen = fac->u.row.len;
   int *rmax = fac->u.row.max;
   int *rbeg = fac->u.row.start;

   n = 0;
   list = &(fac->u.row.list);
   for (ring = list->next; ring != list; ring = ring->next)
   {
      row = ring->idx;
      if (rbeg[row] != n)
      {
         /*@
         fprintf(stderr, "packing rows ...\n");
         */
         do
         {
            row = ring->idx;
            i = rbeg[row];
            assert(rlen[row] <= rmax[row]);
            rbeg[row] = n;
            rmax[row] = rlen[row];
            j = i + rlen[row];
            for (; i < j; ++i, ++n)
            {
               assert(n <= i);
               ridx[n] = ridx[i];
               rval[n] = rval[i];
            }
            ring = ring->next;
         }
         while (ring != list);
         goto terminatePackRows;
      }
      n += rlen[row];
      rmax[row] = rlen[row];
   }

terminatePackRows:
   fac->u.row.max[dim] = 0;
   fac->u.row.used = n;
}


/*
 *      Make row of fac large enough to hold len nonzeros.
 */
void remaxRow(CLUFactor* fc, int row, int len)
{
   fac = fc;
   assert(fac->u.row.max[row] < len);

   if (fac->u.row.elem[row].next == &(fac->u.row.list)) /* last in row file */
   {
      int delta = len - fac->u.row.max[row];

      if (delta > fac->u.row.size - fac->u.row.used)
      {
         packRows(fac);
         if (fac->u.row.size < fac->rowMemMult * fac->u.row.used + len)
            minRowMem(2 * fac->u.row.used + len);
         /* minRowMem(fac->rowMemMult * fac->u.row.used + len); */
      }
      assert(delta <= fac->u.row.size - fac->u.row.used
             && "ERROR: could not allocate memory for row file");

      fac->u.row.used += delta;
      fac->u.row.max[row] = len;
   }

   else                        /* row must be moved to end of row file */
   {
      int i, j, k;
      int *idx;
      double *val;
      Dring *ring;

      if (len > fac->u.row.size - fac->u.row.used)
      {
         packRows(fac);
         if (fac->u.row.size < fac->rowMemMult * fac->u.row.used + len)
            minRowMem(2 * fac->u.row.used + len);
         /* minRowMem(fac->rowMemMult * fac->u.row.used + len);*/
      }
      assert(len <= fac->u.row.size - fac->u.row.used
             && "ERROR: could not allocate memory for row file");

      j = fac->u.row.used;
      i = fac->u.row.start[row];
      k = fac->u.row.len[row] + i;
      fac->u.row.start[row] = j;
      fac->u.row.used += len;

      fac->u.row.max[fac->u.row.elem[row].prev->idx] += fac->u.row.max[row];
      fac->u.row.max[row] = len;
      removeDR(fac->u.row.elem[row]);
      ring = fac->u.row.list.prev;
      init2DR (fac->u.row.elem[row], *ring);

      idx = fac->u.row.idx;
      val = fac->u.row.val;
      for (; i < k; ++i, ++j)
      {
         val[j] = val[i];
         idx[j] = idx[i];
      }
   }
}


/*****************************************************************************/
/*
 *      Perform garbage collection on column file
 */
static void packColumns(CLUFactor* fac)
{
   int n, i, j, col;
   Dring *ring, *list;

   int *cidx = fac->u.col.idx;
   int *clen = fac->u.col.len;
   int *cmax = fac->u.col.max;
   int *cbeg = fac->u.col.start;

   /*@
   fprintf(stderr, "packing columns ...\n");
   */

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
               cidx[n++] = cidx[i];
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
   fac->u.col.max[dim] = 0;
}


/*
 *      Make column col of fac large enough to hold len nonzeros.
 */
void remaxCol(CLUFactor* fc, int col, int len)
{
   fac = fc;
   assert(fac->u.col.max[col] < len);

   if (fac->u.col.elem[col].next == &(fac->u.col.list)) /* last in column file */
   {
      int delta = len - fac->u.col.max[col];

      if (delta > fac->u.col.size - fac->u.col.used)
      {
         packColumns(fac);
         delta = len - fac->u.col.max[col];
         if (fac->u.col.size < fac->colMemMult * fac->u.col.used + len)
            minColMem(2 * fac->u.col.used + len);
         /* minColMem(fac->colMemMult * fac->u.col.used + len); */
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
      Dring *ring;

      if (len > fac->u.col.size - fac->u.col.used)
      {
         packColumns(fac);
         if (fac->u.col.size < fac->colMemMult * fac->u.col.used + len)
            minColMem(2 * fac->u.col.used + len);
         /* minColMem(fac->colMemMult * fac->u.col.used + len); */
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
      for (; i < k; ++i)
         idx[j++] = idx[i];
   }
}


/*****************************************************************************/

int CLUFactorIsConsistent(CLUFactor *fac)
{
   int i, j, k, l;
   Dring *ring;
   Pring *pring;

   /*  Consistency only relevant for real factorizations
    */
   if (fac->stat)
      return 1;

   /*  Test column ring list consistency.
    */
   i = 0;
   for (ring = fac->u.col.list.next; ring != &(fac->u.col.list); ring = ring->next)
   {
      assert(ring->idx >= 0);
      assert(ring->idx < dim);
      assert(ring->prev->next == ring);
      if (ring != fac->u.col.list.next)
      {
         assert
         (
            fac->u.col.start[ring->prev->idx] + fac->u.col.max[ring->prev->idx]
            == fac->u.col.start[ring->idx]
        );
      }
      ++i;
   }
   assert(i == dim);
   assert
   (
      fac->u.col.start[ring->prev->idx] + fac->u.col.max[ring->prev->idx]
      == fac->u.col.used
  );


   /*  Test row ring list consistency.
    */
   i = 0;
   for (ring = fac->u.row.list.next; ring != &(fac->u.row.list); ring = ring->next)
   {
      assert(ring->idx >= 0);
      assert(ring->idx < dim);
      assert(ring->prev->next == ring);
      assert
      (
         fac->u.row.start[ring->prev->idx] + fac->u.row.max[ring->prev->idx]
         == fac->u.row.start[ring->idx]
     );
      ++i;
   }
   assert(i == dim);
   assert
   (
      fac->u.row.start[ring->prev->idx] + fac->u.row.max[ring->prev->idx]
      == fac->u.row.used
  );


   /*  Test consistency of individual svectors.
    */
   for (i = 0; i < dim; ++i)
   {
      assert(fac->u.row.max[i] >= fac->u.row.len[i]);
      assert(fac->u.col.max[i] >= fac->u.col.len[i]);
   }


   /*  Test consistency of column file to row file of U
    */
   for (i = 0; i < dim; ++i)
   {
      for
      (
         j = fac->u.row.start[i] + fac->u.row.len[i] - 1;
         j >= fac->u.row.start[i];
         j--
     )
      {
         k = fac->u.row.idx[j];
         for
         (
            l = fac->u.col.start[k] + fac->u.col.len[k] - 1;
            l >= fac->u.col.start[k];
            l--
        )
         {
            if (fac->u.col.idx[l] == i)
               break;
         }
         assert(fac->u.col.idx[l] == i);
         if (fac->row.perm[i] < 0)
         {
            assert(fac->col.perm[k] < 0);
         }
         else
         {
            assert(fac->col.perm[k] < 0
                    || fac->col.perm[k] > fac->row.perm[i]);
         }
      }
   }

   /*  Test consistency of row file to column file of U
    */
   for (i = 0; i < dim; ++i)
   {
      for
      (
         j = fac->u.col.start[i] + fac->u.col.len[i] - 1;
         j >= fac->u.col.start[i];
         j--
     )
      {
         k = fac->u.col.idx[j];
         for
         (
            l = fac->u.row.start[k] + fac->u.row.len[k] - 1;
            l >= fac->u.row.start[k];
            l--
        )
         {
            if (fac->u.row.idx[l] == i)
               break;
         }
         assert(fac->u.row.idx[l] == i);
         assert(fac->col.perm[i] < 0
                 || fac->row.perm[k] < fac->col.perm[i]);
      }
   }

   /*  Test consistency of nonzero count lists
    */
   if (colNZ && rowNZ)
      for (i = 0; i < dim - stage; ++i)
      {
         for (pring = rowNZ[i].next; pring != &(rowNZ[i]); pring = pring->next)
         {
            assert(fac->row.perm[pring->idx] < 0);
         }
         for (pring = colNZ[i].next; pring != &(colNZ[i]); pring = pring->next)
         {
            assert(fac->col.perm[pring->idx] < 0);
         }
      }

   return 1;
}


/*****************************************************************************/

void* Malloc (int size)
{
   return malloc(size + (size <= 0));
}

void* Realloc (void* old, int size)
{
   return realloc(old, size);
}

void Free (void* old)
{
   free(old);
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
