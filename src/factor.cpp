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
#pragma ident "@(#) $Id: factor.cpp,v 1.11 2001/11/30 14:35:02 bzfbleya Exp $"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "clutypes.h"
#include "clumembers.h"
#include "cluprotos.h"
#include "cring.h"
#include "spxalloc.h"

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

/*****************************************************************************/
/*
 *      Temporary data structures.
 */

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

static CLUFactor::Pring pivots;         /* ring of selected pivot rows */
static CLUFactor::Pring *col,            /* column index handlers for double linked list */
*colNZ,          /* lists for columns to number of nonzeros      */
*row,            /* same for rows */
*rowNZ;         /* same for rows */


   /************************************************************/
   CLUFactor::Temp::Temp(int p_dim)
   {
      spx_alloc(s_max, p_dim);
      spx_alloc(s_cact, p_dim);
      spx_alloc(s_mark, p_dim);
   }
   CLUFactor::Temp::~Temp()
   {
      spx_free(s_mark);
      spx_free(s_cact);
      spx_free(s_max);
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
      spx_realloc(fac->u.row.val, size);
      spx_realloc(fac->u.row.idx, size);

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
      spx_realloc(fac->u.col.idx, size);
      cidx = fac->u.col.idx;
   }
}


/*****************************************************************************/
/*
 *      Make new Lvector to fit len nonzeros.
 *      return index to the 
 */
static void minLMem(CLUFactor *p_fac, int size)
{
   if (size > p_fac->l.size)
   {
      p_fac->l.size = int(0.2 * p_fac->l.size + size);
      spx_realloc(p_fac->l.val, p_fac->l.size);
      spx_realloc(p_fac->l.idx, p_fac->l.size);
      lidx = p_fac->l.idx;
      lval = p_fac->l.val;
   }
}

int makeLvec(CLUFactor* p_fac, int p_len, int p_row)
{
   int* p_lrow = p_fac->l.row;
   int* p_lbeg = p_fac->l.start;
   int first = p_lbeg[p_fac->l.firstUnused];

   if (p_fac->l.firstUnused >= p_fac->l.startSize)
   {
      p_fac->l.startSize += 100;
      spx_realloc(p_fac->l.start, p_fac->l.startSize);
      p_lbeg = p_fac->l.start;
   }

   assert(p_len > 0 && "ERROR: no empty columns allowed in L vectors");

   minLMem(p_fac, first + p_len);
   p_lrow[p_fac->l.firstUnused] = p_row;
   lbeg[++(p_fac->l.firstUnused)] = first + p_len;

   assert(lbeg[p_fac->l.firstUnused] <= p_fac->l.size);
   assert(p_fac->l.firstUnused <= p_fac->l.startSize);
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

static void setPivot(int p_stage, int p_col, int p_row, double val)
{
   /*@
   assert(printf("%5d:   (%3d, %3d) = %g\n", p_stage, p_row, p_col, val) || 1);
   */
   assert(fac->row.perm[p_row] < 0);
   assert(fac->col.perm[p_col] < 0);
   fac->row.orig[p_stage] = p_row;
   fac->col.orig[p_stage] = p_col;
   fac->row.perm[p_row] = p_stage;
   fac->col.perm[p_col] = p_stage;
   diag[p_row] = 1 / val;
   if (diag[p_row] > fac->maxabs)
      fac->maxabs = diag[p_row];
   else if (-diag[p_row] > fac->maxabs)
      fac->maxabs = -diag[p_row];
}


/*****************************************************************************/
/*
 *      Initialize row and column file of working matrix and
 *      mark column singletons.
 */
static void initMatrix(SVector** vec, CLUFactor::Temp& temp )
{
   double x;
   int i, j, l, k, m;
   int tot;
   Dring *rring, *lastrring;
   Dring *cring, *lastcring;
   SVector *psv;
   int *sing = temp.s_mark;

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
   rbeg[dim]       = 0;
   rmax[dim]       = 0;
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
         temp.s_cact[i] = clen[i] = cmax[i] = l;
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

static void colSingletons(CLUFactor::Temp& temp)
{
   int i, j, k, n;
   int len;
   int p_col, p_row, newrow;
   int *idx;
   int *rorig = fac->row.orig;
   int *rperm = fac->row.perm;
   int *sing = temp.s_mark;


   /*  Iteratively update column counts due to removed column singletons
    *  thereby removing new arising columns singletons
    *  and computing the index of the first row singleton (-1)
    *  until no more can be found.
    */
   fac->u.lastColSing = -1;
   for (i = 0; i < stage; ++i)
   {
      p_row = rorig[i];
      assert(p_row >= 0);
      idx = &(ridx[rbeg[p_row]]);
      len = rlen[p_row];

      if (len)
         fac->u.lastColSing = i;

      for (j = 0; j < len; ++j)
      {
         /*  Move pivotal nonzeros to front of column.
          */
         p_col = idx[j];
         n = cbeg[p_col] + clen[p_col] - temp.s_cact[p_col];
         for (k = n; cidx[k] != p_row; ++k)
           ;
         assert(k < cbeg[p_col] + clen[p_col]);
         cidx[k] = cidx[n];
         cidx[n] = p_row;

         n = --(temp.s_cact[p_col]);          /* column nonzeros of ACTIVE matrix */

         if (n == 1)                  /* Here is another singleton */
         {
            newrow = cidx[--clen[p_col] + cbeg[p_col]];

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
            for (k = n; ridx[k] != p_col; --k)
              ;

            /*      Remove singleton from column.
             */
            setPivot(stage, p_col, newrow, rval[k]);
            sing[stage++] = p_col;

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
static void rowSingletons(CLUFactor::Temp& temp)
{
   double pval;
   int i, j, k, l, r;
   int p_row, p_col, len, rs, lk;
   int *idx;
   int *rperm = fac->row.perm;
   int *sing = temp.s_mark;

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
      p_row = sing[rs];
      j = rbeg[p_row];
      p_col = ridx[j];
      pval = rval[j];
      setPivot(rs, p_col, p_row, pval);
      rlen[p_row] = 0;

      /*      Remove pivot column form workingmatrix
       *      thereby building up L vector.
       */
      idx = &(cidx[cbeg[p_col]]);
      i = temp.s_cact[p_col];                /* nr. nonzeros of new L vector */
      lk = makeLvec(fac, i - 1, p_row);
      len = clen[p_col];
      i = (clen[p_col] -= i);         /* remove pivot column from U */

      for (; i < len; ++i)
      {
         r = idx[i];
         if (r != p_row)
         {
            /*      Find pivot column in row.
             */
            l = --(rlen[r]);
            k = rbeg[r] + l;
            for (j = k; ridx[j] != p_col; --j)
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

static void initRings(CLUFactor::Temp& temp)
{
   int i;
   int *rperm = fac->row.perm;
   int *cperm = fac->col.perm;
   CLUFactor::Pring *ring;

   spx_alloc(col,   dim + 1);
   spx_alloc(colNZ, dim + 1);
   spx_alloc(row,   dim + 1);
   spx_alloc(rowNZ, dim + 1);

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
         temp.s_max[i] = -1;
      }
      if (cperm[i] < 0)
      {
         assert(temp.s_cact[i] > 1);
         ring = &(colNZ[temp.s_cact[i]]);
         init2DR(col[i], *ring);
         col[i].idx = i;
         temp.s_mark[i] = 0;
      }
   }
}

static void freeRings(void)
{
   spx_free(col);
   spx_free(colNZ);
   spx_free(row);
   spx_free(rowNZ);
}


/*****************************************************************************/

/*
 *      Eliminate all row singletons from nucleus.
 *      A row singleton may well be column singleton at the same time!
 */
static void eliminateRowSingletons(CLUFactor::Temp& temp)
{
   int i, j, k, l, r;
   int len, lk;
   int pcol, prow;
   double pval;
   int *idx;
   CLUFactor::Pring *sing;

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
      i = temp.s_cact[pcol];
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
            temp.s_max[r] = -1;
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
            temp.s_max[r] = -1;
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
static void eliminateColSingletons(CLUFactor::Temp& temp)
{
   int i, j, k, l, c;
   int pcol, prow;
   CLUFactor::Pring *sing;

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
         l = clen[c] + cbeg[c] - (temp.s_cact[c])--;
         for (k = l; cidx[k] != prow; ++k)
           ;
         cidx[k] = cidx[l];
         cidx[l] = prow;
         l = temp.s_cact[c];
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
         l = clen[c] + cbeg[c] - (temp.s_cact[c])--;
         for (k = l; cidx[k] != prow; ++k)
           ;
         cidx[k] = cidx[l];
         cidx[l] = prow;
         l = temp.s_cact[c];
         removeDR(col[c]);
         init2DR(col[c], colNZ[l]);
         assert(fac->col.perm[c] < 0);
      }
   }

   initDR(colNZ[1]);           /* Remove all column singletons from list */
}


/*
 * No singletons available: Select pivot elements.
 */
static void selectPivots(CLUFactor::Temp& temp)
{
   int ii;
   int i;
   int j;
   int k; 
   int ll = -1; // This value should never be used.
   int l;
   int m;
   int count;
   int num;
   int rw = -1; // This value should never be used.
   int cl = -1; // This value should never be used.
   int len;
   int beg;
   double maxabs;
   double x = 0.0; // This value should never be used.
   int mkwtz;
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
         if ((maxabs = temp.s_max[rw]) < 0)
         {
            maxabs = rval[len];
            if (maxabs < 0)
               maxabs = -maxabs;
            for (i = len - 1; i >= beg; --i)
               if (maxabs < rval[i])
                  maxabs = rval[i];
               else if (maxabs < -rval[i])
                  maxabs = -rval[i];
            temp.s_max[rw] = maxabs;               /* ##### */
         }
         maxabs *= threshold;

         /*  select pivot element with lowest markowitz number in row
          */
         mkwtz = dim + 1;
         for (i = len; i >= beg; --i)
         {
            k = ridx[i];
            j = temp.s_cact[k];
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
         beg = len - temp.s_cact[cl] + 1;
         assert(count == temp.s_cact[cl]);

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
               if (temp.s_max[k] > 0)
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
                  maxabs = temp.s_max[k];
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
                  temp.s_max[k] = maxabs;
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
         CLUFactor::Pring *pr;
         row[rw].pos = ii - rbeg[rw];
         row[rw].mkwtz = mkwtz = (mkwtz - 1) * (count - 1);  // ??? mkwtz originally was long, maybe to avoid an overflow in this instruction?
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
   double pval,
   CLUFactor::Temp& temp
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
      if (temp.s_mark[c])
      {
         /*  count fill elements.
          */
         temp.s_mark[c] = 0;
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
            --(temp.s_cact[c]);
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
      if (temp.s_mark[c])
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
            temp.s_cact[c]++;
         }
      }
      else
         temp.s_mark[c] = 1;
   }

   /*  move row to appropriate list.
    */
   removeDR(row[r]);
   init2DR(row[r], rowNZ[rlen[r]]);
   assert(fac->row.perm[r] < 0);
   temp.s_max[r] = -1;

   return lv;
}

/*
 *      Eliminate pivot element
 */
static void eliminatePivot(int prow, int pos, CLUFactor::Temp& temp )
{
   int i, j, k, l = -1;
   int lv = -1;  // This value should never be used.
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

   /**@todo If this test failes, lv has no value. I suppose that in this
    *       case none of the loops below that uses lv is executed.
    *       But this is unproven.
    */
   if (temp.s_cact[pcol] - 1 > 0)
      lv = makeLvec(fac, temp.s_cact[pcol] - 1, prow);

   /*  init working vector,
    *  remove pivot row from working matrix
    *  and remove columns from list.
    */
   for (i = pbeg; i < pend; ++i)
   {
      j = ridx[i];
      temp.s_mark[j] = 1;
      work[j] = rval[i];
      removeDR(col[j]);
      l = cbeg[j] + clen[j] - temp.s_cact[j];
      for (k = l; cidx[k] != prow; ++k)
        ;
      cidx[k] = cidx[l];
      cidx[l] = prow;
      temp.s_cact[j]--;
   }

   /*  perform L and update loop
    */
   for
   (
      i = clen[pcol] - temp.s_cact[pcol];
      (l = cidx[cbeg[pcol] + i]) != prow;
      ++i
  )
   {
      assert(fac->row.perm[l] < 0);
      assert(lv >= 0);
      updateRow(l, lv++, prow, pcol, pval, temp);
   }

   /*  skip pivot row  */

   l = clen[pcol];
   for (++i; i < l; ++i)
   {
      assert(lv >= 0);
      updateRow(cidx[cbeg[pcol] + i], lv++, prow, pcol, pval, temp);
   }

   /*  remove pivot column from column file.
    */
   clen[pcol] -= temp.s_cact[pcol];

   /*  clear working vector and reinsert columns to lists
    */
   for (i = rbeg[prow], pend = i + plen; i < pend; ++i)
   {
      j = ridx[i];
      work[j] = 0;
      temp.s_mark[j] = 0;
      init2DR(col[j], colNZ[temp.s_cact[j]]);
      assert(fac->col.perm[j] < 0);
   }
}


/*
 *      Factorize nucleus.
 */
static void eliminateNucleus(CLUFactor::Temp& temp)
{
   int r, c;
   CLUFactor::Pring *pivot;

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
         eliminateRowSingletons(temp);
      else if (colNZ[1].next != &(colNZ[1]))   /* column singleton available */
         eliminateColSingletons(temp);
      else
      {
         initDR(pivots);
         selectPivots(temp);

         assert (\
                  pivots.next != &pivots && \
                  "ERROR: no pivot element selected"
               );

         for (pivot = pivots.next; pivot != &pivots; pivot = pivot->next)
         {
            eliminatePivot(pivot->idx, pivot->pos,temp);
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
      spx_free(fc->u.col.val);

   spx_alloc(fc->u.col.val, fc->u.col.size);

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
   int l_dim, vecs, mem;
   int* l_row;
   int* idx;
   double* val;
   int* beg;
   int* l_ridx;
   double* l_rval;
   int* l_rbeg;
   int *rorig, *rrorig;
   int *rperm, *rrperm;

   l_dim = fc->thedim;
   vecs  = fc->l.firstUpdate;
   l_row = fc->l.row;
   idx   = fc->l.idx;
   val   = fc->l.val;
   beg   = fc->l.start;
   mem   = beg[vecs];

   if (fc->l.rval)
   {
      spx_free(fc->l.rval);
      spx_free(fc->l.ridx);
      spx_free(fc->l.rbeg);
      spx_free(fc->l.rorig);
      spx_free(fc->l.rperm);
   }
   spx_alloc(fc->l.rval, mem);
   spx_alloc(fc->l.ridx, mem);
   spx_alloc(fc->l.rbeg, l_dim + 1);
   spx_alloc(fc->l.rorig, l_dim);
   spx_alloc(fc->l.rperm, l_dim);

   l_ridx = fc->l.ridx;
   l_rval = fc->l.rval;
   l_rbeg = fc->l.rbeg;
   rorig  = fc->l.rorig;
   rrorig = fc->row.orig;
   rperm  = fc->l.rperm;
   rrperm = fc->row.perm;

   for (i = l_dim; i--; *l_rbeg++ = 0)
   {
      *rorig++ = *rrorig++;
      *rperm++ = *rrperm++;
   }
   *l_rbeg = 0;

   l_rbeg = fc->l.rbeg + 1;
   for (i = mem; i--;)
      l_rbeg[*idx++]++;
   idx = fc->l.idx;

   for (l = 0, i = l_dim; i--; l_rbeg++)
   {
      j = *l_rbeg;
      *l_rbeg = l;
      l += j;
   }
   assert(l == mem);

   l_rbeg = fc->l.rbeg + 1;
   for (i = j = 0; i < vecs; ++i)
   {
      l = l_row[i];
      assert(idx == &fc->l.idx[fc->l.start[i]]);
      for (; j < beg[i + 1]; j++)
      {
         k = l_rbeg[*idx++]++;
         l_ridx[k] = l;
         l_rval[k] = *val++;
      }
   }
   assert(fc->l.rbeg[l_dim] == mem);
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

   CLUFactor::Temp temp(dim);
   initPerm();

   initMatrix(vec, temp);
   if (fac->stat)
      goto TERMINATE;
   fac->initMaxabs = initMaxabs;

   colSingletons(temp);
   if (fac->stat)
      goto TERMINATE;

   rowSingletons(temp);
   if (fac->stat)
      goto TERMINATE;

   if (stage < dim)
   {
      initRings(temp);
      eliminateNucleus(temp);
      freeRings();
   }

TERMINATE:
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

void dumpCLUFactor(const CLUFactor* p_fac)
{
   int i, j, k;

   /*  Dump U:
    */
   for (i = 0; i < p_fac->thedim; ++i)
   {
      if (p_fac->row.perm[i] >= 0)
         printf("diag[%2d]: [%2d] = %g\n",
                i, p_fac->col.orig[p_fac->row.perm[i]], p_fac->diag[i]);
      for (j = 0; j < p_fac->u.row.len[i]; ++j)
         printf
         (
            "   u[%2d]:  [%2d] = %g\n",
            i,
            p_fac->u.row.idx[p_fac->u.row.start[i] + j],
            p_fac->u.row.val[p_fac->u.row.start[i] + j]
        );
   }

   /*  Dump L:
    */
   for (i = 0; i < p_fac->thedim; ++i)
   {
      for (j = 0; j < p_fac->l.firstUnused; ++j)
         if (p_fac->col.orig[p_fac->row.perm[p_fac->l.row[j]]] == i)
         {
            printf("l[%2d]:\n", i);
            for (k = p_fac->l.start[j]; k < p_fac->l.start[j + 1]; ++k)
               printf
               (
                  "   l[%2d]:  [%2d] = %g\n",
                  k - p_fac->l.start[j],
                  p_fac->l.idx[k],
                  p_fac->l.val[k]
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
void packRows(CLUFactor* p_fac)
{
   int n, i, j, l_row;
   Dring *ring, *list;

   int *l_ridx = p_fac->u.row.idx;
   double *l_rval = p_fac->u.row.val;
   int *l_rlen = p_fac->u.row.len;
   int *l_rmax = p_fac->u.row.max;
   int *l_rbeg = p_fac->u.row.start;

   n = 0;
   list = &(p_fac->u.row.list);
   for (ring = list->next; ring != list; ring = ring->next)
   {
      l_row = ring->idx;
      if (l_rbeg[l_row] != n)
      {
         /*@
         fprintf(stderr, "packing rows ...\n");
         */
         do
         {
            l_row = ring->idx;
            i = l_rbeg[l_row];
            assert(l_rlen[l_row] <= l_rmax[l_row]);
            l_rbeg[l_row] = n;
            l_rmax[l_row] = l_rlen[l_row];
            j = i + l_rlen[l_row];
            for (; i < j; ++i, ++n)
            {
               assert(n <= i);
               l_ridx[n] = l_ridx[i];
               l_rval[n] = l_rval[i];
            }
            ring = ring->next;
         }
         while (ring != list);
         goto terminatePackRows;
      }
      n += l_rlen[l_row];
      l_rmax[l_row] = l_rlen[l_row];
   }

terminatePackRows:
   p_fac->u.row.max[dim] = 0;
   p_fac->u.row.used = n;
}


/*
 *      Make row of fac large enough to hold len nonzeros.
 */
void remaxRow(CLUFactor* fc, int p_row, int len)
{
   fac = fc;
   assert(fac->u.row.max[p_row] < len);

   if (fac->u.row.elem[p_row].next == &(fac->u.row.list)) /* last in row file */
   {
      int delta = len - fac->u.row.max[p_row];

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
      fac->u.row.max[p_row] = len;
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
      i = fac->u.row.start[p_row];
      k = fac->u.row.len[p_row] + i;
      fac->u.row.start[p_row] = j;
      fac->u.row.used += len;

      fac->u.row.max[fac->u.row.elem[p_row].prev->idx] += fac->u.row.max[p_row];
      fac->u.row.max[p_row] = len;
      removeDR(fac->u.row.elem[p_row]);
      ring = fac->u.row.list.prev;
      init2DR (fac->u.row.elem[p_row], *ring);

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
static void packColumns(CLUFactor* p_fac)
{
   int n, i, j, l_col;
   Dring *ring, *list;

   int *l_cidx = p_fac->u.col.idx;
   int *l_clen = p_fac->u.col.len;
   int *l_cmax = p_fac->u.col.max;
   int *l_cbeg = p_fac->u.col.start;

   /*@
   fprintf(stderr, "packing columns ...\n");
   */

   n = 0;
   list = &(p_fac->u.col.list);
   for (ring = list->next; ring != list; ring = ring->next)
   {
      l_col = ring->idx;
      if (l_cbeg[l_col] != n)
      {
         do
         {
            l_col = ring->idx;
            i = l_cbeg[l_col];
            l_cbeg[l_col] = n;
            l_cmax[l_col] = l_clen[l_col];
            j = i + l_clen[l_col];
            for (; i < j; ++i)
               l_cidx[n++] = l_cidx[i];
            ring = ring->next;
         }
         while (ring != list);
         goto terminatePackColumns;
      }
      n += l_clen[l_col];
      l_cmax[l_col] = l_clen[l_col];
   }

terminatePackColumns :
   p_fac->u.col.used = n;
   p_fac->u.col.max[dim] = 0;
}


/*
 *      Make column col of fac large enough to hold len nonzeros.
 */
void remaxCol(CLUFactor* fc, int p_col, int len)
{
   fac = fc;
   assert(fac->u.col.max[p_col] < len);

   if (fac->u.col.elem[p_col].next == &(fac->u.col.list)) /* last in column file */
   {
      int delta = len - fac->u.col.max[p_col];

      if (delta > fac->u.col.size - fac->u.col.used)
      {
         packColumns(fac);
         delta = len - fac->u.col.max[p_col];
         if (fac->u.col.size < fac->colMemMult * fac->u.col.used + len)
            minColMem(2 * fac->u.col.used + len);
         /* minColMem(fac->colMemMult * fac->u.col.used + len); */
      }
      assert(delta <= fac->u.col.size - fac->u.col.used
             && "ERROR: could not allocate memory for column file");

      fac->u.col.used += delta;
      fac->u.col.max[p_col] = len;
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
      i = fac->u.col.start[p_col];
      k = fac->u.col.len[p_col] + i;
      fac->u.col.start[p_col] = j;
      fac->u.col.used += len;

      fac->u.col.max[fac->u.col.elem[p_col].prev->idx] += fac->u.col.max[p_col];
      fac->u.col.max[p_col] = len;
      removeDR(fac->u.col.elem[p_col]);
      ring = fac->u.col.list.prev;
      init2DR (fac->u.col.elem[p_col], *ring);

      idx = fac->u.col.idx;
      for (; i < k; ++i)
         idx[j++] = idx[i];
   }
}


/*****************************************************************************/

int CLUFactorIsConsistent(const CLUFactor *p_fac)
{
   int i, j, k, l;
   Dring *ring;
   CLUFactor::Pring *pring;

   /*  Consistency only relevant for real factorizations
    */
   if (p_fac->stat)
      return 1;

   /*  Test column ring list consistency.
    */
   i = 0;
   for (ring = p_fac->u.col.list.next; ring != &(p_fac->u.col.list); ring = ring->next)
   {
      assert(ring->idx >= 0);
      assert(ring->idx < dim);
      assert(ring->prev->next == ring);
      if (ring != p_fac->u.col.list.next)
      {
         assert
         (
            p_fac->u.col.start[ring->prev->idx] + p_fac->u.col.max[ring->prev->idx]
            == p_fac->u.col.start[ring->idx]
        );
      }
      ++i;
   }
   assert(i == dim);
   assert
   (
      p_fac->u.col.start[ring->prev->idx] + p_fac->u.col.max[ring->prev->idx]
      == p_fac->u.col.used
  );


   /*  Test row ring list consistency.
    */
   i = 0;
   for (ring = p_fac->u.row.list.next; ring != &(p_fac->u.row.list); ring = ring->next)
   {
      assert(ring->idx >= 0);
      assert(ring->idx < dim);
      assert(ring->prev->next == ring);
      assert
      (
         p_fac->u.row.start[ring->prev->idx] + p_fac->u.row.max[ring->prev->idx]
         == p_fac->u.row.start[ring->idx]
     );
      ++i;
   }
   assert(i == dim);
   assert
   (
      p_fac->u.row.start[ring->prev->idx] + p_fac->u.row.max[ring->prev->idx]
      == p_fac->u.row.used
  );


   /*  Test consistency of individual svectors.
    */
   for (i = 0; i < dim; ++i)
   {
      assert(p_fac->u.row.max[i] >= p_fac->u.row.len[i]);
      assert(p_fac->u.col.max[i] >= p_fac->u.col.len[i]);
   }


   /*  Test consistency of column file to row file of U
    */
   for (i = 0; i < dim; ++i)
   {
      for
      (
         j = p_fac->u.row.start[i] + p_fac->u.row.len[i] - 1;
         j >= p_fac->u.row.start[i];
         j--
     )
      {
         k = p_fac->u.row.idx[j];
         for
         (
            l = p_fac->u.col.start[k] + p_fac->u.col.len[k] - 1;
            l >= p_fac->u.col.start[k];
            l--
        )
         {
            if (p_fac->u.col.idx[l] == i)
               break;
         }
         assert(p_fac->u.col.idx[l] == i);
         if (p_fac->row.perm[i] < 0)
         {
            assert(p_fac->col.perm[k] < 0);
         }
         else
         {
            assert(p_fac->col.perm[k] < 0
                    || p_fac->col.perm[k] > p_fac->row.perm[i]);
         }
      }
   }

   /*  Test consistency of row file to column file of U
    */
   for (i = 0; i < dim; ++i)
   {
      for
      (
         j = p_fac->u.col.start[i] + p_fac->u.col.len[i] - 1;
         j >= p_fac->u.col.start[i];
         j--
     )
      {
         k = p_fac->u.col.idx[j];
         for
         (
            l = p_fac->u.row.start[k] + p_fac->u.row.len[k] - 1;
            l >= p_fac->u.row.start[k];
            l--
        )
         {
            if (p_fac->u.row.idx[l] == i)
               break;
         }
         assert(p_fac->u.row.idx[l] == i);
         assert(p_fac->col.perm[i] < 0
                 || p_fac->row.perm[k] < p_fac->col.perm[i]);
      }
   }

   /*  Test consistency of nonzero count lists
    */
   if (colNZ && rowNZ)
      for (i = 0; i < dim - stage; ++i)
      {
         for (pring = rowNZ[i].next; pring != &(rowNZ[i]); pring = pring->next)
         {
            assert(p_fac->row.perm[pring->idx] < 0);
         }
         for (pring = colNZ[i].next; pring != &(colNZ[i]); pring = pring->next)
         {
            assert(p_fac->col.perm[pring->idx] < 0);
         }
      }

   return 1;
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
