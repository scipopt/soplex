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
#pragma ident "@(#) $Id: factor.cpp,v 1.15 2001/12/04 19:28:20 bzfkocht Exp $"

#include <iostream>
#include <assert.h>

#include "clufactor.h"
#include "cring.h"
#include "spxalloc.h"

namespace soplex
{
/************************************************************/
CLUFactor::Temp::Temp()
   : s_mark(0)
   , s_max(0)
   , s_cact(0)
{}

void CLUFactor::Temp::init(int p_dim)
{
   spx_alloc(s_max, p_dim);
   spx_alloc(s_cact, p_dim);
   spx_alloc(s_mark, p_dim);
}
   
void CLUFactor::Temp::clear()
{
   if (s_mark != 0)
      spx_free(s_mark);
   if (s_cact != 0)
      spx_free(s_cact);
   if (s_max != 0)
      spx_free(s_max);
}
   
CLUFactor::Temp::~Temp()
{
   clear();
}

/************************************************************/   

CLUFactor::Pivots CLUFactor::pivots;

/************************************************************/   
void CLUFactor::initPerm()
{
   int* const p_or = row.orig;
   int* const p_pr = row.perm;
   int* const p_oc = col.orig;
   int* const p_pc = col.perm;
   for (int i = 0; i < thedim; ++i)
      p_or[i] = p_pr[i] = p_oc[i] = p_pc[i] = -1;
}

/*****************************************************************************/
   
void CLUFactor::setPivot(int p_stage, int p_col, int p_row, double val)
{
   // std::cout << p_stage << ": (" << p_row ", " 
   //           << p_col << ") = " << val << std::endl;

   assert(row.perm[p_row] < 0);
   assert(col.perm[p_col] < 0);

   row.orig[p_stage] = p_row;
   col.orig[p_stage] = p_col;
   row.perm[p_row]   = p_stage;
   col.perm[p_col]   = p_stage;
   diag[p_row]       = 1.0 / val;

   if (fabs(diag[p_row]) > maxabs)
      maxabs = fabs(diag[p_row]);
}

/*****************************************************************************/
/*
 *      Perform garbage collection on row file
 */
void CLUFactor::packRows()
{
   int n, i, j, l_row;
   Dring *ring, *list;

   int *l_ridx = u.row.idx;
   double *l_rval = u.row.val;
   int *l_rlen = u.row.len;
   int *l_rmax = u.row.max;
   int *l_rbeg = u.row.start;

   n = 0;
   list = &(u.row.list);
   for (ring = list->next; ring != list; ring = ring->next)
   {
      l_row = ring->idx;
      if (l_rbeg[l_row] != n)
      {
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
   u.row.max[thedim] = 0;
   u.row.used = n;
}

/*
 *      Make row of fac large enough to hold len nonzeros.
 */
void CLUFactor::remaxRow(int p_row, int len)
{
   assert(u.row.max[p_row] < len);

   if (u.row.elem[p_row].next == &(u.row.list)) /* last in row file */
   {
      int delta = len - u.row.max[p_row];

      if (delta > u.row.size - u.row.used)
      {
         packRows();
         if (u.row.size < rowMemMult * u.row.used + len)
            minRowMem(2 * u.row.used + len);
         /* minRowMem(rowMemMult * u.row.used + len); */
      }
      assert(delta <= u.row.size - u.row.used
         && "ERROR: could not allocate memory for row file");

      u.row.used += delta;
      u.row.max[p_row] = len;
   }

   else                        /* row must be moved to end of row file */
   {
      int i, j, k;
      int *idx;
      double *val;
      Dring *ring;

      if (len > u.row.size - u.row.used)
      {
         packRows();
         if (u.row.size < rowMemMult * u.row.used + len)
            minRowMem(2 * u.row.used + len);
         /* minRowMem(rowMemMult * u.row.used + len);*/
      }
      assert(len <= u.row.size - u.row.used
         && "ERROR: could not allocate memory for row file");

      j = u.row.used;
      i = u.row.start[p_row];
      k = u.row.len[p_row] + i;
      u.row.start[p_row] = j;
      u.row.used += len;

      u.row.max[u.row.elem[p_row].prev->idx] += u.row.max[p_row];
      u.row.max[p_row] = len;
      removeDR(u.row.elem[p_row]);
      ring = u.row.list.prev;
      init2DR (u.row.elem[p_row], *ring);

      idx = u.row.idx;
      val = u.row.val;
      for (; i < k; ++i, ++j)
      {
         val[j] = val[i];
         idx[j] = idx[i];
      }
   }
}

/*************************************************************************/
/*
 *      Perform garbage collection on column file
 */
void CLUFactor::packColumns()
{
   int n, i, j, l_col;
   Dring *ring, *list;
      
   int *l_cidx = u.col.idx;
   int *l_clen = u.col.len;
   int *l_cmax = u.col.max;
   int *l_cbeg = u.col.start;
      
   n = 0;
   list = &(u.col.list);
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
         
      u.col.used = n;
   u.col.max[thedim] = 0;
}

/*
    *      Make column col of fac large enough to hold len nonzeros.
    */
void CLUFactor::remaxCol(int p_col, int len)
{
   assert(u.col.max[p_col] < len);
      
   if (u.col.elem[p_col].next == &(u.col.list)) /* last in column file */
   {
      int delta = len - u.col.max[p_col];
            
      if (delta > u.col.size - u.col.used)
      {
         packColumns();
         delta = len - u.col.max[p_col];
         if (u.col.size < colMemMult * u.col.used + len)
            minColMem(2 * u.col.used + len);
         /* minColMem(colMemMult * u.col.used + len); */
      }
      assert(delta <= u.col.size - u.col.used
         && "ERROR: could not allocate memory for column file");
            
      u.col.used += delta;
      u.col.max[p_col] = len;
   }
      
   else                        /* column must be moved to end of column file */
   {
      int i, j, k;
      int *idx;
      Dring *ring;
            
      if (len > u.col.size - u.col.used)
      {
         packColumns();
         if (u.col.size < colMemMult * u.col.used + len)
            minColMem(2 * u.col.used + len);
         /* minColMem(colMemMult * u.col.used + len); */
      }
      assert(len <= u.col.size - u.col.used
         && "ERROR: could not allocate memory for column file");
            
      j = u.col.used;
      i = u.col.start[p_col];
      k = u.col.len[p_col] + i;
      u.col.start[p_col] = j;
      u.col.used += len;
            
      u.col.max[u.col.elem[p_col].prev->idx] += u.col.max[p_col];
      u.col.max[p_col] = len;
      removeDR(u.col.elem[p_col]);
      ring = u.col.list.prev;
      init2DR (u.col.elem[p_col], *ring);
            
      idx = u.col.idx;
      for (; i < k; ++i)
         idx[j++] = idx[i];
   }
}

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

 /*****************************************************************************/
 /*
  *      Initialize row and column file of working matrix and
  *      mark column singletons.
  */
void CLUFactor::initFactorMatrix(SVector** vec, 
   const double eps, 
   int& stage )
{
   double x;
   int i, j, ll, k, m;
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
   for (i = 0; i < thedim; ++i)
      u.row.max[i] = u.row.len[i] = 0;

   tot = 0;

   for (i = 0; i < thedim; ++i)
   {
      psv = vec[i];
      ll = psv->size();
      if (ll > 1)
      {
         tot += ll;
         for (j = 0; j < ll; ++j)
            ++(u.row.max[psv->index(j)]);
      }
      else if (ll == 0)
      {
         stat = SLinSolver::SINGULAR;
         return;
      }
   }



   /*  Resize nonzero memory if necessary
    */
   minRowMem(int(rowMemMult * tot));
   minColMem(int(colMemMult * tot));
   minLMem  (int(lMemMult * tot));


   /*  Initialize:
    *  - row ring lists
    *  - row vectors in file
    *  - column ring lists
    */
   u.row.start[0] = 0;

   rring = u.row.elem;
   lastrring = &(u.row.list);
   lastrring->idx = thedim;
   lastrring->next = rring;

   cring = u.col.elem;
   lastcring = &(u.col.list);
   lastcring->idx = thedim;
   lastcring->next = cring;

   m = 0;
   for (i = 0; i < thedim; ++i)
   {
      u.row.start[i] = m;
      m += u.row.max[i];

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
   u.row.start[thedim]       = 0;
   u.row.max[thedim]       = 0;
   u.row.used = m;

   lastrring->next = &(u.row.list);
   lastrring->next->prev = lastrring;

   lastcring->next = &(u.col.list);
   lastcring->next->prev = lastcring;

   /*  Copy matrix to row and column file
    *  excluding and marking column singletons!
    */
   m = 0;
   stage = 0;

   initMaxabs = 0;
   for (i = 0; i < thedim; ++i)
   {
      psv = vec[i];
      ll = psv->size();
      u.col.start[i] = m;
      if (ll > 1)                               /* exclude column singletons */
      {
         int kk, lll;
         for (j = lll = 0; j < ll; ++j)
         {
            x = psv->value(j);
            if (isNonZero(x, eps))
            {
               k = psv->index(j);
               kk = u.row.start[k] + (u.row.len[k]++);
               u.col.idx[m++] = k;
               u.row.idx[kk] = i;
               u.row.val[kk] = x;
               ++lll;
               if (x > initMaxabs)
                  initMaxabs = x;
               else if (-x > initMaxabs)
                  initMaxabs = -x;
            }
         }
         ll = lll;
         --m;
      }
      if (ll > 1)
      {
         ++m;
         temp.s_cact[i] = u.col.len[i] = u.col.max[i] = ll;
      }
      else if (ll <= 0)       /* singular */
      {
         stat = SLinSolver::SINGULAR;
         return;
      }
      else                    /* singleton */
      {
         u.col.len[i] = u.col.max[i] = 0;
         for (j = 0; isZero(psv->value(j), eps); ++j)
            ;
         if (row.perm[psv->index(j)] >= 0)
         {
            stat = SLinSolver::SINGULAR;
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

   u.col.used = m;
}



/*****************************************************************************/
/*
 *      Remove column singletons
 */

void CLUFactor::colSingletons( int& stage )
{
   int i, j, k, n;
   int len;
   int p_col, p_row, newrow;
   int *idx;
   int *rorig = row.orig;
   int *rperm = row.perm;
   int *sing = temp.s_mark;


   /*  Iteratively update column counts due to removed column singletons
    *  thereby removing new arising columns singletons
    *  and computing the index of the first row singleton (-1)
    *  until no more can be found.
    */
   u.lastColSing = -1;
   for (i = 0; i < stage; ++i)
   {
      p_row = rorig[i];
      assert(p_row >= 0);
      idx = &(u.row.idx[u.row.start[p_row]]);
      len = u.row.len[p_row];

      if (len)
         u.lastColSing = i;

      for (j = 0; j < len; ++j)
      {
         /*  Move pivotal nonzeros to front of column.
          */
         p_col = idx[j];
         n = u.col.start[p_col] + u.col.len[p_col] - temp.s_cact[p_col];
         for (k = n; u.col.idx[k] != p_row; ++k)
            ;
         assert(k < u.col.start[p_col] + u.col.len[p_col]);
         u.col.idx[k] = u.col.idx[n];
         u.col.idx[n] = p_row;

         n = --(temp.s_cact[p_col]);          /* column nonzeros of ACTIVE matrix */

         if (n == 1)                  /* Here is another singleton */
         {
            newrow = u.col.idx[--u.col.len[p_col] + u.col.start[p_col]];

            /*      Ensure, matrix not singular
             */
            if (rperm[newrow] >= 0)
            {
               stat = SLinSolver::SINGULAR;
               return;
            }

            /*      Find singleton in row.
             */
            n = u.row.start[newrow] + (--(u.row.len[newrow]));
            for (k = n; u.row.idx[k] != p_col; --k)
               ;

            /*      Remove singleton from column.
             */
            setPivot(stage, p_col, newrow, u.row.val[k]);
            sing[stage++] = p_col;

            /*      Move pivot element to diag.
             */
            u.row.val[k] = u.row.val[n];
            u.row.idx[k] = u.row.idx[n];
         }
         else if (n == 0)
         {
            stat = SLinSolver::SINGULAR;
            return;
         }
      }
   }

   assert(stage <= thedim);
}
 

/*****************************************************************************/
/*
 *      Remove row singletons
 */
void CLUFactor::rowSingletons( int& stage )
{
   double pval;
   int i, j, k, ll, r;
   int p_row, p_col, len, rs, lk;
   int *idx;
   int *rperm = row.perm;
   int *sing = temp.s_mark;

   /*  Mark row singletons
    */
   rs = stage;
   for (i = 0; i < thedim; ++i)
   {
      if (rperm[i] < 0 && u.row.len[i] == 1)
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
      j = u.row.start[p_row];
      p_col = u.row.idx[j];
      pval = u.row.val[j];
      setPivot(rs, p_col, p_row, pval);
      u.row.len[p_row] = 0;

      /*      Remove pivot column form workingmatrix
       *      thereby building up L vector.
       */
      idx = &(u.col.idx[u.col.start[p_col]]);
      i = temp.s_cact[p_col];                /* nr. nonzeros of new L vector */
      lk = makeLvec(i - 1, p_row);
      len = u.col.len[p_col];
      i = (u.col.len[p_col] -= i);         /* remove pivot column from U */

      for (; i < len; ++i)
      {
         r = idx[i];
         if (r != p_row)
         {
            /*      Find pivot column in row.
             */
            ll = --(u.row.len[r]);
            k = u.row.start[r] + ll;
            for (j = k; u.row.idx[j] != p_col; --j)
               ;
            assert(k >= u.row.start[r]);

            /*      Initialize L vector
             */
            l.idx[lk] = r;
            l.val[lk] = u.row.val[j] / pval;
            ++lk;

            /*      Remove pivot column from row.
             */
            u.row.idx[j] = u.row.idx[k];
            u.row.val[j] = u.row.val[k];

            /*      Check new row length.
             */
            if (ll == 1)
               sing[stage++] = r;
            else if (ll == 0)
            {
               stat = SLinSolver::SINGULAR;
               return;
            }
         }
      }
   }

   u.lastRowSing = stage - 1;
}


/*****************************************************************************/
/*
 *      Init nonzero number Ring lists
 *      and required entries of arrays max and mark
 */

void CLUFactor::initFactorRings( const int stage )
{
   int i;
   int *rperm = row.perm;
   int *cperm = col.perm;
   CLUFactor::Pring *ring;

   spx_alloc(pivots.pivot_col,   thedim + 1);
   spx_alloc(pivots.pivot_colNZ, thedim + 1);
   spx_alloc(pivots.pivot_row,   thedim + 1);
   spx_alloc(pivots.pivot_rowNZ, thedim + 1);

   for (i = thedim - stage; i >= 0; --i)
   {
      initDR(pivots.pivot_colNZ[i]);
      initDR(pivots.pivot_rowNZ[i]);
   }

   for (i = 0; i < thedim; ++i)
   {
      if (rperm[i] < 0)
      {
         assert(u.row.len[i] > 1);
         ring = &(pivots.pivot_rowNZ[u.row.len[i]]);
         init2DR(pivots.pivot_row[i], *ring);
         pivots.pivot_row[i].idx = i;
         temp.s_max[i] = -1;
      }
      if (cperm[i] < 0)
      {
         assert(temp.s_cact[i] > 1);
         ring = &(pivots.pivot_colNZ[temp.s_cact[i]]);
         init2DR(pivots.pivot_col[i], *ring);
         pivots.pivot_col[i].idx = i;
         temp.s_mark[i] = 0;
      }
   }
}

void CLUFactor::freeFactorRings(void)
{
   spx_free(pivots.pivot_col);
   spx_free(pivots.pivot_colNZ);
   spx_free(pivots.pivot_row);
   spx_free(pivots.pivot_rowNZ);
}
   

/*****************************************************************************/
   
   /*
    *      Eliminate all row singletons from nucleus.
    *      A row singleton may well be column singleton at the same time!
    */
void CLUFactor::eliminateRowSingletons( int& stage )
{
   int i, j, k, ll, r;
   int len, lk;
   int pcol, prow;
   double pval;
   int *idx;
   CLUFactor::Pring *sing;

   for (sing = pivots.pivot_rowNZ[1].prev; sing != &(pivots.pivot_rowNZ[1]); sing = sing->prev)
   {
      prow = sing->idx;
      i = u.row.start[prow];
      pcol = u.row.idx[i];
      pval = u.row.val[i];
      setPivot(stage++, pcol, prow, pval);
      u.row.len[prow] = 0;
      removeDR(pivots.pivot_col[pcol]);

      /*      Eliminate pivot column and build L vector.
       */
      i = temp.s_cact[pcol];
      if (i > 1)
      {
         idx = &(u.col.idx[u.col.start[pcol]]);
         len = u.col.len[pcol];
         lk = makeLvec(i - 1, prow);
         i = u.col.len[pcol] -= i;

         for (; (r = idx[i]) != prow; ++i)
         {
            /*      Find pivot column in row.
             */
            ll = --(u.row.len[r]);
            k = u.row.start[r] + ll;
            for (j = k; u.row.idx[j] != pcol; --j)
               ;
            assert(j >= u.row.start[r]);

            /*      Initialize L vector
             */
            l.idx[lk] = r;
            l.val[lk] = u.row.val[j] / pval;
            ++lk;

            /*      Remove pivot column from row.
             */
            u.row.idx[j] = u.row.idx[k];
            u.row.val[j] = u.row.val[k];

            /*      Move column to appropriate nonzero ring.
             */
            removeDR(pivots.pivot_row[r]);
            init2DR (pivots.pivot_row[r], pivots.pivot_rowNZ[ll]);
            assert(row.perm[r] < 0);
            temp.s_max[r] = -1;
         }

         /* skip pivot element */
         assert(i < len && "ERROR: pivot column does not contain pivot row");

         for (++i; i < len; ++i)
         {
            /*      Find pivot column in row.
             */
            r = idx[i];
            ll = --(u.row.len[r]);
            k = u.row.start[r] + ll;
            for (j = k; u.row.idx[j] != pcol; --j)
               ;
            assert(j >= u.row.start[r]);

            /*      Initialize L vector
             */
            l.idx[lk] = r;
            l.val[lk] = u.row.val[j] / pval;
            ++lk;

            /*      Remove pivot column from row.
             */
            u.row.idx[j] = u.row.idx[k];
            u.row.val[j] = u.row.val[k];

            /*      Move column to appropriate nonzero ring.
             */
            removeDR(pivots.pivot_row[r]);
            init2DR (pivots.pivot_row[r], pivots.pivot_rowNZ[ll]);
            assert(row.perm[r] < 0);
            temp.s_max[r] = -1;
         }
      }
      else
         u.col.len[pcol] -= i;
   }

   initDR(pivots.pivot_rowNZ[1]);           /* Remove all row singletons from list */
}



/*
    *      Eliminate all column singletons from nucleus.
    *      A column singleton must not be row singleton at the same time!
    */
void CLUFactor::eliminateColSingletons( int& stage)
{
   int i, j, k, l, c;
   int pcol, prow;
   CLUFactor::Pring *sing;

   for (sing = pivots.pivot_colNZ[1].prev; sing != &(pivots.pivot_colNZ[1]); sing = sing->prev)
   {
      /*      Find pivot value
       */
      pcol = sing->idx;
      j = --(u.col.len[pcol]) + u.col.start[pcol];     /* remove pivot column */
      prow = u.col.idx[j];
      removeDR(pivots.pivot_row[prow]);

      j = --(u.row.len[prow]) + u.row.start[prow];
      for (i = j; (c = u.row.idx[i]) != pcol; --i)
      {
         l = u.col.len[c] + u.col.start[c] - (temp.s_cact[c])--;
         for (k = l; u.col.idx[k] != prow; ++k)
            ;
         u.col.idx[k] = u.col.idx[l];
         u.col.idx[l] = prow;
         l = temp.s_cact[c];
         removeDR(pivots.pivot_col[c]);
         init2DR(pivots.pivot_col[c], pivots.pivot_colNZ[l]);
         assert(col.perm[c] < 0);
      }

      /*      remove pivot element from pivot row
       */
      setPivot(stage++, pcol, prow, u.row.val[i]);
      u.row.idx[i] = u.row.idx[j];
      u.row.val[i] = u.row.val[j];

      j = u.row.start[prow];
      for (--i; i >= j; --i)
      {
         c = u.row.idx[i];
         l = u.col.len[c] + u.col.start[c] - (temp.s_cact[c])--;
         for (k = l; u.col.idx[k] != prow; ++k)
            ;
         u.col.idx[k] = u.col.idx[l];
         u.col.idx[l] = prow;
         l = temp.s_cact[c];
         removeDR(pivots.pivot_col[c]);
         init2DR(pivots.pivot_col[c], pivots.pivot_colNZ[l]);
         assert(col.perm[c] < 0);
      }
   }

   initDR(pivots.pivot_colNZ[1]);           /* Remove all column singletons from list */
}


/*
    * No singletons available: Select pivot elements.
    */
void CLUFactor::selectPivots( double threshold, const int stage)
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

   candidates = thedim - stage - 1;
   if (candidates > 4)
      candidates = 4;

   num = 0;
   count = 2;

   for (;;)
   {
      ii = -1;

      if (pivots.pivot_rowNZ[count].next != &(pivots.pivot_rowNZ[count]))
      {
         rw = pivots.pivot_rowNZ[count].next->idx;
         beg = u.row.start[rw];
         len = u.row.len[rw] + beg - 1;

         /*  set maxabs to maximum absolute value in row
          *  (compute it if necessary).
          */
         if ((maxabs = temp.s_max[rw]) < 0)
         {
            maxabs = u.row.val[len];
            if (maxabs < 0)
               maxabs = -maxabs;
            for (i = len - 1; i >= beg; --i)
               if (maxabs < u.row.val[i])
                  maxabs = u.row.val[i];
               else if (maxabs < - u.row.val[i])
                  maxabs = - u.row.val[i];
            temp.s_max[rw] = maxabs;               /* ##### */
         }
         maxabs *= threshold;

         /*  select pivot element with lowest markowitz number in row
          */
         mkwtz = thedim + 1;
         for (i = len; i >= beg; --i)
         {
            k = u.row.idx[i];
            j = temp.s_cact[k];
            x = u.row.val[i];
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

      else if (pivots.pivot_colNZ[count].next != &(pivots.pivot_colNZ[count]))
      {
         cl = pivots.pivot_colNZ[count].next->idx;
         beg = u.col.start[cl];
         len = u.col.len[cl] + beg - 1;
         beg = len - temp.s_cact[cl] + 1;
         assert(count == temp.s_cact[cl]);

         /*  select pivot element with lowest markowitz number in column
          */
         mkwtz = thedim + 1;
         for (i = len; i >= beg; --i)
         {
            k = u.col.idx[i];
            j = u.row.len[k];
            if (j < mkwtz)
            {
               /*  ensure that element (cl,k) is stable.
                */
               if (temp.s_max[k] > 0)
               {
                  /*  case 1: maxabs is known
                   */
                  for (m = u.row.start[k], l = m + u.row.len[k] - 1; l >= m; --l)
                  {
                     if (u.row.idx[l] == cl)
                     {
                        x = u.row.val[l];
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
                  m = u.row.start[k];
                  maxabs = u.row.val[m];
                  maxabs = (maxabs > 0) ? maxabs : -maxabs;
                  for (l = m + u.row.len[k] - 1; l >= m; --l)
                  {
                     if (maxabs < u.row.val[l])
                        maxabs = u.row.val[l];
                     else if (maxabs < - u.row.val[l])
                        maxabs = - u.row.val[l];
                     if (u.row.idx[l] == cl)
                     {
                        x = u.row.val[l];
                        ll = l;
                        break;
                     }
                  }
                  for (--l; l > m; --l)
                  {
                     if (maxabs < u.row.val[l])
                        maxabs = u.row.val[l];
                     else if (maxabs < - u.row.val[l])
                        maxabs = - u.row.val[l];
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

      removeDR(pivots.pivot_col[cl]);
      initDR(pivots.pivot_col[cl]);

      if (ii >= 0)
      {
         /*  Initialize selected pivot element
          */
         CLUFactor::Pring *pr;
         pivots.pivot_row[rw].pos = ii - u.row.start[rw];
         pivots.pivot_row[rw].mkwtz = mkwtz = (mkwtz - 1) * (count - 1);  // ??? mkwtz originally was long, maybe to avoid an overflow in this instruction?
         for (pr = pivots.pivots.next; pr->idx >= 0; pr = pr->next)
         {
            if (pr->idx == rw || pr->mkwtz >= mkwtz)
               break;
         }
         pr = pr->prev;
         if (pr->idx != rw)
         {
            removeDR(pivots.pivot_row[rw]);
            init2DR (pivots.pivot_row[rw], *pr);
         }
         num++;
         if (num >= candidates)
            break;
      }
   }

   /*
     while(pivots.pivots.next->mkwtz < pivots.pivots.prev->mkwtz)
     {
     Pring   *pr;
     pr = pivots.pivots.prev;
     removeDR(*pr);
     init2DR (*pr, rowNZ[u.row.len[pr->idx]]);
     }
   */

   assert(row.perm[rw] < 0);
   assert(col.perm[cl] < 0);
}


/*
    *      Perform L and update loop for row r
    */
int CLUFactor::updateRow   (int r,
   int lv,
   int prow,
   int pcol,
   double pval,
   double eps )
{
   int fill;
   double x, lx;
   int c, i, j, k, ll, m, n;

   n = u.row.start[r];
   m = --(u.row.len[r]) + n;

   /*  compute L vector entry and
    *  and remove pivot column form row file
    */
   for (j = m; u.row.idx[j] != pcol; --j)
      ;
   lx = u.row.val[j] / pval;
   l.val[lv] = lx;
   l.idx[lv] = r;
   ++lv;

   u.row.idx[j] = u.row.idx[m];
   u.row.val[j] = u.row.val[m];


   /*  update loop (I) and
    *  computing expected fill
    */
   fill = u.row.len[prow];
   for (j = m - 1; j >= n; --j)
   {
      c = u.row.idx[j];
      if (temp.s_mark[c])
      {
         /*  count fill elements.
          */
         temp.s_mark[c] = 0;
         --fill;

         /*  update row values
          */
         x = u.row.val[j] -= work[c] * lx;
         if (isZero(x, eps))
         {
            /* Eliminate zero from row r
             */
            --u.row.len[r];
            --m;
            u.row.val[j] = u.row.val[m];
            u.row.idx[j] = u.row.idx[m];

            /* Eliminate zero from column c
             */
            --(temp.s_cact[c]);
            k = --(u.col.len[c]) + u.col.start[c];
            for (i = k; u.col.idx[i] != r; --i)
               ;
            u.col.idx[i] = u.col.idx[k];
         }
      }
   }


   /*  create space for fill in row file
    */
   ll = u.row.len[r];
   if (ll + fill > u.row.max[r])
      remaxRow(r, ll + fill);
   ll += u.row.start[r];

   /*  fill creating update loop (II)
    */
   for (j = u.row.start[prow], m = j + u.row.len[prow]; j < m; ++j)
   {
      c = u.row.idx[j];
      if (temp.s_mark[c])
      {
         x = - work[c] * lx;
         if (isNonZero(x, eps))
         {
            /* produce fill element in row r
             */
            u.row.val[ll] = x;
            u.row.idx[ll] = c;
            ll++;
            u.row.len[r]++;

            /* produce fill element in column c
             */
            if (u.col.len[c] >= u.col.max[c])
               remaxCol(c, u.col.len[c] + 1);
            u.col.idx[u.col.start[c] + (u.col.len[c])++] = r;
            temp.s_cact[c]++;
         }
      }
      else
         temp.s_mark[c] = 1;
   }

   /*  move row to appropriate list.
    */
   removeDR(pivots.pivot_row[r]);
   init2DR(pivots.pivot_row[r], pivots.pivot_rowNZ[u.row.len[r]]);
   assert(row.perm[r] < 0);
   temp.s_max[r] = -1;

   return lv;
}

/*
    *      Eliminate pivot element
    */
void CLUFactor::eliminatePivot(int prow, int pos, double eps, int& stage )
{
   int i, j, k, l = -1;
   int lv = -1;  // This value should never be used.
   int pcol;
   double pval;
   int pbeg = u.row.start[prow];
   int plen = --(u.row.len[prow]);
   int pend = pbeg + plen;


   /*  extract pivot element   */
   i = pbeg + pos;
   pcol = u.row.idx[i];
   pval = u.row.val[i];
   removeDR(pivots.pivot_col[pcol]);
   initDR(pivots.pivot_col[pcol]);

   /*  remove pivot from pivot row     */
   u.row.idx[i] = u.row.idx[pend];
   u.row.val[i] = u.row.val[pend];

   /*  set pivot element and construct L vector */
   setPivot(stage++, pcol, prow, pval);

   /**@todo If this test failes, lv has no value. I suppose that in this
    *       case none of the loops below that uses lv is executed.
    *       But this is unproven.
    */
   if (temp.s_cact[pcol] - 1 > 0)
      lv = makeLvec(temp.s_cact[pcol] - 1, prow);

   /*  init working vector,
    *  remove pivot row from working matrix
    *  and remove columns from list.
    */
   for (i = pbeg; i < pend; ++i)
   {
      j = u.row.idx[i];
      temp.s_mark[j] = 1;
      work[j] = u.row.val[i];
      removeDR(pivots.pivot_col[j]);
      l = u.col.start[j] + u.col.len[j] - temp.s_cact[j];
      for (k = l; u.col.idx[k] != prow; ++k)
         ;
      u.col.idx[k] = u.col.idx[l];
      u.col.idx[l] = prow;
      temp.s_cact[j]--;
   }

   /*  perform L and update loop
    */
   for
      (
         i = u.col.len[pcol] - temp.s_cact[pcol];
         (l = u.col.idx[u.col.start[pcol] + i]) != prow;
         ++i
         )
   {
      assert(row.perm[l] < 0);
      assert(lv >= 0);
      updateRow(l, lv++, prow, pcol, pval, eps);
   }

   /*  skip pivot row  */

   l = u.col.len[pcol];
   for (++i; i < l; ++i)
   {
      assert(lv >= 0);
      updateRow(u.col.idx[u.col.start[pcol] + i], lv++, prow, pcol, pval, eps);
   }

   /*  remove pivot column from column file.
    */
   u.col.len[pcol] -= temp.s_cact[pcol];

   /*  clear working vector and reinsert columns to lists
    */
   for (i = u.row.start[prow], pend = i + plen; i < pend; ++i)
   {
      j = u.row.idx[i];
      work[j] = 0;
      temp.s_mark[j] = 0;
      init2DR(pivots.pivot_col[j], pivots.pivot_colNZ[temp.s_cact[j]]);
      assert(col.perm[j] < 0);
   }
}


/*
    *      Factorize nucleus.
    */
void CLUFactor::eliminateNucleus( const double eps, 
   const double threshold, 
   int& stage)
{
   int r, c;
   CLUFactor::Pring *pivot;

   pivots.pivots.mkwtz = -1;
   pivots.pivots.idx = -1;
   pivots.pivots.pos = -1;

   while (stage < thedim - 1)
   {
#ifdef DEBUG
      int i;
      CLUFactorIsConsistent(fac);
      for (i = 0; i < thedim; ++i)
         if (col.perm[i] < 0)
         {
            assert(s_mark[i] == 0);
         }
#endif

      if (pivots.pivot_rowNZ[1].next != &(pivots.pivot_rowNZ[1]))        /* row singleton available */
         eliminateRowSingletons( stage );
      else if (pivots.pivot_colNZ[1].next != &(pivots.pivot_colNZ[1]))   /* column singleton available */
         eliminateColSingletons( stage );
      else
      {
         initDR(pivots.pivots);
         selectPivots( threshold, stage );

         assert ( pivots.pivots.next != &pivots.pivots &&  "ERROR: no pivot element selected" );

         for (pivot = pivots.pivots.next; pivot != &pivots.pivots; pivot = pivot->next)
         {
            eliminatePivot(pivot->idx, pivot->pos, eps, stage);
         }
      }

      if (pivots.pivot_rowNZ->next != pivots.pivot_rowNZ || pivots.pivot_colNZ->next != pivots.pivot_colNZ)
      {
         stat = SLinSolver::SINGULAR;
         return;
      }
   }

   if (stage < thedim)
   {
      /*      Eliminate remaining element.
       *      Note, that this must be both, column and row singleton.
       */
      assert(pivots.pivot_rowNZ[1].next != &(pivots.pivot_rowNZ[1]) && "ERROR: one row must be left");
      assert(pivots.pivot_colNZ[1].next != &(pivots.pivot_colNZ[1]) && "ERROR: one col must be left");
      r = pivots.pivot_rowNZ[1].next->idx;
      c = pivots.pivot_colNZ[1].next->idx;
      u.row.len[r] = 0;
      u.col.len[c]--;
      setPivot(stage, c, r, u.row.val[u.row.start[r]]);
   }
}

/*****************************************************************************/

int CLUFactor::setupColVals( CLUFactor* fc )
{
   int i, j, k, n;
   int* idx;
   double* val;
   double* cval;
   double maxabs;

   if (u.col.val)
      spx_free(u.col.val);

   spx_alloc(u.col.val, u.col.size);

   cval = u.col.val;

   for (i = fc->thedim - 1; i >= 0; --i)
      fc->u.col.len[i] = 0;

   maxabs = 0;
   for (n = fc->thedim, i = fc->thedim - 1; i >= 0; --i)
   {
      j = fc->u.row.start[i];
      idx = & fc->u.row.idx[j];
      val = & fc->u.row.val[j];
      j = fc->u.row.len[i];
      n += j;
      while (j--)
      {
         k = fc->u.col.start[*idx] + fc->u.col.len[*idx]++;
         fc->u.col.idx[k] = i;
         cval[k] = *val;
         if (*val > maxabs)
            maxabs = *val;
         else if (-*val > maxabs)
            maxabs = -*val;
         ++idx;
         ++val;
      }
   }
   fc->maxabs = maxabs;
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

 
void CLUFactor::factor( 
   SVector** vec,   /* Array of column vector pointers   */
   double threshold,           /* pivoting threshold                */
   double eps           /* epsilon for zero detection        */
                   )
{
#warning remove this assignmnent
   //fac = this;
   int stage = 0;
   stat = SLinSolver::OK;

   l.start[0] = 0;
   l.firstUpdate = 0;
   l.firstUnused = 0;

   temp.init(thedim);
   initPerm();

   initFactorMatrix(vec, eps, stage);
   if (stat)
      goto TERMINATE;
   //   initMaxabs = initMaxabs;

   colSingletons(stage);
   if (stat != SLinSolver::OK)
      goto TERMINATE;

   rowSingletons(stage);
   if (stat != SLinSolver::OK)
      goto TERMINATE;

   if (stage < thedim)
   {
      initFactorRings(stage);
      eliminateNucleus(eps,threshold,stage);
      freeFactorRings();
   }

 TERMINATE:
   l.firstUpdate = l.firstUnused;

   if (stat == SLinSolver::OK)
   {
#ifdef WITH_L_ROWS
      setupRowVals(this);
#endif
      nzCnt = setupColVals(this);
   }

   /* assert(dump()); */
}

void CLUFactor::dump() const
{
   int i, j, k;

   /*  Dump U:
    */
   for (i = 0; i < thedim; ++i)
   {
      if (row.perm[i] >= 0)
         std::cout << "diag[" << i << "]: [" << col.orig[row.perm[i]] 
                   << "] = " << diag[i] << std::endl;

      for (j = 0; j < u.row.len[i]; ++j)
         std::cout << "   u[" << i << "]: [" 
                   << u.row.idx[u.row.start[i] + j] << "] = "
                   << u.row.val[u.row.start[i] + j] << std::endl;
   }

   /*  Dump L:
    */
   for (i = 0; i < thedim; ++i)
   {
      for (j = 0; j < l.firstUnused; ++j)
         if (col.orig[row.perm[l.row[j]]] == i)
         {
            std::cout << "l[" << i << "]" << std::endl;

            for (k = l.start[j]; k < l.start[j + 1]; ++k)
               std::cout << "   l[" << k - l.start[j]
                         << "]:  [" << l.idx[k]
                         << "] = " << l.val[k] << std::endl;
            break;
         }
   }
   return;
}

/*****************************************************************************/
/*
 *      Ensure that row memory is at least size.
    */
void CLUFactor::minRowMem(int size)
{
   if (u.row.size < size)
   {
      u.row.size = size;
      spx_realloc(u.row.val, size);
      spx_realloc(u.row.idx, size);
   }
}

/*****************************************************************************/
/*
 *      Ensure that column memory is at least size.
    */
void CLUFactor::minColMem(int size)
{
   if (u.col.size < size)
   {
      u.col.size = size;
      spx_realloc(u.col.idx, size);
   }
}

void CLUFactor::minLMem(int size)
{
   if (size > l.size)
   {
      l.size = int(0.2 * l.size + size);
      spx_realloc(l.val, l.size);
      spx_realloc(l.idx, l.size);
   }
}

int CLUFactor::makeLvec(int p_len, int p_row)
{
   int* p_lrow = l.row;
   int* p_lbeg = l.start;
   int first   = p_lbeg[l.firstUnused];

   if (l.firstUnused >= l.startSize)
   {
      l.startSize += 100;
      spx_realloc(l.start, l.startSize);
      p_lbeg = l.start;
   }

   assert(p_len > 0 && "ERROR: no empty columns allowed in L vectors");

   minLMem(first + p_len);
   p_lrow[l.firstUnused] = p_row;
   l.start[++(l.firstUnused)] = first + p_len;

   assert(l.start[l.firstUnused] <= l.size);
   assert(l.firstUnused <= l.startSize);
   return first;
}



/*****************************************************************************/

bool CLUFactor::isConsistent() const
{
   int              i, j, k, ll;
   Dring            *ring;
   CLUFactor::Pring *pring;

   /*  Consistency only relevant for real factorizations
    */
   if (stat)
      return true;

   /*  Test column ring list consistency.
    */
   i = 0;
   for (ring = u.col.list.next; ring != &(u.col.list); ring = ring->next)
   {
      assert(ring->idx >= 0);
      assert(ring->idx < thedim);
      assert(ring->prev->next == ring);
      if (ring != u.col.list.next)
      {
         assert( u.col.start[ring->prev->idx] + u.col.max[ring->prev->idx]
            == u.col.start[ring->idx] );
      }
      ++i;
   }
   assert(i == thedim);
   assert(u.col.start[ring->prev->idx] + u.col.max[ring->prev->idx]
      == u.col.used);


   /*  Test row ring list consistency.
    */
   i = 0;
   for (ring = u.row.list.next; ring != &(u.row.list); ring = ring->next)
   {
      assert(ring->idx >= 0);
      assert(ring->idx < thedim);
      assert(ring->prev->next == ring);
      assert(u.row.start[ring->prev->idx] + u.row.max[ring->prev->idx]
         == u.row.start[ring->idx]);
      ++i;
   }
   assert(i == thedim);
   assert(u.row.start[ring->prev->idx] + u.row.max[ring->prev->idx]
      == u.row.used);


   /*  Test consistency of individual svectors.
    */
   for (i = 0; i < thedim; ++i)
   {
      assert(u.row.max[i] >= u.row.len[i]);
      assert(u.col.max[i] >= u.col.len[i]);
   }


   /*  Test consistency of column file to row file of U
    */
   for (i = 0; i < thedim; ++i)
   {
      for(j = u.row.start[i] + u.row.len[i] - 1;
          j >= u.row.start[i];
          j--)
      {
         k = u.row.idx[j];
         for(ll = u.col.start[k] + u.col.len[k] - 1;
             ll >= u.col.start[k];
             ll-- )
         {
            if (u.col.idx[ll] == i)
               break;
         }
         assert(u.col.idx[ll] == i);
         if (row.perm[i] < 0)
         {
            assert(col.perm[k] < 0);
         }
         else
         {
            assert(col.perm[k] < 0 
               || col.perm[k] > row.perm[i]);
         }
      }
   }

   /*  Test consistency of row file to column file of U
    */
   for (i = 0; i < thedim; ++i)
   {
      for(j = u.col.start[i] + u.col.len[i] - 1;
          j >= u.col.start[i];
          j--)
      {
         k = u.col.idx[j];
         for( ll = u.row.start[k] + u.row.len[k] - 1;
              ll >= u.row.start[k];
              ll--)
         {
            if (u.row.idx[ll] == i)
               break;
         }
         assert(u.row.idx[ll] == i);
         assert(col.perm[i] < 0
            || row.perm[k] < col.perm[i]);
      }
   }

   //       impossible without global stage variable
   //       /*  Test consistency of nonzero count lists
   //        */
   //       if (colNZ && rowNZ)
   //          for (i = 0; i < thedim - stage; ++i)
   //             {
   //                for (pring = rowNZ[i].next; pring != &(rowNZ[i]); pring = pring->next)
   //                   {
   //                      assert(row.perm[pring->idx] < 0);
   //                   }
   //                for (pring = colNZ[i].next; pring != &(colNZ[i]); pring = pring->next)
   //                   {
   //                      assert(col.perm[pring->idx] < 0);
   //                   }
   //             }

   return true;
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
