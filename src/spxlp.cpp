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
#pragma ident "@(#) $Id: spxlp.cpp,v 1.7 2001/12/25 16:03:24 bzfkocht Exp $"

#include <stdio.h>

#include "spxlp.h"
#include "spxmessage.h"

namespace soplex
{
const double SPxLP::infinity = 1e+100;

void SPxLP::getRow(int i, LPRow& row) const
{
   row.lhs() = lhs(i);
   row.rhs() = rhs(i);
   row.rowVector() = rowVector(i);
}

void SPxLP::getRows(int start, int end, LPRowSet& p_set) const
{
   int i;
   p_set.clear();
   for (i = 0; start <= end; ++i, ++start)
      p_set.add(lhs(start), rowVector(start), rhs(start));
}

void SPxLP::getCol(int i, LPCol& col) const
{
   col.upper() = upper(i);
   col.lower() = lower(i);
   col.obj() = spxSense() * obj(i);
   col.colVector() = colVector(i);
}

void SPxLP::getCols(int start, int end, LPColSet& p_set) const
{
   int i;
   p_set.clear();
   for (i = 0; start <= end; ++i, ++start)
      p_set.add(obj(start), lower(start), colVector(start), upper(start));
}

void SPxLP::getObj(Vector& p_obj) const
{
   p_obj = LPColSet::obj();
   if (spxSense() == MINIMIZE)
      p_obj *= -1;
}

void SPxLP::doAddRow(const LPRow& row)
{
   int idx = nRows();
   const SVector& vec = row.rowVector();

   LPRowSet::add(row);

   // now insert nonzeros to column file also
   for (int j = vec.size() - 1; j >= 0; --j)
   {
      double val = vec.value(j);
      int i = vec.index(j);
      if (i >= nCols())      // create new columns if required
      {
         LPCol empty;
         for (int k = nCols(); k <= i; ++k)
            doAddCol(empty);
         //              LPColSet::add(empty);
      }

      assert(i < nCols());
      LPColSet::add2(i, 1, &idx, &val);
   }

   addedRows(1);
}

void SPxLP::doAddCol(const LPCol& col)
{
   int idx = nCols();
   const SVector& vec = col.colVector();

   LPColSet::add(col);
   LPColSet::obj(idx) *= thesense;

   // now insert nonzeros to row file also
   for (int j = vec.size() - 1; j >= 0; --j)
   {
      double val = vec.value(j);
      int i = vec.index(j);
      if (i >= nRows())              // create new rows if required
      {
         LPRow empty;
         for (int k = nRows(); k <= i; ++k)
            doAddRow(empty);
         //              LPRowSet::add(empty);
      }

      assert(i < nRows());
      LPRowSet::add2(i, 1, &idx, &val);
   }

   addedCols(1);
}

void SPxLP::added2Set(SVSet& p_set, const SVSet& p_add, int n)
{
   int i, j, end, tot;
   DataArray < int > moreArray(p_set.num());
   int* more = moreArray.get_ptr();

   for (i = p_set.num() - 1; i >= 0; --i)
      more[i] = 0;

   for (tot = 0, i = p_add.num() - n, end = p_add.num(); i < end; ++i)
   {
      const SVector& vec = p_add[i];
      tot += vec.size();
      for (j = vec.size() - 1; j >= 0; --j)
         more[ vec.index(j) ]++;
   }

   if (p_set.memMax() < tot)
      p_set.memRemax(tot);

   for (i = p_set.num() - 1; i >= 0; --i)
   {
      j = p_set[i].size();
      p_set.xtend(p_set[i], j + more[i]);
      p_set[i].set_size( j + more[i] );
      more[i] = j;
   }

   for (i = p_add.num() - n; i < p_add.num(); ++i)
   {
      const SVector& vec = p_add[i];
      for (j = vec.size() - 1; j >= 0; --j)
      {
         int n = vec.index(j);
         int m = more[n] ++;
         SVector& l_xtend = p_set[n];
         l_xtend.index(m) = i;
         l_xtend.value(m) = vec.value(j);
      }
   }
}

void SPxLP::doAddRows(const LPRowSet& p_set)
{
   int i, j, k, ii, idx;
   SVector* col;
   DataArray < int > newCols(nCols());
   DataArray < double > newColVals(nCols());
   int oldRowNumber = nRows();

   if (&p_set != this)
      LPRowSet::add(p_set);
   assert(LPRowSet::isConsistent());
   assert(LPColSet::isConsistent());

   // count additional nonzeros per column
   for (i = nCols() - 1; i >= 0; --i)
   {
      newCols[i] = 0;
      newColVals[i] = 0;
   }
   for (i = p_set.num() - 1; i >= 0; --i)
   {
      const SVector& vec = p_set.rowVector(i);
      for (j = vec.size() - 1; j >= 0; --j)
      {
         ii = vec.index(j);
         if (ii >= nCols()) // create new columns if required
         {
            LPCol empty;
            newCols.reSize(ii + 1);
            newColVals.reSize(ii + 1);
            for (k = nCols(); k <= ii; ++k)
            {
               newCols[k] = 0;
               newColVals[k] = 0;
               doAddCol(empty);
               //                  LPColSet::add(empty);
            }

         }
         assert(ii < nCols());
         newCols[ii]++;
      }
   }

   // extend columns as required, with dummy elements
   for (i = newCols.size() - 1; i >= 0; --i)
   {
      if (newCols[i] >= 0)
         LPColSet::add2(i, newCols[i],
                        newCols.get_ptr(), newColVals.get_ptr());
   }

   // insert new elements to column file
   for (i = nRows() - 1; i >= oldRowNumber; --i)
   {
      const SVector& vec = rowVector(i);
      for (j = vec.size() - 1; j >= 0; --j)
      {
         k = vec.index(j);
         col = &colVector_w(k);
         idx = col->size() - newCols[k];
         assert(newCols[k] >= 0);
         newCols[k]--;
         col->index(idx) = i;
         col->value(idx) = vec.value(j);
      }
   }

   addedRows(p_set.num());
}

void SPxLP::doAddCols(const LPColSet& p_set)
{
   int i, j;
   int oldColNumber = nCols();
   DataArray < int > newRows(nRows());
   DataArray < double > newRowVals(nRows());

   if (&p_set != this)
      LPColSet::add(p_set);
   assert(LPColSet::isConsistent());

   // count additional nonzeros per row
   for (i = nRows() - 1; i >= 0; --i)
   {
      newRows[i] = 0;
      newRowVals[i] = 0;
   }
   for (i = p_set.num() - 1; i >= 0; --i)
   {
      const SVector& vec = p_set.colVector(i);
      for (j = vec.size() - 1; j >= 0; --j)
      {
         int l = vec.index(j);
         if (l >= nRows())  // create new rows if required
         {
            LPRow empty;
            newRows.reSize(l + 1);
            newRowVals.reSize(l + 1);
            for (int k = nRows(); k <= l; ++k)
            {
               newRows[k] = 0;
               newRowVals[k] = 0;
               doAddRow(empty);
               //                  LPRowSet::add(empty);
            }

         }
         assert(l < nRows());
         newRows[l]++;
      }
   }

   // extend rows as required, with dummy elements
   for (i = newRows.size() - 1; i >= 0; --i)
   {
      if (newRows[i] >= 0)
      {
         int len = newRows[i] + rowVector(i).size();
         LPRowSet::xtend(i, len);
         rowVector_w(i).set_size( len );
      }
   }

   // insert new elements to row file
   for (i = oldColNumber; i < nCols(); ++i)
   {
      LPColSet::obj(i) *= thesense;
      const SVector& vec = colVector(i);
      for (j = vec.size() - 1; j >= 0; --j)
      {
         int k = vec.index(j);
         SVector& col = rowVector_w(k);
         int idx = col.size() - newRows[k];
         assert(newRows[k] >= 0);
         newRows[k]--;
         col.index(idx) = i;
         col.value(idx) = vec.value(j);
      }
   }
   assert(SPxLP::isConsistent());

   addedCols(p_set.num());
}

void SPxLP::addRows(SPxRowId id[], const LPRowSet& p_set)
{
   int i = nRows();
   addRows(p_set);
   for (int j = 0; i < nRows(); ++i, ++j)
      id[j] = rId(i);
}

void SPxLP::addCols(SPxColId id[], const LPColSet& p_set)
{
   int i = nCols();
   addCols(p_set);
   for (int j = 0; i < nCols(); ++i, ++j)
      id[j] = cId(i);
}

void SPxLP::doRemoveRow(int j)
{
   const SVector& vec = rowVector(j);
   int i;

   // remove row vector from column file
   for (i = vec.size() - 1; i >= 0; --i)
   {
      SVector& remvec = colVector_w(vec.index(i));
      remvec.remove(remvec.number(j));
   }

   int idx = nRows() - 1;
   if (j != idx)              // move last row to removed position
   {
      const SVector& l_vec = rowVector(idx);
      for (i = l_vec.size() - 1; i >= 0; --i)
      {
         SVector& movevec = colVector_w(l_vec.index(i));
         movevec.index(movevec.number(idx)) = j;
      }
   }

   LPRowSet::remove(j);
}

void SPxLP::doRemoveCol(int j)
{
   const SVector& vec = colVector(j);
   int i;

   // remove column vector from row file
   for (i = vec.size() - 1; i >= 0; --i)
   {
      SVector& remvec = rowVector_w(vec.index(i));
      remvec.remove(remvec.number(j));
   }

   int idx = nCols() - 1;
   if (j != idx)              // move last column to removed position
   {
      const SVector& l_vec = colVector(idx);
      for (i = l_vec.size() - 1; i >= 0; --i)
      {
         SVector& movevec = rowVector_w(l_vec.index(i));
         movevec.index(movevec.number(idx)) = j;
      }
   }

   LPColSet::remove(j);
}

void SPxLP::doRemoveRows(int perm[])
{
   int j = nCols();

   LPRowSet::remove(perm);
   for (int i = 0; i < j; ++i)
   {
      SVector& vec = colVector_w(i);
      for (int k = vec.size() - 1; k >= 0; --k)
      {
         int idx = vec.index(k);
         if (perm[idx] < 0)
            vec.remove(k);
         else
            vec.index(k) = perm[idx];
      }
   }
}

void SPxLP::doRemoveCols(int perm[])
{
   int j = nRows();

   LPColSet::remove(perm);
   for (int i = 0; i < j; ++i)
   {
      SVector& vec = rowVector_w(i);
      for (int k = vec.size() - 1; k >= 0; --k)
      {
         int idx = vec.index(k);
         if (perm[idx] < 0)
            vec.remove(k);
         else
            vec.index(k) = perm[idx];
      }
   }
}

void SPxLP::removeRows(SPxRowId id[], int n, int perm[])
{
   if (perm == 0)
   {
      DataArray < int > p(nRows());
      removeRows(id, n, p.get_ptr());
      return;
   }
   for (int i = nRows() - 1; i >= 0; --i)
      perm[i] = i;
   while (n--)
      perm[number(id[n])] = -1;
   removeRows(perm);
}

void SPxLP::removeRows(int nums[], int n, int perm[])
{
   if (perm == 0)
   {
      DataArray < int > p(nRows());
      removeRows(nums, n, p.get_ptr());
      return;
   }
   for (int i = nRows() - 1; i >= 0; --i)
      perm[i] = i;
   while (n--)
      perm[nums[n]] = -1;
   removeRows(perm);
}

void SPxLP::removeRowRange(int start, int end, int* perm)
{
   if (perm == 0)
   {
      int i = end - start + 1;
      DataArray < int > p(i);
      while (--i >= 0)
         p[i] = start + i;
      removeRows(p.get_ptr(), end - start + 1);
      return;
   }
   int i;
   for (i = 0; i < start; ++i)
      perm[i] = i;
   for (; i <= end; ++i)
      perm[i] = -1;
   for (; i < nRows(); ++i)
      perm[i] = i;
   removeRows(perm);
}

void SPxLP::removeCols(SPxColId id[], int n, int perm[])
{
   if (perm == 0)
   {
      DataArray < int > p(nCols());
      removeCols(id, n, p.get_ptr());
      return;
   }
   for (int i = nCols() - 1; i >= 0; --i)
      perm[i] = i;
   while (n--)
      perm[number(id[n])] = -1;
   removeCols(perm);
}

void SPxLP::removeCols(int nums[], int n, int perm[])
{
   if (perm == 0)
   {
      DataArray < int > p(nCols());
      removeCols(nums, n, p.get_ptr());
      return;
   }
   for (int i = nCols() - 1; i >= 0; --i)
      perm[i] = i;
   while (n--)
      perm[nums[n]] = -1;
   removeCols(perm);
}

void SPxLP::removeColRange(int start, int end, int* perm)
{
   if (perm == 0)
   {
      int i = end - start + 1;
      DataArray < int > p(i);
      while (--i >= 0)
         p[i] = start + i;
      removeCols(p.get_ptr(), end - start + 1);
      return;
   }
   int i;
   for (i = 0; i < start; ++i)
      perm[i] = i;
   for (; i <= end; ++i)
      perm[i] = -1;
   for (; i < nCols(); ++i)
      perm[i] = i;
   removeCols(perm);
}

void SPxLP::clear()
{
   LPRowSet::clear();
   LPColSet::clear();
   thesense = MAXIMIZE;
}

void SPxLP::changeObj(const Vector& newObj)
{
   assert(maxObj().dim() == newObj.dim());
   LPColSet::obj() = newObj;
   LPColSet::obj() *= spxSense();
   assert(isConsistent());
}

void SPxLP::changeObj(int i, double newVal)
{
   LPColSet::obj(i) = spxSense() * newVal;
   assert(isConsistent());
}

void SPxLP::changeLower(const Vector& newLower)
{
   assert(lower().dim() == newLower.dim());
   LPColSet::lower() = newLower;
   assert(isConsistent());
}

void SPxLP::changeLower(int i, double newLower)
{
   LPColSet::lower(i) = newLower;
   assert(isConsistent());
}

void SPxLP::changeUpper(const Vector& newUpper)
{
   assert(upper().dim() == newUpper.dim());
   LPColSet::upper() = newUpper;
   assert(isConsistent());
}

void SPxLP::changeUpper(int i, double newUpper)
{
   LPColSet::upper(i) = newUpper;
   assert(isConsistent());
}

void SPxLP::changeLhs(const Vector& newLhs)
{
   assert(lhs().dim() == newLhs.dim());
   LPRowSet::lhs() = newLhs;
   assert(isConsistent());
}

void SPxLP::changeBounds(const Vector& newLower, const Vector& newUpper)
{
   changeLower(newLower);
   changeUpper(newUpper);
   assert(isConsistent());
}

void SPxLP::changeBounds(int i, double newLower, double newUpper)
{
   changeLower(i, newLower);
   changeUpper(i, newUpper);
   assert(isConsistent());
}

void SPxLP::changeLhs(int i, double newLhs)
{
   LPRowSet::lhs(i) = newLhs;
   assert(isConsistent());
}

void SPxLP::changeRhs(const Vector& newRhs)
{
   assert(rhs().dim() == newRhs.dim());
   LPRowSet::rhs() = newRhs;
   assert(isConsistent());
}

void SPxLP::changeRhs(int i, double newRhs)
{
   LPRowSet::rhs(i) = newRhs;
   assert(isConsistent());
}

void SPxLP::changeRange(const Vector& newLhs, const Vector& newRhs)
{
   changeLhs(newLhs);
   changeRhs(newRhs);
   assert(isConsistent());
}

void SPxLP::changeRange(int i, double newLhs, double newRhs)
{
   changeLhs(i, newLhs);
   changeRhs(i, newRhs);
   assert(isConsistent());
}


void SPxLP::changeRow(int n, const LPRow& newRow)
{
   int j;
   SVector& row = rowVector_w(n);
   for (j = row.size() - 1; j >= 0; --j)
   {
      SVector& col = colVector_w(row.index(j));
      col.remove(col.number(n));
   }
   row.clear();

   changeLhs(n, newRow.lhs());
   changeRhs(n, newRow.rhs());
   const SVector& newrow = newRow.rowVector();
   for (j = newrow.size() - 1; j >= 0; --j)
   {
      int idx = newrow.index(j);
      double val = newrow.value(j);
      LPRowSet::add2(n, 1, &idx, &val);
      LPColSet::add2(idx, 1, &n, &val);
   }
   assert(isConsistent());
}

void SPxLP::changeCol(int n, const LPCol& newCol)
{
   int j;
   SVector& col = colVector_w(n);
   for (j = col.size() - 1; j >= 0; --j)
   {
      SVector& row = rowVector_w(col.index(j));
      row.remove(row.number(n));
   }
   col.clear();

   changeUpper(n, newCol.upper());
   changeLower(n, newCol.lower());
   changeObj (n, newCol.obj());
   const SVector& newcol = newCol.colVector();
   for (j = newcol.size() - 1; j >= 0; --j)
   {
      int idx = newcol.index(j);
      double val = newcol.value(j);
      LPColSet::add2(n, 1, &idx, &val);
      LPRowSet::add2(idx, 1, &n, &val);
   }
   assert(isConsistent());
}

void SPxLP::changeElement(int i, int j, double val)
{
   SVector& row = rowVector_w(i);
   SVector& col = colVector_w(j);

   if (val != 0)
   {
      if (row.number(j) >= 0)
      {
         row.value(row.number(j)) = val;
         col.value(col.number(i)) = val;
      }
      else
      {
         LPRowSet::add2(i, 1, &j, &val);
         LPColSet::add2(j, 1, &i, &val);
      }
   }
   else if (row.number(j) >= 0)
   {
      row.remove(row.number(j));
      col.remove(col.number(i));
   }
   assert(isConsistent());
}

int SPxLP::isConsistent() const
{
   int i, j, n;

   for (i = nCols() - 1; i >= 0; --i)
   {
      const SVector& v = colVector(i);
      for (j = v.size() - 1; j >= 0; --j)
      {
         const SVector& w = rowVector(v.index(j));
         n = w.number(i);
         if (n < 0)
            return SPXinconsistent("SPxLP");
         if (v.value(j) != w.value(n))
            return SPXinconsistent("SPxLP");
      }
   }

   for (i = nRows() - 1; i >= 0; --i)
   {
      const SVector& v = rowVector(i);
      for (j = v.size() - 1; j >= 0; --j)
      {
         const SVector& w = colVector(v.index(j));
         n = w.number(i);
         if (n < 0)
            return SPXinconsistent("SPxLP");
         if (v.value(j) != w.value(n))
            return SPXinconsistent("SPxLP");
      }
   }
   return LPRowSet::isConsistent() && LPColSet::isConsistent();
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
