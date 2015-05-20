/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  spxlpbase.h
 * @brief Saving LPs in a form suitable for SoPlex.
 */
#ifndef _SPXLPBASE_H_
#define _SPXLPBASE_H_

/* undefine SOPLEX_DEBUG flag from including files; if SOPLEX_DEBUG should be defined in this file, do so below */
#ifdef SOPLEX_DEBUG
#define SOPLEX_DEBUG_SPXLPBASE
#undef SOPLEX_DEBUG
#endif

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <typeinfo>

#include "spxdefines.h"
#include "basevectors.h"
#include "dataarray.h"
#include "datakey.h"
#include "spxid.h"
#include "lprowbase.h"
#include "lpcolbase.h"
#include "lprowsetbase.h"
#include "lpcolsetbase.h"
#include "nameset.h"
#include "didxset.h"
#include "spxfileio.h"

namespace soplex
{
class SPxSolver;

/**@brief   Saving LPs in a form suitable for SoPlex.
 * @ingroup Algo
 *
 *  Class SPxLPBase provides the data structures required for saving a linear program in the form
 *  \f[
 *  \begin{array}{rl}
 *      \hbox{max}  & c^T x              \\
 *      \hbox{s.t.} & l_r \le Ax \le u_r \\
 *                  & l_c \le x \le u_c
 *  \end{array}
 *  \f]
 *  suitable for solving with SoPlex. This includes:
 *  - SVSetBase%s for both columns and rows
 *  - objective Vector
 *  - upper and lower bound Vectors for variables (\f$l_c\f$ and \f$u_c\f$)
 *  - upper and lower bound Vectors for inequalities (\f$l_r\f$ and \f$u_r\f$)
 *
 *  Note, that the optimization sense is not saved directly. Instead, the objective function are multiplied by -1 to
 *  transform the LP to our standard form maximizing the objective function. However, the sense of the loaded LP can be
 *  retreived with method #spxSense().
 *
 *  Further, equality constraints are modelled by \f$l_r = u_r\f$.  Analogously, fixed variables have \f$l_c = u_c\f$.
 *
 *  #SPxLPBase%s are saved as an SVSet, both for columns and rows. Note that this is redundant but eases the access.
 */
template < class R >
class SPxLPBase : protected LPRowSetBase<R>, protected LPColSetBase<R>
{
   template < class S > friend class SPxLPBase;
   friend class SPxBasis;
   friend class SPxScaler;
   friend class SPxEquiliSC;
   friend class SPxGeometSC;
   friend class SPxMainSM;

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Types */
   //@{

   /// Optimization sense.
   enum SPxSense
   {
      MAXIMIZE = 1,
      MINIMIZE = -1
   };

   //@}

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   SPxSense thesense;   ///< optimization sense.
   R offset;            ///< offset computed, e.g., in simplification step

   //@}

public:

   // message handler
   SPxOut* spxout;

public:

   void setOutstream(SPxOut& newOutstream)
   {
      spxout = &newOutstream;
   }
   // ------------------------------------------------------------------------------------------------------------------
   /**@name Inquiry */
   //@{

   /// Returns number of rows in LP.
   int nRows() const
   {
      return LPRowSetBase<R>::num();
   }

   /// Returns number of columns in LP.
   int nCols() const
   {
      return LPColSetBase<R>::num();
   }

   /// Returns number of nonzeros in LP.
   int nNzos() const
   {

      int n = 0;
      for( int i = 0; i < nCols(); ++i )
         n += colVector(i).size();

      return n;
   }

   /// Absolute smallest non-zero element in LP.
   R minAbsNzo() const
   {

      R mini = infinity;

      for( int i = 0; i < nCols(); ++i )
      {
         R m = colVector(i).minAbs();

         if( m < mini )
            mini = m;
      }

      assert(mini >= R(0));

      return mini;
   }

   /// Absolute biggest non-zero element in LP.
   R maxAbsNzo() const
   {

      R maxi = R(0);

      for( int i = 0; i < nCols(); ++i )
      {
         R m = colVector(i).maxAbs();

         if( m > maxi )
            maxi = m;
      }

      assert(maxi >= R(0));

      return maxi;
   }

   /// Gets \p i 'th row.
   void getRow(int i, LPRowBase<R>& row) const
   {
      row.setLhs(lhs(i));
      row.setRhs(rhs(i));
      row.setObj(rowObj(i));
      row.setRowVector(DSVectorBase<R>(rowVector(i)));
   }

   /// Gets row with identifier \p id.
   void getRow(const SPxRowId& id, LPRowBase<R>& row) const
   {
      getRow(number(id), row);
   }

   /// Gets rows \p start, ... \p end.
   void getRows(int start, int end, LPRowSetBase<R>& set) const
   {

      set.clear();
      for( int i = start; i <= end; i++ )
         set.add(lhs(i), rowVector(i), rhs(i), rowObj(i));
   }

   /// Gets row vector of row \p i.
   const SVectorBase<R>& rowVector(int i) const
   {
      return LPRowSetBase<R>::rowVector(i);
   }

   /// Gets row vector of row with identifier \p id.
   const SVectorBase<R>& rowVector(const SPxRowId& id) const
   {
      return LPRowSetBase<R>::rowVector(id);
   }

   /// Returns right hand side vector.
   const VectorBase<R>& rhs() const
   {
      return LPRowSetBase<R>::rhs();
   }

   ///
   const R& rhs(int i) const
   {
      return LPRowSetBase<R>::rhs(i);
   }

   /// Returns right hand side of row with identifier \p id.
   const R& rhs(const SPxRowId& id) const
   {
      return LPRowSetBase<R>::rhs(id);
   }

   /// Returns left hand side vector.
   const VectorBase<R>& lhs() const
   {
      return LPRowSetBase<R>::lhs();
   }

   ///
   const R& lhs(int i) const
   {
      return LPRowSetBase<R>::lhs(i);
   }

   /// Returns left hand side of row with identifier \p id.
   const R& lhs(const SPxRowId& id) const
   {
      return LPRowSetBase<R>::lhs(id);
   }

   /// Gets row objective function vector.
   void getRowObj(VectorBase<R>& prowobj) const
   {
      prowobj = LPRowSetBase<R>::obj();
      if( spxSense() == MINIMIZE )
         prowobj *= -1.0;
   }

   ///
   R rowObj(int i) const
   {
      if( spxSense() == MINIMIZE )
         return -maxRowObj(i);
      else
         return maxRowObj(i);
   }

   /// Returns row objective function value of row with identifier \p id.
   R rowObj(const SPxRowId& id) const
   {
      if( spxSense() == MINIMIZE )
         return -maxRowObj(id);
      else
         return maxRowObj(id);
   }

   ///
   const VectorBase<R>& maxRowObj() const
   {
      return LPRowSetBase<R>::obj();
   }

   ///
   const R& maxRowObj(int i) const
   {
      return LPRowSetBase<R>::obj(i);
   }

   /// Returns row objective function value of row with identifier \p id.
   const R& maxRowObj(const SPxRowId& id) const
   {
      return LPRowSetBase<R>::obj(id);
   }

   /// Returns the inequality type of the \p i'th LPRow.
   typename LPRowBase<R>::Type rowType(int i) const
   {
      return LPRowSetBase<R>::type(i);
   }

   /// Returns the inequality type of the row with identifier \p key.
   typename LPRowBase<R>::Type rowType(const SPxRowId& id) const
   {
      return LPRowSetBase<R>::type(id);
   }

   /// Gets \p i 'th column.
   void getCol(int i, LPColBase<R>& col) const
   {

      col.setUpper(upper(i));
      col.setLower(lower(i));
      col.setObj(obj(i));
      col.setColVector(colVector(i));
   }

   /// Gets column with identifier \p id.
   void getCol(const SPxColId& id, LPColBase<R>& col) const
   {
      getCol(number(id), col);
   }

   /// Gets columns \p start, ..., \p end.
   void getCols(int start, int end, LPColSetBase<R>& set) const
   {

      set.clear();
      for( int i = start; i <= end; i++ )
         set.add(obj(i), lower(i), colVector(i), upper(i));
   }

   /// Returns column vector of column \p i.
   const SVectorBase<R>& colVector(int i) const
   {
      return LPColSetBase<R>::colVector(i);
   }

   /// Returns column vector of column with identifier \p id.
   const SVectorBase<R>& colVector(const SPxColId& id) const
   {
      return LPColSetBase<R>::colVector(id);
   }

   /// Gets objective vector.
   void getObj(VectorBase<R>& pobj) const
   {

      pobj = LPColSetBase<R>::maxObj();
      if( spxSense() == MINIMIZE )
         pobj *= -1.0;
   }

   /// Returns objective value of column \p i.
   R obj(int i) const
   {
      R res = maxObj(i);
      if( spxSense() == MINIMIZE )
         res *= -1;
      return res;
   }

   /// Returns objective value of column with identifier \p id.
   R obj(const SPxColId& id) const
   {
      R res = maxObj(id);
      if( spxSense() == MINIMIZE )
         res *= -1;
      return res;
   }

   /// Returns objective vector for maximization problem.
   /** Methods #maxObj() return the objective vector or its elements, after transformation to a maximization
    *  problem. Since this is how SPxLPBase internally stores any LP these methods are generally faster. The following
    *  condition holds: #obj() = #spxSense() * maxObj().
    */
   const VectorBase<R>& maxObj() const
   {
      return LPColSetBase<R>::maxObj();
   }

   /// Returns objective value of column \p i for maximization problem.
   const R& maxObj(int i) const
   {
      return LPColSetBase<R>::maxObj(i);
   }

   /// Returns objective value of column with identifier \p id for maximization problem.
   const R& maxObj(const SPxColId& id) const
   {
      return LPColSetBase<R>::maxObj(id);
   }

   /// Returns upper bound vector.
   const VectorBase<R>& upper() const
   {
      return LPColSetBase<R>::upper();
   }

   /// Returns upper bound of column \p i.
   const R& upper(int i) const
   {
      return LPColSetBase<R>::upper(i);
   }

   /// Returns upper bound of column with identifier \p id.
   const R& upper(const SPxColId& id) const
   {
      return LPColSetBase<R>::upper(id);
   }

   /// Returns lower bound vector.
   const VectorBase<R>& lower() const
   {
      return LPColSetBase<R>::lower();
   }

   /// Returns lower bound of column \p i.
   const R& lower(int i) const
   {
      return LPColSetBase<R>::lower(i);
   }

   /// Returns lower bound of column with identifier \p id.
   const R& lower(const SPxColId& id) const
   {
      return LPColSetBase<R>::lower(id);
   }

   /// Returns the optimization sense.
   SPxSense spxSense() const
   {
      return thesense;
   }

   /// Returns the objective function value offset
   const R& objOffset() const
   {
      return offset;
   }

   /// Returns the row number of the row with identifier \p id.
   int number(const SPxRowId& id) const
   {
      return LPRowSetBase<R>::number(id);
   }

   /// Returns the column number of the column with identifier \p id.
   int number(const SPxColId& id) const
   {
      return LPColSetBase<R>::number(id);
   }

   /// Returns the row or column number for identifier \p id.
   int number(const SPxId& id) const
   {
      return (id.type() == SPxId::COL_ID)
         ? LPColSetBase<R>::number(id)
         : LPRowSetBase<R>::number(id);
   }

   /// Returns the row identifier for row \p n.
   SPxRowId rId(int n) const
   {
      return SPxRowId(LPRowSetBase<R>::key(n));
   }

   /// Returns the column identifier for column \p n.
   SPxColId cId(int n) const
   {
      return SPxColId(LPColSetBase<R>::key(n));
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Extension */
   //@{

   ///
   virtual void addRow(const LPRowBase<R>& row)
   {
      doAddRow(row);
   }

   ///
   virtual void addRow(const R& lhsValue, const SVectorBase<R>& rowVec, const R& rhsValue)
   {
      doAddRow(lhsValue, rowVec, rhsValue);
   }

   ///
   template < class S >
   void addRow(const S* lhsValue, const S* rowValues, const int* rowIndices, int rowSize, const S* rhsValue)
   {
      assert(lhsValue != 0);
      assert(rowSize <= 0 || rowValues != 0);
      assert(rowSize <= 0 || rowIndices != 0);
      assert(rhsValue != 0);

      int idx = nRows();
      int oldColNumber = nCols();

      LPRowSetBase<R>::add(lhsValue, rowValues, rowIndices, rowSize, rhsValue);

      // now insert nonzeros to column file also
      for( int j = rowSize - 1; j >= 0; --j )
      {
         const S& val = rowValues[j];
         int i = rowIndices[j];

         // create new columns if required
         if( i >= nCols() )
         {
            LPColBase<R> empty;
            for( int k = nCols(); k <= i; ++k )
               LPColSetBase<R>::add(empty);
         }

         assert(i < nCols());
         LPColSetBase<R>::add2(i, 1, &idx, &val);
      }

      addedRows(1);
      addedCols(nCols() - oldColNumber);
   }

   /// Adds \p row to LPRowSetBase.
   virtual void addRow(SPxRowId& id, const LPRowBase<R>& row)
   {
      addRow(row);
      id = rId(nRows() - 1);
   }

   ///
   virtual void addRows(const LPRowSetBase<R>& pset)
   {
      doAddRows(pset);
   }

   ///
   template < class S >
   void addRows(const S* lhsValues, const S* rowValues, const int* rowIndices, const int* rowStarts, const int* rowLengths, const int numRows, const int numValues, const S* rhsValues)
   {
      assert(lhsValues != 0);
      assert(numValues <= 0 || rowValues != 0);
      assert(numValues <= 0 || rowIndices != 0);
      assert(numValues <= 0 || rowStarts != 0);
      assert(numValues <= 0 || rowLengths != 0);
      assert(rhsValues != 0);

      int i, j, k, idx;
      SVectorBase<R>* col;
      DataArray < int > newCols(nCols());
      int oldRowNumber = nRows();
      int oldColNumber = nCols();

      LPRowSetBase<R>::memRemax(oldRowNumber + numRows);
      for( i = 0; i < numRows; i++ )
      {
         assert(numValues <= 0 || rowStarts[i] + rowLengths[i] <= numValues);
         if( numValues <= 0 )
            LPRowSetBase<R>::add(&(lhsValues[i]), (S*)0, (int*)0, 0, &(rhsValues[i]));
         else
            LPRowSetBase<R>::add(&(lhsValues[i]), &(rowValues[rowStarts[i]]), &(rowIndices[rowStarts[i]]), rowLengths[i], &(rhsValues[i]));
      }

      assert(LPRowSetBase<R>::isConsistent());
      assert(LPColSetBase<R>::isConsistent());

      // count additional nonzeros per column
      for( i = nCols() - 1; i >= 0; --i )
         newCols[i] = 0;
      if( numValues > 0 )
      {
         for( i = 0; i < numRows; i++ )
         {
            for( j = rowStarts[i]; j < rowStarts[i] + rowLengths[i]; j++ )
            {
               ///@todo implement the addition of new columns as in doAddRows()
               assert(rowIndices[j] >= 0);
               assert(rowIndices[j] < oldColNumber);
               newCols[rowIndices[j]]++;
            }
         }
      }

      // extend columns as required (backward because of memory efficiency reasons)
      for( i = nCols() - 1; i >= 0; --i )
      {
         if( newCols[i] > 0 )
         {
            int len = newCols[i] + colVector(i).size();
            LPColSetBase<R>::xtend(i, len);

            /* preset the sizes: beware that this can irritate a consistency check call from xtend(). We need to set the
             * sizes here, because a possible garbage collection called from xtend might destroy the sizes again. */
            colVector_w(i).set_size( len );
         }
      }

      // insert new elements to column file
      for( i = nRows() - 1; i >= oldRowNumber; --i )
      {
         const SVectorBase<R>& vec = rowVector(i);

         for( j = vec.size() - 1; j >= 0; --j )
         {
            k = vec.index(j);
            col = &colVector_w(k);
            idx = col->size() - newCols[k];
            assert(newCols[k] > 0);
            assert(idx >= 0);
            newCols[k]--;
            col->index(idx) = i;
            col->value(idx) = vec.value(j);
         }
      }

#ifndef NDEBUG
      for( i = 0; i < nCols(); ++i )
         assert( newCols[i] == 0 );
#endif

      assert(SPxLPBase<R>::isConsistent());

      assert( numRows == nRows() - oldRowNumber );
      addedRows( nRows() - oldRowNumber );
      addedCols( nCols() - oldColNumber );
   }

   /// adds all LPRowBase%s of \p pset to LPRowSetBase.
   virtual void addRows(SPxRowId id[], const LPRowSetBase<R>& set)
   {

      int i = nRows();
      addRows(set);
      for( int j = 0; i < nRows(); ++i, ++j )
         id[j] = rId(i);
   }

   ///
   virtual void addCol(const LPColBase<R>& col)
   {
      doAddCol(col);
   }

   ///
   virtual void addCol(const R& objValue, const R& lowerValue, const SVectorBase<R>& colVec, const R& upperValue)
   {
      doAddCol(objValue, lowerValue, colVec, upperValue);
   }

   ///
   template < class S >
   void addCol(const S* objValue, const S* lowerValue, const S* colValues, const int* colIndices, int colSize, const S* upperValue)
   {
      int idx = nCols();
      int oldRowNumber = nRows();

      LPColSetBase<R>::add(objValue, lowerValue, colValues, colIndices, colSize, upperValue);
      if( thesense != MAXIMIZE )
         LPColSetBase<R>::maxObj_w(idx) *= -1;

      // now insert nonzeros to column file also
      for( int j = colSize - 1; j >= 0; --j )
      {
         const S& val = colValues[j];
         int i = colIndices[j];

         // create new rows if required
         if( i >= nRows() )
         {
            LPRowBase<R> empty;
            for( int k = nRows(); k <= i; ++k )
               LPRowSetBase<R>::add(empty);
         }

         assert(i < nRows());
         LPRowSetBase<R>::add2(i, 1, &idx, &val);
      }

      addedCols(1);
      addedRows(nRows() - oldRowNumber);
   }

   /// Adds \p col to LPColSetVBase.
   virtual void addCol(SPxColId& id, const LPColBase<R>& col)
   {
      addCol(col);
      id = cId(nCols() - 1);
   }

   ///
   virtual void addCols(const LPColSetBase<R>& pset)
   {
      doAddCols(pset);
   }

   ///
   template < class S >
   void addCols(const S* objValue, const S* lowerValues, const S* colValues, const int* colIndices, const int* colStarts, const int* colLengths, const int numCols, const int numValues, const S* upperValues)
   {
      assert(lowerValues != 0);
      assert(numValues <= 0 || colValues != 0);
      assert(numValues <= 0 || colIndices != 0);
      assert(numValues <= 0 || colStarts != 0);
      assert(numValues <= 0 || colLengths != 0);
      assert(upperValues != 0);

      int i, j, k, idx;
      SVectorBase<R>* row;
      DataArray < int > newRows(nRows());
      int oldColNumber = nCols();
      int oldRowNumber = nRows();
      idx = nCols();

      LPColSetBase<R>::memRemax(oldColNumber + numCols);
      for( i = 0; i < numCols; i++ )
      {
         assert(numValues <= 0 || colStarts[i] + colLengths[i] <= numValues);
         if( numValues <= 0 )
            LPColSetBase<R>::add(&(objValue[i]), &(lowerValues[i]), (S*)0, (int*)0, 0, &(upperValues[i]));
         else
            LPColSetBase<R>::add(&(objValue[i]), &(lowerValues[i]), &(colValues[colStarts[i]]), &(colIndices[colStarts[i]]), colLengths[i], &(upperValues[i]));

         if( thesense != MAXIMIZE )
            LPColSetBase<R>::maxObj_w(idx + i) *= -1;
      }

      assert(LPColSetBase<R>::isConsistent());
      assert(LPRowSetBase<R>::isConsistent());

      // count additional nonzeros per rows
      for( i = nRows() - 1; i >= 0; --i )
         newRows[i] = 0;
      for( i = numValues - 1; i >= 0; --i )
      {
         ///@todo implement the addition of new rows as in doAddCols()
         assert(colIndices[i] >= 0);
         assert(colIndices[i] < oldRowNumber);
         newRows[colIndices[i]]++;
      }

      // extend rows as required (backward because of memory efficiency reasons)
      for( i = nRows() - 1; i >= 0; --i )
      {
         if( newRows[i] > 0 )
         {
            int len = newRows[i] + rowVector(i).size();
            LPRowSetBase<R>::xtend(i, len);

            /* preset the sizes: beware that this can irritate a consistency check call from xtend(). We need to set the
             * sizes here, because a possible garbage collection called from xtend might destroy the sizes again. */
            rowVector_w(i).set_size( len );
         }
      }

      // insert new elements to row file
      for( i = nCols() - 1; i >= oldColNumber; --i )
      {
         const SVectorBase<R>& vec = colVector(i);

         for( j = vec.size() - 1; j >= 0; --j )
         {
            k = vec.index(j);
            row = &rowVector_w(k);
            idx = row->size() - newRows[k];
            assert(newRows[k] > 0);
            assert(idx >= 0);
            newRows[k]--;
            row->index(idx) = i;
            row->value(idx) = vec.value(j);
         }
      }

#ifndef NDEBUG
      for( i = 0; i < nRows(); ++i )
         assert( newRows[i] == 0 );
#endif

      assert(SPxLPBase<R>::isConsistent());

      assert( numCols == nCols() - oldColNumber );
      addedCols( nCols() - oldColNumber );
      addedRows( nRows() - oldRowNumber );
   }

   /// Adds all LPColBase%s of \p set to LPColSetBase.
   virtual void addCols(SPxColId id[], const LPColSetBase<R>& set)
   {

      int i = nCols();
      addCols(set);
      for( int j = 0; i < nCols(); ++i, ++j )
         id[j] = cId(i);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Shrinking */
   //@{

   /// Removes \p i 'th row.
   virtual void removeRow(int i)
   {
      if( i < 0 )
         return;

      doRemoveRow(i);
   }

   /// Removes row with identifier \p id.
   virtual void removeRow(SPxRowId id)
   {
      removeRow(number(id));
   }

   /// Removes multiple rows.
   /** This method removes all LPRowBase%s from the SPxLPBase with an index \p i such that \p perm[i] < 0. Upon
    *  completion, \p perm[i] >= 0 indicates the new index where the \p i'th LPRow has been moved to due to this
    *  removal. Note that \p perm must point to an array of at least #nRows() ints.
    */
   virtual void removeRows(int perm[])
   {
      doRemoveRows(perm);
   }

   ///
   virtual void removeRows(SPxRowId id[], int n, int perm[] = 0)
   {

      if( perm == 0 )
      {
         DataArray < int > p(nRows());
         removeRows(id, n, p.get_ptr());
         return;
      }

      for( int i = nRows() - 1; i >= 0; --i )
         perm[i] = i;

      while( n-- )
         perm[number(id[n])] = -1;

      removeRows(perm);
   }

   /// Removes \p n LPRowBase%s.
   /** Removing multiple rows with one method invocation is available in two flavours. An array \p perm can be passed as
    *  third argument or not. If given, \p perm must be an array at least of size #nRows(). It is used to return the
    *  permutations resulting from this removal: \p perm[i] < 0 indicates, that the element to index \p i has been
    *  removed.  Otherwise, \p perm[i] is the new index of the element with index \p i before the removal.
    */
   virtual void removeRows(int nums[], int n, int perm[] = 0)
   {

      if( perm == 0 )
      {
         DataArray < int > p(nRows());
         removeRows(nums, n, p.get_ptr());
         return;
      }

      for( int i = nRows() - 1; i >= 0; --i )
         perm[i] = i;

      while( n-- )
         perm[nums[n]] = -1;

      removeRows(perm);
   }

   /// Removes rows from \p start to \p end (including both).
   virtual void removeRowRange(int start, int end, int perm[] = 0)
   {

      if( perm == 0 )
      {
         int i = end - start + 1;
         DataArray < int > p(i);

         while( --i >= 0 )
            p[i] = start + i;

         removeRows(p.get_ptr(), end - start + 1);
         return;
      }

      int i;
      for( i = 0; i < start; ++i )
         perm[i] = i;
      for( ; i <= end; ++i )
         perm[i] = -1;
      for( ; i < nRows(); ++i )
         perm[i] = i;

      removeRows(perm);
   }

   /// Removes \p i 'th column.
   virtual void removeCol(int i)
   {
      if( i < 0 )
         return;

      doRemoveCol(i);
   }

   /// Removes column with identifier \p id.
   virtual void removeCol(SPxColId id)
   {
      removeCol(number(id));
   }

   /// Removes multiple columns.
   /** This method removes all LPColBase%s from the SPxLPBase with an index \p i such that \p perm[i] < 0. Upon
    *  completion, \p perm[i] >= 0 indicates the new index where the \p i 'th LPColBase has been moved to due to this
    *  removal. Note, that \p perm must point to an array of at least #nCols() ints.
    */
   virtual void removeCols(int perm[])
   {
      doRemoveCols(perm);
   }

   ///
   virtual void removeCols(SPxColId id[], int n, int perm[] = 0)
   {

      if( perm == 0 )
      {
         DataArray < int > p(nCols());
         removeCols(id, n, p.get_ptr());
         return;
      }

      for( int i = nCols() - 1; i >= 0; --i )
         perm[i] = i;

      while( n-- )
         perm[number(id[n])] = -1;

      removeCols(perm);
   }

   /// Removes \p n LPCols.
   /** Removing multiple columns with one method invocation is available in two flavours. An array \p perm can be passed
    *  as third argument or not. If given, \p perm must be an array at least of size #nCols(). It is used to return the
    *  permutations resulting from this removal: \p perm[i] < 0 indicates, that the element to index \p i has been
    *  removed.  Otherwise, \p perm[i] is the new index of the element with index \p i before the removal.
    */
   virtual void removeCols(int nums[], int n, int perm[] = 0)
   {

      if( perm == 0 )
      {
         DataArray < int > p(nCols());
         removeCols(nums, n, p.get_ptr());
         return;
      }

      for( int i = nCols() - 1; i >= 0; --i )
         perm[i] = i;

      while( n-- )
         perm[nums[n]] = -1;

      removeCols(perm);
   }

   /// Removes columns from \p start to \p end (including both).
   virtual void removeColRange(int start, int end, int perm[] = 0)
   {

      if( perm == 0 )
      {
         int i = end - start + 1;
         DataArray < int > p(i);

         while( --i >= 0 )
            p[i] = start + i;

         removeCols(p.get_ptr(), end - start + 1);
         return;
      }

      int i;
      for( i = 0; i < start; ++i )
         perm[i] = i;
      for( ; i <= end; ++i )
         perm[i] = -1;
      for( ; i < nCols(); ++i )
         perm[i] = i;

      removeCols(perm);
   }

   /// clears the LP.
   virtual void clear()
   {

      LPRowSetBase<R>::clear();
      LPColSetBase<R>::clear();
      thesense = MAXIMIZE;
      offset = 0;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name IO */
   //@{

   /// Reads LP in LP format from input stream \p in.
   virtual bool readLPF(std::istream& in, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// Reads an LP in MPS format from input stream \p in.
   virtual bool readMPS(std::istream& in, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0);

   /// Reads LP in LP or MPS format from input stream \p in.
   /**@param in       input stream.
    * @param rowNames contains after the call the names of the constraints (rows) in the same order as the rows in the
    *                 LP.  Constraints without a name (only possible with LPF files) are automatically assigned a name.
    *                 Maybe 0 if the names are not needed.
    * @param colNames contains after the call the names of the variables (columns) in the same order as the columns in
    *                 the LP.  Maybe 0 if the names are not needed.
    * @param intVars contains after the call the indices of those variables that where marked as beeing integer in the
    *                 file.  Maybe 0 if the information is not needed.
    * @todo Make sure the Id's in the NameSet%s are the same as in the LP.
    */
   virtual bool read(std::istream& in, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars  = 0)
   {
      bool ok;
      char c;

      in.get(c);
      in.putback(c);

      /* MPS starts either with a comment mark '*' or with the keyword 'NAME' at the first column.  LPF starts either
       * with blanks, a comment mark '\' or with the keyword "MAX" or "MIN" in upper or lower case.  There is no
       * possible valid LPF file starting with a '*' or 'N'.
       */
      ok = ((c == '*') || (c == 'N'))
         ? readMPS(in, rowNames, colNames, intVars)
         : readLPF(in, rowNames, colNames, intVars);

      return ok;
   }

   /// Reads LP from a file.
   virtual bool readFile(const char* filename, NameSet* rowNames = 0, NameSet* colNames = 0, DIdxSet* intVars = 0)
   {

      spxifstream file(filename);

      if( !file )
         return false;

      return read(file, rowNames, colNames, intVars);
   }

   /** Writes a file in LP format to \p out. If \p rowNames and \p colNames are \c NULL, default names are used for the
    *  constraints and variables. If \p intVars is not \c NULL, the variables contained in it are marked as integer in
    *  the output.
    */
   virtual void writeLPF(std::ostream&  out, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* p_intvars = 0) const;

   /// Writes a file in MPS format to \p out.
   virtual void writeMPS(std::ostream&  out, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* p_intvars = 0) const;

   /// Write loaded LP to \p filename.
   virtual void writeFile(const char* filename, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* p_intvars = 0) const
   {

      std::ofstream tmp(filename);
      size_t len_f = strlen(filename);

      if( len_f > 4 && filename[len_f-1] == 's' && filename[len_f-2] == 'p' && filename[len_f-3] == 'm' && filename[len_f-4] == '.' )
      {
         writeMPS(tmp, rowNames, colNames, p_intvars);
      }
      else
      {
         writeLPF(tmp, rowNames, colNames, p_intvars);
      }
   }

   /** prints problem statistics */
   void printProblemStatistics(std::ostream& os)
   {
      int countLower = 0;
      int countUpper = 0;
      int countBoxed = 0;
      int countFreeCol = 0;

      int countLhs = 0;
      int countRhs = 0;
      int countRanged = 0;
      int countFreeRow = 0;

      for( int i = 0; i < nCols(); i++ )
      {
         bool hasLower = false;
         bool hasUpper = false;

         if( lower(i) > -infinity )
         {
            countLower++;
            hasLower = true;
         }

         if( upper(i) < infinity )
         {
            countUpper++;
            hasUpper = true;
         }

         if( hasUpper && hasLower )
            countBoxed++;

         if( !hasUpper && !hasLower )
            countFreeCol++;
      }

      for( int i = 0; i < nRows(); i++)
      {
         bool hasRhs = false;
         bool hasLhs = false;

         if( lhs(i) > -infinity )
         {
            countLhs++;
            hasLhs = true;
         }

         if( rhs(i) < infinity )
         {
            countRhs++;
            hasRhs = true;
         }

         if( hasRhs && hasLhs )
            countRanged++;

         if( !hasRhs && !hasLhs )
            countFreeRow++;
      }

      os << "  Columns           : " << nCols() << "\n"
         << "              boxed : " << countBoxed << "\n"
         << "        lower bound : " << countLower << "\n"
         << "        upper bound : " << countUpper << "\n"
         << "               free : " << countFreeCol << "\n"
         << "  Rows              : " << nRows() << "\n"
         << "             ranged : " << countRanged << "\n"
         << "                lhs : " << countLhs << "\n"
         << "                rhs : " << countRhs << "\n"
         << "               free : " << countFreeRow << "\n"
         << "  Nonzeros          : " << nNzos() << "\n"
         << "         per column : " << Real(nNzos()) / Real(nCols()) << "\n"
         << "            per row : " << Real(nNzos()) / Real(nRows()) << "\n"
         << "           sparsity : " << Real(nNzos()) / Real(nCols()) / Real(nRows()) << "\n"
         << "    min. abs. value : " << Real(minAbsNzo()) << "\n"
         << "    max. abs. value : " << Real(maxAbsNzo()) << "\n";
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Manipulation */
   //@{

   /// Changes objective vector to \p newObj.
   virtual void changeObj(const VectorBase<R>& newObj)
   {

      assert(maxObj().dim() == newObj.dim());
      LPColSetBase<R>::maxObj_w() = newObj;
      if( spxSense() == MINIMIZE )
         LPColSetBase<R>::maxObj_w() *= -1;
      assert(isConsistent());
   }

   /// changes \p i 'th objective vector element to \p newVal.
   virtual void changeObj(int i, const R& newVal)
   {
      LPColSetBase<R>::maxObj_w(i) = newVal;
      if( spxSense() == MINIMIZE )
         LPColSetBase<R>::maxObj_w(i) *= -1;
      assert(isConsistent());
   }

   /// changes \p i 'th objective vector element to \p newVal.
   template < class S >
   void changeObj(int i, const S* newVal)
   {
      LPColSetBase<R>::maxObj_w(i) = *newVal;
      if( spxSense() == MINIMIZE )
         LPColSetBase<R>::maxObj_w(i) *= -1;
      assert(isConsistent());
   }

   /// Changes objective value of column with identifier \p id to \p newVal.
   virtual void changeObj(SPxColId id, const R& newVal)
   {
      changeObj(number(id), newVal);
   }

   /// Changes objective vector to \p newObj.
   virtual void changeMaxObj(const VectorBase<R>& newObj)
   {
      assert(maxObj().dim() == newObj.dim());
      LPColSetBase<R>::maxObj_w() = newObj;
      assert(isConsistent());
   }

   /// changes \p i 'th objective vector element to \p newVal.
   virtual void changeMaxObj(int i, const R& newVal)
   {
      LPColSetBase<R>::maxObj_w(i) = newVal;
      assert(isConsistent());
   }

   /// changes \p i 'th objective vector element to \p newVal.
   template < class S >
   void changeMaxObj(int i, const S* newVal)
   {
      LPColSetBase<R>::maxObj_w(i) = *newVal;
      assert(isConsistent());
   }

   /// Changes objective value of column with identifier \p id to \p newVal.
   virtual void changeMaxObj(SPxColId id, const R& newVal)
   {
      changeMaxObj(number(id), newVal);
   }

   /// Changes vector of lower bounds to \p newLower.
   virtual void changeLower(const VectorBase<R>& newLower)
   {

      assert(lower().dim() == newLower.dim());
      LPColSetBase<R>::lower_w() = newLower;
      assert(isConsistent());
   }

   /// changes \p i 'th lower bound to \p newLower.
   virtual void changeLower(int i, const R& newLower)
   {

      LPColSetBase<R>::lower_w(i) = newLower;
      assert(isConsistent());
   }

   /// changes \p i 'th lower bound to \p newLower.
   template < class S >
   void changeLower(int i, const S* newLower)
   {
      LPColSetBase<R>::lower_w(i) = *newLower;
      assert(isConsistent());
   }

   /// changes lower bound of column with identifier \p id to \p newLower.
   virtual void changeLower(SPxColId id, const R& newLower)
   {
      changeLower(number(id), newLower);
   }

   /// Changes vector of upper bounds to \p newUpper.
   virtual void changeUpper(const VectorBase<R>& newUpper)
   {

      assert(upper().dim() == newUpper.dim());
      LPColSetBase<R>::upper_w() = newUpper;
      assert(isConsistent());
   }

   /// Changes \p i 'th upper bound to \p newUpper.
   virtual void changeUpper(int i, const R& newUpper)
   {

      LPColSetBase<R>::upper_w(i) = newUpper;
      assert(isConsistent());
   }

   /// Changes \p i 'th upper bound to \p newUpper.
   template < class S >
   void changeUpper(int i, const S* newUpper)
   {
      LPColSetBase<R>::upper_w(i) = *newUpper;
      assert(isConsistent());
   }

   /// Changes upper bound of column with identifier \p id to \p newLower.
   virtual void changeUpper(SPxColId id, const R& newUpper)
   {
      changeUpper(number(id), newUpper);
   }

   /// Changes variable bounds to \p newLower and \p newUpper.
   virtual void changeBounds(const VectorBase<R>& newLower, const VectorBase<R>& newUpper)
   {

      changeLower(newLower);
      changeUpper(newUpper);
      assert(isConsistent());
   }

   /// Changes bounds of column \p i to \p newLower and \p newUpper.
   virtual void changeBounds(int i, const R& newLower, const R& newUpper)
   {

      changeLower(i, newLower);
      changeUpper(i, newUpper);
      assert(isConsistent());
   }

   /// Changes bounds of column \p i to \p newLower and \p newUpper.
   template < class S >
   void changeBounds(int i, const S* newLower, const S* newUpper)
   {
      LPColSetBase<R>::lower_w(i) = *newLower;
      LPColSetBase<R>::upper_w(i) = *newUpper;
      assert(isConsistent());
   }

   /// Changes bounds of column with identifier \p id.
   virtual void changeBounds(SPxColId id, const R& newLower, const R& newUpper)
   {
      changeBounds(number(id), newLower, newUpper);
   }

   /// Changes left hand side vector for constraints to \p newLhs.
   virtual void changeLhs(const VectorBase<R>& newLhs)
   {
      assert(lhs().dim() == newLhs.dim());
      LPRowSetBase<R>::lhs_w() = newLhs;
      assert(isConsistent());
   }

   /// Changes \p i 'th left hand side value to \p newLhs.
   virtual void changeLhs(int i, const R& newLhs)
   {
      LPRowSetBase<R>::lhs_w(i) = newLhs;
      assert(isConsistent());
   }

   /// Changes \p i 'th left hand side value to \p newLhs.
   template < class S >
   void changeLhs(int i, const S* newLhs)
   {
      LPRowSetBase<R>::lhs_w(i) = *newLhs;
      assert(isConsistent());
   }

   /// Changes left hand side value for row with identifier \p id.
   virtual void changeLhs(SPxRowId id, const R& newLhs)
   {
      changeLhs(number(id), newLhs);
   }

   /// Changes right hand side vector for constraints to \p newRhs.
   virtual void changeRhs(const VectorBase<R>& newRhs)
   {
      assert(rhs().dim() == newRhs.dim());
      LPRowSetBase<R>::rhs_w() = newRhs;
      assert(isConsistent());
   }

   /// Changes \p i 'th right hand side value to \p newRhs.
   virtual void changeRhs(int i, const R& newRhs)
   {
      LPRowSetBase<R>::rhs_w(i) = newRhs;
      assert(isConsistent());
   }

   /// Changes right hand side value for row with identifier \p id.
   virtual void changeRhs(SPxRowId id, const R& newRhs)
   {
      changeRhs(number(id), newRhs);
   }

   /// Changes left and right hand side vectors.
   virtual void changeRange(const VectorBase<R>& newLhs, const VectorBase<R>& newRhs)
   {

      changeLhs(newLhs);
      changeRhs(newRhs);
      assert(isConsistent());
   }

   /// Changes left and right hand side of row \p i.
   virtual void changeRange(int i, const R& newLhs, const R& newRhs)
   {

      changeLhs(i, newLhs);
      changeRhs(i, newRhs);
      assert(isConsistent());
   }

   /// Changes left and right hand side of row \p i.
   template < class S >
   void changeRange(int i, const S* newLhs, const S* newRhs)
   {
      LPRowSetBase<R>::lhs_w(i) = *newLhs;
      LPRowSetBase<R>::rhs_w(i) = *newRhs;
      assert(isConsistent());
   }

   /// Changes left and right hand side of row with identifier \p id.
   virtual void changeRange(SPxRowId id, const R& newLhs, const R& newRhs)
   {
      changeRange(number(id), newLhs, newRhs);
   }

   /// Changes row objective function vector to \p newRowObj.
   virtual void changeRowObj(const VectorBase<R>& newRowObj)
   {
      assert(maxRowObj().dim() == newRowObj.dim());
      LPRowSetBase<R>::obj_w() = newRowObj;
      if( spxSense() == MINIMIZE )
         LPRowSetBase<R>::obj_w() *= -1;
      assert(isConsistent());
   }

   /// Changes \p i 'th row objective function value to \p newRowObj.
   virtual void changeRowObj(int i, const R& newRowObj)
   {
      LPRowSetBase<R>::obj_w(i) = newRowObj;
      if( spxSense() == MINIMIZE )
         LPRowSetBase<R>::obj_w(i) *= -1;
      assert(isConsistent());
   }

   /// Changes row objective function value for row with identifier \p id.
   virtual void changeRowObj(SPxRowId id, const R& newRowObj)
   {
      changeRowObj(number(id), newRowObj);
   }

   /// Clears row objective function values for all rows
   virtual void clearRowObjs()
   {
      LPRowSetBase<R>::obj_w().clear();
   }

   /// Replaces \p i 'th row of LP with \p newRow.
   virtual void changeRow(int n, const LPRowBase<R>& newRow)
   {
      if( n < 0 )
         return;

      int j;
      SVectorBase<R>& row = rowVector_w(n);
      for( j = row.size() - 1; j >= 0; --j )
      {
         SVectorBase<R>& col = colVector_w(row.index(j));
         col.remove(col.number(n));
      }

      row.clear();

      changeLhs(n, newRow.lhs());
      changeRhs(n, newRow.rhs());
      changeRowObj(n, newRow.obj());

      const SVectorBase<R>& newrow = newRow.rowVector();
      for( j = newrow.size() - 1; j >= 0; --j )
      {
         int idx = newrow.index(j);
         R val = newrow.value(j);
         LPRowSetBase<R>::add2(n, 1, &idx, &val);
         LPColSetBase<R>::add2(idx, 1, &n, &val);
      }

      assert(isConsistent());
   }

   /// Replaces row with identifier \p id with \p newRow.
   virtual void changeRow(SPxRowId id, const LPRowBase<R>& newRow)
   {
      changeRow(number(id), newRow);
   }

   /// Replaces \p i 'th column of LP with \p newCol.
   virtual void changeCol(int n, const LPColBase<R>& newCol)
   {
      if( n < 0 )
         return;

      int j;
      SVectorBase<R>& col = colVector_w(n);
      for( j = col.size() - 1; j >= 0; --j )
      {
         SVectorBase<R>& row = rowVector_w(col.index(j));
         row.remove(row.number(n));
      }

      col.clear();

      changeUpper(n, newCol.upper());
      changeLower(n, newCol.lower());
      changeObj(n, newCol.obj());

      const SVectorBase<R>& newcol = newCol.colVector();
      for( j = newcol.size() - 1; j >= 0; --j )
      {
         int idx = newcol.index(j);
         R val = newcol.value(j);
         LPColSetBase<R>::add2(n, 1, &idx, &val);
         LPRowSetBase<R>::add2(idx, 1, &n, &val);
      }

      assert(isConsistent());
   }

   /// Replaces column with identifier \p id with \p newCol.
   virtual void changeCol(SPxColId id, const LPColBase<R>& newCol)
   {
      changeCol(number(id), newCol);
   }

   /// Changes LP element (\p i, \p j) to \p val.
   virtual void changeElement(int i, int j, const R& val)
   {
      if( i < 0 || j < 0 )
         return;

      SVectorBase<R>& row = rowVector_w(i);
      SVectorBase<R>& col = colVector_w(j);

      if( val != R(0) )
      {
         if( row.number(j) >= 0 )
         {
            row.value(row.number(j)) = val;
            col.value(col.number(i)) = val;
         }
         else
         {
            LPRowSetBase<R>::add2(i, 1, &j, &val);
            LPColSetBase<R>::add2(j, 1, &i, &val);
         }
      }
      else if( row.number(j) >= 0 )
      {
         row.remove(row.number(j));
         col.remove(col.number(i));
      }

      assert(isConsistent());
   }

   /// Changes LP element (\p i, \p j) to \p val.
   template < class S >
   void changeElement(int i, int j, const S* val)
   {
      if( i < 0 || j< 0 )
         return;

      SVectorBase<R>& row = rowVector_w(i);
      SVectorBase<R>& col = colVector_w(j);

      if( mpq_get_d(*val) != R(0) )
      {
         if( row.number(j) >= 0 )
         {
            row.value(row.number(j)) = *val;
            col.value(col.number(i)) = *val;
         }
         else
         {
            LPRowSetBase<R>::add2(i, 1, &j, val);
            LPColSetBase<R>::add2(j, 1, &i, val);
         }
      }
      else if( row.number(j) >= 0 )
      {
         row.remove(row.number(j));
         col.remove(col.number(i));
      }

      assert(isConsistent());
   }

   /// Changes LP element identified by (\p rid, \p cid) to \p val.
   virtual void changeElement(SPxRowId rid, SPxColId cid, const R& val)
   {
      changeElement(number(rid), number(cid), val);
   }

   /// Changes optimization sense to \p sns.
   virtual void changeSense(SPxSense sns)
   {
      if( sns != thesense )
      {
         LPColSetBase<R>::maxObj_w() *= -1;
         LPRowSetBase<R>::obj_w() *= -1;
      }
      thesense = sns;
   }

   virtual void changeObjOffset(const R& o)
   {
      offset = o;
   }

   /// Computes activity of the rows for a given primal vector; activity does not need to be zero
   /// @throw SPxInternalCodeException if the dimension of primal vector does not match number of columns or if the
   ///        dimension of the activity vector does not match the number of rows
   virtual void computePrimalActivity(const VectorBase<R>& primal, VectorBase<R>& activity) const
   {
      if( primal.dim() != nCols() )
      {
         throw SPxInternalCodeException("XSPXLP01 Primal vector for computing row activity has wrong dimension");
      }

      if( activity.dim() != nRows() )
      {
         throw SPxInternalCodeException("XSPXLP03 Activity vector computing row activity has wrong dimension");
      }

      int c;
      for( c = 0; c < nCols() && primal[c] == 0; c++ )
         ;

      if( c >= nCols() )
      {
         activity.clear();
         return;
      }

      activity = colVector(c);
      activity *= primal[c];
      c++;

      for( ; c < nCols(); c++ )
      {
         if( primal[c] != 0 )
            activity.multAdd(primal[c], colVector(c));
      }
   }

   /// Updates activity of the rows for a given primal vector; activity does not need to be zero
   /// @throw SPxInternalCodeException if the dimension of primal vector does not match number of columns or if the
   ///        dimension of the activity vector does not match the number of rows
   virtual void addPrimalActivity(const SVectorBase<R>& primal, VectorBase<R>& activity) const
   {
      if( activity.dim() != nRows() )
      {
         throw SPxInternalCodeException("XSPXLP03 Activity vector computing row activity has wrong dimension");
      }

      for( int i = primal.size() - 1; i >= 0; i-- )
      {
         assert(primal.index(i) >= 0);
         assert(primal.index(i) < nCols());
         activity.multAdd(primal.value(i), colVector(primal.index(i)));
      }
   }

   /// Computes "dual" activity of the columns for a given dual vector, i.e., y^T A; activity does not need to be zero
   /// @throw SPxInternalCodeException if dimension of dual vector does not match number of rows or if the dimension of
   ///        the activity vector does not match the number of columns
   virtual void computeDualActivity(const VectorBase<R>& dual, VectorBase<R>& activity) const
   {
      if( dual.dim() != nRows() )
      {
         throw SPxInternalCodeException("XSPXLP02 Dual vector for computing dual activity has wrong dimension");
      }

      if( activity.dim() != nCols() )
      {
         throw SPxInternalCodeException("XSPXLP04 Activity vector computing dual activity has wrong dimension");
      }

      int r;
      for( r = 0; r < nRows() && dual[r] == 0; r++ )
         ;

      if( r >= nRows() )
      {
         activity.clear();
         return;
      }

      activity = rowVector(r);
      activity *= dual[r];
      r++;

      for( ; r < nRows(); r++ )
      {
         if( dual[r] != 0 )
            activity.multAdd(dual[r], rowVector(r));
      }
   }

   /// Updates "dual" activity of the columns for a given dual vector, i.e., y^T A; activity does not need to be zero
   /// @throw SPxInternalCodeException if dimension of dual vector does not match number of rows or if the dimension of
   ///        the activity vector does not match the number of columns
   virtual void addDualActivity(const SVectorBase<R>& dual, VectorBase<R>& activity) const
   {
      if( activity.dim() != nCols() )
      {
         throw SPxInternalCodeException("XSPXLP04 Activity vector computing dual activity has wrong dimension");
      }

      for( int i = dual.size() - 1; i >= 0; i-- )
      {
         assert(dual.index(i) >= 0);
         assert(dual.index(i) < nRows());
         activity.multAdd(dual.value(i), rowVector(dual.index(i)));
      }
   }

   /// Updates "dual" activity of the columns for a given dual vector, i.e., y^T A; activity does not need to be zero
   /// @throw SPxInternalCodeException if dimension of dual vector does not match number of rows or if the dimension of
   ///        the activity vector does not match the number of columns
   virtual void subDualActivity(const VectorBase<R>& dual, VectorBase<R>& activity) const
   {
      if( dual.dim() != nRows() )
      {
         throw SPxInternalCodeException("XSPXLP02 Dual vector for computing dual activity has wrong dimension");
      }

      if( activity.dim() != nCols() )
      {
         throw SPxInternalCodeException("XSPXLP04 Activity vector computing dual activity has wrong dimension");
      }

      for( int r = 0; r < nRows(); r++ )
      {
         if( dual[r] != 0 )
            activity.multSub(dual[r], rowVector(r));
      }
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Miscellaneous */
   //@{

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS

      for( int i = nCols() - 1; i >= 0; --i )
      {
         const SVectorBase<R>& v = colVector(i);

         for( int j = v.size() - 1; j >= 0; --j )
         {
            const SVectorBase<R>& w = rowVector(v.index(j));
            int n = w.number(i);

            if( n < 0 )
               return MSGinconsistent("SPxLPBase");

            if( v.value(j) != w.value(n) )
               return MSGinconsistent("SPxLPBase");
         }
      }

      for( int i = nRows() - 1; i >= 0; --i )
      {
         const SVectorBase<R>& v = rowVector(i);

         for( int j = v.size() - 1; j >= 0; --j )
         {
            const SVectorBase<R>& w = colVector(v.index(j));
            int n = w.number(i);

            if( n < 0 )
               return MSGinconsistent("SPxLPBase");

            if( v.value(j) != w.value(n) )
               return MSGinconsistent("SPxLPBase");
         }
      }

      return LPRowSetBase<R>::isConsistent() && LPColSetBase<R>::isConsistent();
#else
      return true;
#endif
   }

   //@}

protected:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Protected write access */
   //@{

   /// Returns right hand side of row \p i.
   R& rhs_w(int i)
   {
      return LPRowSetBase<R>::rhs_w(i);
   }

   /// Returns left hand side of row \p i.
   R& lhs_w(int i)
   {
      return LPRowSetBase<R>::lhs_w(i);
   }

   /// Returns objective function value of row \p i.
   R& maxRowObj_w(int i)
   {
      return LPRowSetBase<R>::obj_w(i);
   }

   /// Returns objective value of column \p i for maximization problem.
   R& maxObj_w(int i)
   {
      return LPColSetBase<R>::maxObj_w(i);
   }

   /// Returns upper bound of column \p i.
   R& upper_w(int i)
   {
      return LPColSetBase<R>::upper_w(i);
   }

   /// Returns lower bound of column \p i.
   R& lower_w(int i)
   {
      return LPColSetBase<R>::lower_w(i);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Protected helpers */
   //@{

   /// Returns the LP as an LPRowSetBase.
   const LPRowSetBase<R>* lprowset() const
   {
      return static_cast<const LPRowSetBase<R>*>(this);
   }

   /// Returns the LP as an LPColSetBase.
   const LPColSetBase<R>* lpcolset() const
   {
      return static_cast<const LPColSetBase<R>*>(this);
   }

   /// Internal helper method.
   virtual void doRemoveRow(int j)
   {

      const SVectorBase<R>& vec = rowVector(j);

      // remove row vector from column file
      for( int i = vec.size() - 1; i >= 0; --i )
      {
         SVectorBase<R>& remvec = colVector_w(vec.index(i));
         remvec.remove(remvec.number(j));
      }

      // move last row to removed position
      int idx = nRows() - 1;
      if( j != idx )
      {
         const SVectorBase<R>& l_vec = rowVector(idx);
         for( int i = l_vec.size() - 1; i >= 0; --i )
         {
            SVectorBase<R>& movevec = colVector_w(l_vec.index(i));
            movevec.index(movevec.number(idx)) = j;
         }
      }

      LPRowSetBase<R>::remove(j);
   }

   /// Internal helper method.
   virtual void doRemoveRows(int perm[])
   {

      int j = nCols();

      LPRowSetBase<R>::remove(perm);

      for( int i = 0; i < j; ++i )
      {
         SVectorBase<R>& vec = colVector_w(i);
         for( int k = vec.size() - 1; k >= 0; --k )
         {
            int idx = vec.index(k);
            if( perm[idx] < 0 )
               vec.remove(k);
            else
               vec.index(k) = perm[idx];
         }
      }
   }

   /// Internal helper method.
   virtual void doRemoveCol(int j)
   {

      const SVectorBase<R>& vec = colVector(j);
      int i;

      // remove column vector from row file
      for( i = vec.size() - 1; i >= 0; --i )
      {
         SVectorBase<R>& remvec = rowVector_w(vec.index(i));
         remvec.remove(remvec.number(j));
      }

      // move last column to removed position
      int idx = nCols() - 1;
      if( j != idx )
      {
         const SVectorBase<R>& l_vec = colVector(idx);
         for( i = l_vec.size() - 1; i >= 0; --i )
         {
            SVectorBase<R>& movevec = rowVector_w(l_vec.index(i));
            movevec.index(movevec.number(idx)) = j;
         }
      }

      LPColSetBase<R>::remove(j);
   }

   /// Internal helper method.
   virtual void doRemoveCols(int perm[])
   {

      int j = nRows();

      LPColSetBase<R>::remove(perm);

      for( int i = 0; i < j; ++i )
      {
         SVectorBase<R>& vec = rowVector_w(i);

         for( int k = vec.size() - 1; k >= 0; --k )
         {
            int idx = vec.index(k);
            if( perm[idx] < 0 )
               vec.remove(k);
            else
               vec.index(k) = perm[idx];
         }
      }
   }

   /// Called after the last \p n rows have just been added.
   virtual void addedRows(int)
   {}

   /// Called after the last \p n columns have just been added.
   virtual void addedCols(int)
   {}

   ///
   void added2Set(SVSetBase<R>& set, const SVSetBase<R>& addset, int n)
   {

      if( n == 0 )
         return;

      DataArray<int> moreArray(set.num());
      int* more = moreArray.get_ptr();

      for( int i = set.num() - 1; i >= 0; --i )
         more[i] = 0;

      int tot = 0;
      int end = addset.num();

      for( int i = addset.num() - n; i < end; ++i )
      {
         const SVectorBase<R>& vec = addset[i];

         tot += vec.size();
         for( int j = vec.size() - 1; j >= 0; --j )
            more[vec.index(j)]++;
      }

      if( set.memMax() < tot )
         set.memRemax(tot);

      for( int i = set.num() - 1; i >= 0; --i )
      {
         int j = set[i].size();
         set.xtend(set[i], j + more[i]);
         set[i].set_size( j + more[i] );
         more[i] = j;
      }

      for( int i = addset.num() - n; i < addset.num(); ++i)
      {
         const SVectorBase<R>& vec = addset[i];

         for( int j = vec.size() - 1; j >= 0; --j )
         {
            int k = vec.index(j);
            int m = more[k]++;
            SVectorBase<R>& l_xtend = set[k];
            l_xtend.index(m) = i;
            l_xtend.value(m) = vec.value(j);
         }
      }
   }

   //@}


private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Private helpers */
   //@{

   /// Returns the LP as an LPRowSet.
   SVectorBase<R>& colVector_w(int i)
   {
      return LPColSetBase<R>::colVector_w(i);
   }

   ///
   SVectorBase<R>& rowVector_w(int i)
   {
      return LPRowSetBase<R>::rowVector_w(i);
   }

   ///
   void doAddRow (const LPRowBase<R>& row)
   {

      int idx = nRows();
      int oldColNumber = nCols();
      const SVectorBase<R>& vec = row.rowVector();

      LPRowSetBase<R>::add(row);

      // now insert nonzeros to column file also
      for( int j = vec.size() - 1; j >= 0; --j )
      {
         R val = vec.value(j);
         int i = vec.index(j);

         // create new columns if required
         if( i >= nCols() )
         {
            LPColBase<R> empty;
            for( int k = nCols(); k <= i; ++k )
               LPColSetBase<R>::add(empty);
         }

         assert(i < nCols());
         LPColSetBase<R>::add2(i, 1, &idx, &val);
      }

      addedRows(1);
      addedCols(nCols() - oldColNumber);
   }

   ///
   void doAddRow (const R& lhsValue, const SVectorBase<R>& rowVec, const R& rhsValue)
   {
      int idx = nRows();
      int oldColNumber = nCols();

      LPRowSetBase<R>::add(lhsValue, rowVec, rhsValue);

      // now insert nonzeros to column file also
      for( int j = rowVec.size() - 1; j >= 0; --j )
      {
         R val = rowVec.value(j);
         int i = rowVec.index(j);

         // create new columns if required
         if( i >= nCols() )
         {
            LPColBase<R> empty;
            for( int k = nCols(); k <= i; ++k )
               LPColSetBase<R>::add(empty);
         }

         assert(i < nCols());
         LPColSetBase<R>::add2(i, 1, &idx, &val);
      }

      addedRows(1);
      addedCols(nCols() - oldColNumber);
   }

   ///
   void doAddRows(const LPRowSetBase<R>& set)
   {

      int i, j, k, ii, idx;
      SVectorBase<R>* col;
      DataArray < int > newCols(nCols());
      int oldRowNumber = nRows();
      int oldColNumber = nCols();

      if( &set != this )
         LPRowSetBase<R>::add(set);

      assert(LPRowSetBase<R>::isConsistent());
      assert(LPColSetBase<R>::isConsistent());

      // count additional nonzeros per column
      for( i = nCols() - 1; i >= 0; --i )
         newCols[i] = 0;
      for( i = set.num() - 1; i >= 0; --i )
      {
         const SVectorBase<R>& vec = set.rowVector(i);

         for( j = vec.size() - 1; j >= 0; --j )
         {
            // create new columns if required
            ii = vec.index(j);
            if( ii >= nCols() )
            {
               LPColBase<R> empty;
               newCols.reSize(ii + 1);
               for( k = nCols(); k <= ii; ++k )
               {
                  newCols[k] = 0;
                  LPColSetBase<R>::add(empty);
               }
            }

            assert(ii < nCols());
            newCols[ii]++;
         }
      }

      // extend columns as required (backward because of memory efficiency reasons)
      for( i = nCols() - 1; i >= 0; --i )
      {
         if( newCols[i] > 0 )
         {
            int len = newCols[i] + colVector(i).size();
            LPColSetBase<R>::xtend(i, len);

            /* preset the sizes: beware that this can irritate a consistency check call from xtend(). We need to set the
             * sizes here, because a possible garbage collection called from xtend might destroy the sizes again. */
            colVector_w(i).set_size( len );
         }
      }

      // insert new elements to column file
      for( i = nRows() - 1; i >= oldRowNumber; --i )
      {
         const SVectorBase<R>& vec = rowVector(i);

         for( j = vec.size() - 1; j >= 0; --j )
         {
            k = vec.index(j);
            col = &colVector_w(k);
            idx = col->size() - newCols[k];
            assert(newCols[k] > 0);
            assert(idx >= 0);
            newCols[k]--;
            col->index(idx) = i;
            col->value(idx) = vec.value(j);
         }
      }

#ifndef NDEBUG
      for( i = 0; i < nCols(); ++i )
         assert( newCols[i] == 0 );
#endif

      assert(SPxLPBase<R>::isConsistent());

      assert( set.num() == nRows() - oldRowNumber );
      addedRows( nRows() - oldRowNumber );
      addedCols( nCols() - oldColNumber );
   }

   ///
   void doAddCol (const LPColBase<R>& col)
   {

      int idx = nCols();
      int oldRowNumber = nRows();
      const SVectorBase<R>& vec = col.colVector();

      LPColSetBase<R>::add(col);
      if( thesense != MAXIMIZE )
         LPColSetBase<R>::maxObj_w(idx) *= -1;

      // now insert nonzeros to row file also
      for( int j = vec.size() - 1; j >= 0; --j )
      {
         R val = vec.value(j);
         int i = vec.index(j);

         // create new rows if required
         if( i >= nRows() )
         {
            LPRowBase<R> empty;
            for( int k = nRows(); k <= i; ++k )
               LPRowSetBase<R>::add(empty);
         }

         assert(i < nRows());
         LPRowSetBase<R>::add2(i, 1, &idx, &val);
      }

      addedCols(1);
      addedRows(nRows() - oldRowNumber);
   }

   ///
   void doAddCol (const R& objValue, const R& lowerValue, const SVectorBase<R>& colVec, const R& upperValue)
   {
      int idx = nCols();
      int oldRowNumber = nRows();

      LPColSetBase<R>::add(objValue, lowerValue, colVec, upperValue);
      if( thesense != MAXIMIZE )
         LPColSetBase<R>::maxObj_w(idx) *= -1;

      // now insert nonzeros to row file also
      for( int j = colVec.size() - 1; j >= 0; --j )
      {
         R val = colVec.value(j);
         int i = colVec.index(j);

         // create new rows if required
         if( i >= nRows() )
         {
            LPRowBase<R> empty;
            for( int k = nRows(); k <= i; ++k )
               LPRowSetBase<R>::add(empty);
         }

         assert(i < nRows());
         LPRowSetBase<R>::add2(i, 1, &idx, &val);
      }

      addedCols(1);
      addedRows(nRows() - oldRowNumber);
   }

   ///
   void doAddCols(const LPColSetBase<R>& set)
   {

      int i, j;
      int oldColNumber = nCols();
      int oldRowNumber = nRows();
      DataArray < int > newRows(nRows());

      if( &set != this )
         LPColSetBase<R>::add(set);

      assert(LPColSetBase<R>::isConsistent());
      assert(LPRowSetBase<R>::isConsistent());

      // count additional nonzeros per row
      for( i = nRows() - 1; i >= 0; --i )
         newRows[i] = 0;

      for( i = set.num() - 1; i >= 0; --i )
      {
         const SVectorBase<R>& vec = set.colVector(i);

         for( j = vec.size() - 1; j >= 0; --j )
         {
            // create new rows if required
            int l = vec.index(j);
            if( l >= nRows() )
            {
               LPRowBase<R> empty;
               newRows.reSize(l + 1);
               for( int k = nRows(); k <= l; ++k )
               {
                  newRows[k] = 0;
                  LPRowSetBase<R>::add(empty);
               }

            }

            assert(l < nRows());
            newRows[l]++;
         }
      }

      // extend rows as required
      for( i = 0; i < nRows(); ++i )
      {
         if( newRows[i] > 0 )
         {
            int len = newRows[i] + rowVector(i).size();
            LPRowSetBase<R>::xtend(i, len);
            rowVector_w(i).set_size( len );
         }
      }

      // insert new elements to row file
      for( i = oldColNumber; i < nCols(); ++i )
      {
         LPColSetBase<R>::maxObj_w(i) *= thesense;
         const SVectorBase<R>& vec = colVector(i);

         for( j = vec.size() - 1; j >= 0; --j )
         {
            int k = vec.index(j);
            SVectorBase<R>& row = rowVector_w(k);
            int idx = row.size() - newRows[k];
            assert(newRows[k] > 0);
            newRows[k]--;
            row.index(idx) = i;
            row.value(idx) = vec.value(j);
         }
      }

#ifndef NDEBUG
      for( i = 0; i < nRows(); ++i )
         assert( newRows[i] == 0 );
#endif

      assert(SPxLPBase<R>::isConsistent());

      assert(set.num() == nCols() - oldColNumber);
      addedCols(nCols() - oldColNumber);
      addedRows(nRows() - oldRowNumber);
   }

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Constructors / Destructors */
   //@{

   /// Default constructor.
   SPxLPBase<R>()
   {
      SPxLPBase<R>::clear(); // clear is virtual.

      assert(isConsistent());
   }

   /// Destructor.
   virtual ~SPxLPBase<R>()
   {}

   /// Copy constructor.
   SPxLPBase<R>(const SPxLPBase<R>& old)
      : LPRowSetBase<R>(old)
      , LPColSetBase<R>(old)
      , thesense(old.thesense)
      , offset(old.offset)
      , spxout(old.spxout)
   {
      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   SPxLPBase<R>(const SPxLPBase<S>& old)
      : LPRowSetBase<R>(old)
      , LPColSetBase<R>(old)
      , thesense(old.thesense == SPxLPBase<S>::MINIMIZE ? SPxLPBase<R>::MINIMIZE : SPxLPBase<R>::MAXIMIZE)
      , offset(old.offset)
      , spxout(old.spxout)
   {
      assert(isConsistent());
   }

   /// Assignment operator.
   SPxLPBase<R>& operator=(const SPxLPBase<R>& old)
   {
      if( this != &old )
      {
         LPRowSetBase<R>::operator=(old);
         LPColSetBase<R>::operator=(old);
         thesense = old.thesense;
         offset = old.offset;

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   SPxLPBase<R>& operator=(const SPxLPBase<S>& old)
   {
      if( this != (const SPxLPBase<R>*)(&old) )
      {
         LPRowSetBase<R>::operator=(old);
         LPColSetBase<R>::operator=(old);
         thesense = (old.thesense) == SPxLPBase<S>::MINIMIZE ? SPxLPBase<R>::MINIMIZE : SPxLPBase<R>::MAXIMIZE;
         offset = R(old.offset);

         assert(isConsistent());
      }

      return *this;
   }

   //@}
};

} // namespace soplex

/* reset the SOPLEX_DEBUG flag to its original value */
#undef SOPLEX_DEBUG
#ifdef SOPLEX_DEBUG_SPXLPBASE
#define SOPLEX_DEBUG
#undef SOPLEX_DEBUG_SPXLPBASE
#endif

#endif // _SPXLPBASE_H_
