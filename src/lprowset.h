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
#pragma ident "@(#) $Id: lprowset.h,v 1.12 2002/01/19 16:05:25 bzfkocht Exp $"

/**@file  lprowset.h
 * @brief Set of LP columns.
 */
#ifndef _LPROWSET_H_
#define _LPROWSET_H_


#include <assert.h>

#include "lprow.h"
#include "dvector.h"
#include "svset.h"
#include "datakey.h"

namespace soplex
{
/**@brief   Set of LP rows.
   @ingroup Algebra
   
   Class #LPRowSet implements a set of #LPRow%s. Unless for memory
   limitations, any number of #LPRow%s may be #add%ed to an #LPRowSet. Single
   or multiple #LPRow%s may be #add%ed to an #LPRowSet, where each method
   #add() comes with two different signatures. One with and one without a
   parameter, used for returning the #Key%s assigned to the new #LPRow%s
   by the set. See #DataSet::Key for a more detailed description of the
   concept of #Key%s. For the concept of renumbering #LPRow%s within an
   #LPRowSet after removal of some #LPRow%s see #DataSet.
   
   @see        DataSet, DataSet::Key
*/
class LPRowSet : protected SVSet
{
private:
   DVector left;  ///< vector of left hand sides (lower bounds) of #LPRow%s.
   DVector right; ///< vector of right hand sides (upper bounds) of #LPRow%s.

public:

   /**@name Inquiry */
   //@{
   /// returns the number of #LPRow%s in #LPRowSet.
   int num() const
   {
      return SVSet::num();
   }

   /// returns the maximum number of #LPRow%s that fit.
   int max() const
   {
      return SVSet::max();
   }

   ///
   const Vector& lhs() const
   {
      return left;
   }
   /// returns the vector of #lhs values.
   Vector& lhs()
   {
      return left;
   }

   ///
   double lhs(int i) const
   {
      return left[i];
   }
   /// returns the #lhs of the \p i 'th #LPRow.
   double& lhs(int i)
   {
      return left[i];
   }

   ///
   double lhs(const DataKey& k) const
   {
      return left[number(k)];
   }
   /// returns the #lhs of the #LPRow with #DataKey \k in #LPRowSet.
   double& lhs(const DataKey& k)
   {
      return left[number(k)];
   }

   ///
   const Vector& rhs() const
   {
      return right;
   }
   /// returns the vector of #rhs values.
   Vector& rhs()
   {
      return right;
   }

   ///
   double rhs(int i) const
   {
      return right[i];
   }
   /// returns the #rhs of the \p i 'th #LPRow.
   double& rhs(int i)
   {
      return right[i];
   }

   ///
   double rhs(const DataKey& k) const
   {
      return right[number(k)];
   }
   /// returns the #rhs of the #LPRow with #DataKey \k in #LPRowSet.
   double& rhs(const DataKey& k)
   {
      return right[number(k)];
   }

   /// returns a writable #rowVector of the \p i 'th #LPRow.
   SVector& rowVector_w(int i)
   {
      return operator[](i);
   }

   /// returns the #rowVector of the \p i 'th #LPRow.
   const SVector& rowVector(int i) const
   {
      return operator[](i);
   }

   /// returns a writable #rowVector of the #LPRow# with #DataKey \p k.
   SVector& rowVector_w(const DataKey& k)
   {
      return operator[](k);
   }

   /// returns the #rowVector of the #LPRow# with #DataKey \p k.
   const SVector& rowVector(const DataKey& k) const
   {
      return operator[](k);
   }

   /// returns the inequalitiy type of the \p i 'th #LPRow.
   LPRow::Type type(int i) const
   {
      if (rhs(i) >= LPRow::infinity)
         return LPRow::GREATER_EQUAL;
      if (lhs(i) <= -LPRow::infinity)
         return LPRow::LESS_EQUAL;
      if (lhs(i) == rhs(i))
         return LPRow::EQUAL;
      return LPRow::RANGE;
   }

   /// returns the inequality type of the #LPRow with #DataKey \p k.
   LPRow::Type type(const DataKey& k) const
   {
      return type(number(k));
   }

   /// changes the inequality type of row \p i to \p type.
   void setType(int i, LPRow::Type type);

   /// returns the value of the \p i'th #LPRow.
   double value(int i) const
   {
      if (rhs(i) < LPRow::infinity)
         return rhs(i);
      else
      {
         assert(lhs(i) > -LPRow::infinity);
         return lhs(i);
      }
   }

   /// returns the value of the #LPRow with #DataKey \p k.
   /** The \em value of a row depends on its type: if the inequality is of
       type "greater or equal", the value is the #lhs of the row. Otherwise,
       the value is the #rhs.
   */
   double value(const DataKey& k) const
   {
      return value(number(k));
   }

   /// returns the #DataKey of the \p i 'th #LPRow in #LPRowSet.
   DataKey key(int i) const
   {
      return SVSet::key(i);
   }

   /// returns the number of the #LPRow# with #DataKey \p k in #LPRowSet.
   int number(const DataKey& k) const
   {
      return SVSet::number(k);
   }

   /// does #DataKey \p k belong to #LPRowSet ?
   int has(const DataKey& k) const
   {
      return SVSet::has(k);
   }
   //@}


   /**@name Extension
      Extension methods come with two signatures, one of them providing a
      parameter to return the assigned #DataKey(s). See #DataSet for a more
      detailed description. All extension methods will automatically rearrange
      or allocate more memory if required.
   */
   //@{
   ///
   void add(const LPRow& row)
   {
      DataKey k;
      add(k, row);
   }
   /// adds \p row to #LPRowSet.
   void add(DataKey& pkey, const LPRow& prow)
   {
      add(pkey, prow.lhs(), prow.rowVector(), prow.rhs());
   }

   ///
   void add(double plhs, const SVector& prowVector, double prhs)
   {
      DataKey k;
      add(k, plhs, prowVector, prhs);
   }
   /// adds #LPRow consisting of left hand side \p lhs, row vector \p rowVector, and right hand side \p rhs to #LPRowSet.
   void add(DataKey& key, double lhs, const SVector& rowVector, double rhs);

   ///
   void add(const LPRowSet& set);
   /// adds all #LPRow%s of \p set to #LPRowSet.
   void add(DataKey key[], const LPRowSet& set);

   /// extends row \p n to fit \p newmax nonzeros.
   void xtend(int n, int newmax)
   {
      SVSet::xtend(rowVector_w(n), newmax);
   }

   /// extend row with #DataKey \p key to fit \p newmax nonzeros.
   void xtend(const DataKey& pkey, int pnewmax)
   {
      SVSet::xtend(rowVector_w(pkey), pnewmax);
   }

   /// adds \p n nonzero (\p idx, \p val)-pairs to #rowVector with #DataKey \p k.
   void add2(const DataKey& k, int n, int idx[], double val[])
   {
      SVSet::add2(rowVector_w(k), n, idx, val);
   }

   /// adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th #rowVector.
   void add2(int i, int n, int idx[], double val[])
   {
      SVSet::add2(rowVector_w(i), n, idx, val);
   }

   ///
   SVector& create(int pnonzeros = 0, double plhs = 0, double prhs = 1)
   {
      DataKey k;
      return create(k, pnonzeros, plhs, prhs);
   }
   /// creates new #LPRow with specified parameters and returns a reference to its row vector.
   SVector& create(DataKey& nkey, int nonzeros = 0, double lhs = 0, double rhs = 1);
   //@}


   /**@name Shrinking
       See \Ref{DataSet} for a description of the renumbering of the remaining
       #LPRow#s in a #LPRowSet# after the call of a removal method.
    */
   //@{
   /// removes \p i 'th #LPRow.
   void remove(int i);
   /// removes #LPRow with #DataKey \p k.
   void remove(const DataKey& k)
   {
      remove(number(k));
   }


   /// removes multiple #LPRow%s.
   void remove(int perm[]);

   /// removes \p n #LPRow%s with #DataKey%s \p keys.
   void remove(DataKey keys[], int n)
   {
      DataArray<int> perm(num());
      remove(keys, n, perm.get_ptr());
   }

   /// removes \p n #LPRow%s with row numbers given by \p nums.
   void remove(int nums[], int n)
   {
      DataArray<int> perm(num());
      remove(nums, n, perm.get_ptr());
   }

   /// removes \p n #LPRow%s with #DataKey%s \p keys, and stores permutation of row indices in \p perm.
   void remove(DataKey keys[], int n, int* perm);

   /// removes \p n #LPRow%s with row numbers given by \p nums, and stores permutation of row indices in \p perm.
   void remove(int nums[], int n, int* perm);

   /// removes all #LPRow%s.
   void clear();
   //@}


   /**@name Memory Management
       For a description of the memory management methods, see the
       documentation of #SVSet, which has been used for implementating
       #LPRowSet.
    */
   //@{
   /// reallocates memory to be able to store \newmax #LPRow%s.
   void reMax(int newmax = 0)
   {
      SVSet::reMax(newmax);
      left.reSize (max());
      right.reSize(max());
   }

   /// returns number of used nonzero entries.
   int memSize() const
   {
      return SVSet::memSize();
   }

   /// returns length of nonzero memory.
   int memMax() const
   {
      return SVSet::memMax();
   }

   /// reallocates memory to be able to store \newmax nonzeros.
   void memRemax(int newmax)
   {
      SVSet::memRemax(newmax);
   }

   /// garbage collection in nonzero memory.
   void memPack()
   {
      SVSet::memPack();
   }
   //@}

   /// check consistency.
   bool isConsistent() const;

   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   /** The user can specify the initial maximum number of rows \p max
       and the initial maximum number of nonzero entries \p memmax. If these
       parameters are omitted, a default size is used. However, one can add
       an arbitrary number of rows to the #LPRowSet, which may result in
       automated memory realllocation.
   */
   LPRowSet(int pmax = -1, int pmemmax = -1)
      : SVSet(pmax, pmemmax), left(0), right(0)
   { }

   /// assignment operator.
   LPRowSet& operator=(const LPRowSet& rs)
   {
      if (this != &rs)
      {
         SVSet::operator=(rs);
         left  = rs.left;
         right = rs.right;
      }
      return *this;
   }
   //@}
};
} // namespace soplex
#endif // _LPROWSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
