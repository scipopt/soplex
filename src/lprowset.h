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
#pragma ident "@(#) $Id: lprowset.h,v 1.6 2001/11/21 09:30:13 bzfkocht Exp $"

#ifndef _LPROWSET_H_
#define _LPROWSET_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */
#include "lprow.h"
#include "dvector.h"
#include "svset.h"
#include "datakey.h"

namespace soplex
{

//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** set of LP rows.
    Class #LPRowSet# implements a set of #LPRow#s. Unless for memory
    limitations, any number of #LPRow#s may be #add#ed to an #LPRowSet#. Single
    or multiple #LPRow#s may be #add#ed to an #LPRowSet#, where each method
    #add()# comes with two different signatures. One with an one without a
    parameter, used for returning the #Key#s assigned to the new #LPRow#s
    by the set. See \Ref{DataSet::Key} for a more detailed description of the
    concept of #Key#s. For the concept of renumbering #LPRow#s within an
    #LPRowSet# after removal of some #LPRow#s see \Ref{DataSet}.
 */
class LPRowSet : protected SVSet
{
private:
   ///
   typedef DataKey Key;

   DVector left;
   DVector right;

public:

   /**@name Inquiry */
   //@{
   /// number of #LPRow#s in #LPRowSet#.
   int num() const
   {
      return SVSet::num();
   }
   /// maximum number of #LPRow#s that fit.
   int max() const
   {
      return SVSet::max();
   }

   ///
   const Vector& lhs() const
   {
      return left;
   }
   /// vector of #lhs# values.
   Vector& lhs()
   {
      return left;
   }

   ///
   double lhs(int i) const
   {
      return left[i];
   }
   /// #lhs# of #i#-th #LPRow#.
   double& lhs(int i)
   {
      return left[i];
   }

   ///
   double lhs(const Key& k) const
   {
      return left[number(k)];
   }
   /// #lhs# of #k#-th #LPRow# in #LPRowSet#.
   double& lhs(const Key& k)
   {
      return left[number(k)];
   }

   ///
   const Vector& rhs() const
   {
      return right;
   }
   /// vector of #rhs# values.
   Vector& rhs()
   {
      return right;
   }

   ///
   double rhs(int i) const
   {
      return right[i];
   }
   /// #rhs# of #i#-th #LPRow#.
   double& rhs(int i)
   {
      return right[i];
   }

   ///
   double rhs(const Key& k) const
   {
      return right[number(k)];
   }
   /// #rhs# of #k#-th #LPRow#.
   double& rhs(const Key& k)
   {
      return right[number(k)];
   }

   ///
   SVector& rowVector(int i)
   {
      return operator[](i);
   }
   /// #rowVector# of #i#-th #LPRow#.
   const SVector& rowVector(int i) const
   {
      return operator[](i);
   }

   ///
   const SVector& rowVector(const Key& k) const
   {
      return operator[](k);
   }
   /// #rowVector# of #k#-th #LPRow#.
   SVector& rowVector(const Key& k)
   {
      return operator[](k);
   }

   ///
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
   /// inequality type of #k#-th #LPRow#.
   LPRow::Type type(const Key& k) const
   {
      return type(number(k));
   }

   /// change type of row i
   void setType(int i, LPRow::Type type);

   ///
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
   /// value of #k#-th #LPRow#.
   double value(const Key& k) const
   {
      return value(number(k));
   }

   /// return #Key# of #i#-th #LPRow# in #LPRowSet#.
   Key key(int i) const
   {
      return SVSet::key(i);
   }
   /// return number of #k#-th #LPRow# in #LPRowSet#.
   int number(const Key& k) const
   {
      return SVSet::number(k);
   }
   /// does #Key k# belong to #LPRowSet#?.
   int has(const Key& k) const
   {
      return SVSet::has(k);
   }
   //@}


   /**@name Extension
       Extension methods come with two signatures, one of which providing a
       parameter to return the assigned #Key#(s). See \Ref{DataSet} for a more
       detailed description. All extension methods will automatically rearrange
       or allocate more memory if required.
    */
   //@{
   ///
   void add(const LPRow& row)
   {
      Key k;
      add(k, row);
   }
   /// add #row# to #LPRowSet#.
   void add(Key& pkey, const LPRow& prow)
   {
      add(pkey, prow.lhs(), prow.rowVector(), prow.rhs());
   }

   ///
   void add(double plhs, const SVector& prowVector, double prhs)
   {
      Key k;
      add(k, plhs, prowVector, prhs);
   }
   /// add #LPRow# consisting of #lhs#, #rowVector# and #rhs# to #LPRowSet#.
   void add(Key& key, double lhs, const SVector& rowVector, double rhs);

   ///
   void add(const LPRowSet& set);
   /// add all #LPRow#s of #set# to #LPRowSet#.
   void add(Key key[], const LPRowSet& set);

   /// extend row #n# to fit #newmax# nonzeros.
   void xtend(int n, int newmax)
   {
      SVSet::xtend(rowVector(n), newmax);
   }
   /// extend row #key# to fit #newmax# nonzeros.
   void xtend(const Key& pkey, int pnewmax)
   {
      SVSet::xtend(rowVector(pkey), pnewmax);
   }
   ///
   void add2(const Key& k, int n, int idx[], double val[])
   {
      SVSet::add2(rowVector(k), n, idx, val);
   }
   /// add #n# nonzero (#idx#, #val#) to #i#-th #rowVector#..
   void add2(int i, int n, int idx[], double val[])
   {
      SVSet::add2(rowVector(i), n, idx, val);
   }

   ///
   SVector& create(int pnonzeros = 0, double plhs = 0, double prhs = 1)
   {
      Key k;
      return create(k, pnonzeros, plhs, prhs);
   }
   /** Create new #LPRow# with specified parameters and return a reference
       to its row vector.
    */
   SVector& create(Key& nkey, int nonzeros = 0, double lhs = 0, double rhs = 1);
   //@}


   /**@name Shrinking
       See \Ref{DataSet} for a description of the renumbering of the remaining
       #LPRow#s in a #LPRowSet# after the call of a removal method.
    */
   //@{
   /// remove #i#-th #LPRow#.
   void remove(int i);
   /// remove #k#-th #LPRow#.
   void remove(const Key& k)
   {
      remove(number(k));
   }


   /// remove multiple elements.
   void remove(int perm[]);

   ///
   void remove(Key keys[], int n)
   {
      DataArray<int> perm(num());
      remove(keys, n, perm.get_ptr());
   }

   ///
   void remove(int nums[], int n)
   {
      DataArray<int> perm(num());
      remove(nums, n, perm.get_ptr());
   }

   ///
   void remove(Key keys[], int n, int* perm);

   /// remove #n# #LPRow#s.
   void remove(int nums[], int n, int* perm);

   /// remove all #LPRow#s.
   void clear();
   //@}


   /**@name Memory Management
       For a description of the memory management methods, see the
       documentation of #SVSet#, which has benn used for implementating
       #LPRowSet#.
    */
   //@{
   ///
   void reMax(int newmax = 0)
   {
      SVSet::reMax(newmax);
      left.reSize (max());
      right.reSize(max());
   }
   /// used nonzero memory.
   int memSize() const
   {
      return SVSet::memSize();
   }

   /// length of nonzero memory.
   int memMax() const
   {
      return SVSet::memMax();
   }

   /// reset length of nonzero memory.
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


   /**@name Miscellaneous */
   //@{
   ///
   LPRowSet& operator=(const LPRowSet& rs)
   {
      SVSet::operator=(rs);
      left = rs.left;
      right = rs.right;
      return *this;
   }

   ///
   LPRowSet(int pmax = -1, int pmemmax = -1)
      : SVSet(pmax, pmemmax), left(0), right(0)
   { }

   /// check consistency.

   int isConsistent() const;
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
