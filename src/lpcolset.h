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
#pragma ident "@(#) $Id: lpcolset.h,v 1.8 2001/12/28 14:55:12 bzfkocht Exp $"

/**@file  lpcolset.h
 * @brief Set of LP columns.
 */
#ifndef _LPCOLSET_H_
#define _LPCOLSET_H_

#include <assert.h>

#include "lpcol.h"
#include "dvector.h"
#include "svset.h"
#include "datakey.h"

namespace soplex
{
/**@brief   Set of LP columns.
   @ingroup Algebra

   Class #LPColSet implements a set of #LPCol%s. Unless for memory
   limitations, any number of #LPCol%s may be #add%ed to an #LPColSet. Single
   or multiple #LPCol%s may be #add%ed to an #LPColSet, where each method
   #add() comes with two different signatures. One with and one without a
   parameter, used for returning the #Key%s assigned to the new #LPCol%s
   by the set. See #DataSet::Key for a more detailed description of the
   concept of #Key%s. For the concept of renumbering #LPCol%s within an
   #LPColSet after removal of some #LPCol%s see #DataSet.
   
   @see        DataSet, DataSet::Key
 */
class LPColSet : protected SVSet
{
private:
   /**@todo  get rid of this typedef! */
   typedef DataKey Key;

   DVector low;     ///< vector of lower bounds.
   DVector up;      ///< vector of upper bounds.
   DVector object;  ///< vector of objective coefficients.

public:

   /**@name Inquiry */
   //@{
   /// returns the number of #LPCol%s currently in #LPColSet.
   int num() const
   {
      return SVSet::num();
   }

   /// returns maximum number of #LPCol%s currently fitting into #LPColSet.
   int max() const
   {
      return SVSet::max();
   }

   ///
   const Vector& obj() const
   {
      return object;
   }
   /// returns vector of #obj%ective values.
   Vector& obj()
   {
      return object;
   }

   ///
   double obj(int i) const
   {
      return object[i];
   }
   /// returns #obj%jective value of \p i 'th #LPCol in #LPColSet.
   double& obj(int i)
   {
      return object[i];
   }

   ///
   double obj(const Key& k) const
   {
      return object[number(k)];
   }
   /// returns #obj%jective value of #LPCol with #Key \p k in #LPColSet.
   double& obj(const Key& k)
   {
      return object[number(k)];
   }

   ///
   const Vector& lower() const
   {
      return low;
   }
   /// returns vector of #lower bound values.
   Vector& lower()
   {
      return low;
   }

   ///
   double lower(int i) const
   {
      return low[i];
   }
   /// returns #lower bound of \p i 'th #LPCol in #LPColSet.
   double& lower(int i)
   {
      return low[i];
   }

   ///
   double lower(const Key& k) const
   {
      return low[number(k)];
   }
   /// returns #lower bound of #LPCol# with #Key \p k in #LPColSet.
   double& lower(const Key& k)
   {
      return low[number(k)];
   }

   ///
   const Vector& upper() const
   {
      return up;
   }
   /// returns vector of #upper bound values.
   Vector& upper()
   {
      return up;
   }

   ///
   double upper(int i) const
   {
      return up[i];
   }
   /// returns #upper bound of \p i 'th #LPCol in #LPColSet.
   double& upper(int i)
   {
      return up[i];
   }

   ///
   double upper(const Key& k) const
   {
      return up[number(k)];
   }
   /// returns #upper bound of #LPCol with #Key \p k in #LPColSet.
   double& upper(const Key& k)
   {
      return up[number(k)];
   }

   ///
   SVector& colVector_w(int i)
   {
      return operator[](i);
   }
   /// returns #colVector of \p i 'th #LPCol in #LPColSet.
   const SVector& colVector(int i) const
   {
      return operator[](i);
   }

   /// returns writeable #colVector of #LPCol with #Key \p k in #LPColSet.
   SVector& colVector_w(const Key& k)
   {
      return operator[](k);
   }

   /// returns #colVector of #LPCol with #Key \p k in #LPColSet.
   const SVector& colVector(const Key& k) const
   {
      return operator[](k);
   }

   /// returns #Key of \p i 'th #LPCol in #LPColSet.
   Key key(int i) const
   {
      return SVSet::key(i);
   }

   /// returns number of #LPCol# with #Key \p k in #LPColSet.
   int number(const Key& k) const
   {
      return SVSet::number(k);
   }

   /// does #Key \p k belong to #LPColSet ?
   int has(const Key& k) const
   {
      return SVSet::has(k);
   }
   //@}


   /**@name Extension
      All extension methods come with two signatures, one of which providing a
      parameter to return the assigned #Key%(s). See #DataSet for a more
      detailed description. All extension methods are designed to
      automatically reallocate memory if required.
   */
   //@{
   ///
   void add(const LPCol& pcol)
   {
      Key k;
      add(k, pcol);
   }
   /// adds p pcol to #LPColSet.
   void add(Key& pkey, const LPCol& pcol)
   {
      add(pkey, pcol.obj(), pcol.lower(),
           pcol.colVector(), pcol.upper());
   }

   ///
   void add(double pobj, double plower, const SVector& pcolVector, double pupper)
   {
      Key k;
      add(k, pobj, plower, pcolVector, pupper);
   }
   /// adds #LPCol consisting of objective value \p obj, lower bound \p lower, column vector \p colVector and upper bound \p upper to #LPColSet.
   void add (Key& key,
             double obj,
             double lower,
             const SVector& colVector,
             double upper);

   ///
   void add(const LPColSet& set);
   /// adds all #LPCol%s of \p set to #LPColSet.
   void add(Key key[], const LPColSet& set);

   ///
   void add2(const Key& k, int n, int idx[], double val[])
   {
      SVSet::add2(colVector_w(k), n, idx, val);
   }
   /// adds \p n nonzero (\p idx, \p val)-pairs to \p i 'th #colVector.
   void add2(int i, int n, int idx[], double val[])
   {
      SVSet::add2(colVector_w(i), n, idx, val);
   }

   ///
   SVector& create(int pnonzeros = 0, double pobj = 1, double plw = 0, double pupp = 1)
   {
      Key k;
      return create(k, pnonzeros, pobj, plw, pupp);
   }
   /// creates new #LPCol with specified arguments and returns a reference to its column vector.
   SVector& create(Key& nkey, int nonzeros = 0, double obj = 1, double low = 0, double up = 1);
   //@}


   /**@name Shrinking
      See #DataSet for a description of the renumbering of the remaining
      #LPCol%s in a #LPColSet after the call of a removal method.
   */
   //@{
   /// removes \p i 'th #LPCol.
   void remove(int i);

   /// removes #LPCol with #Key \p k.
   void remove(const Key& k)
   {
      remove(number(k));
   }

   /// removes multiple elements.
   void remove(int perm[]);

   /// removes all #LPCol%s associated to the \p n #Key%s \p keys.
   void remove(Key keys[], int n)
   {
      DataArray < int > perm(num());
      remove(keys, n, perm.get_ptr());
   }

   /// removes #LPCol%s with numbers \p nums, where \p n is the length of the array \p nums
   void remove(int nums[], int n)
   {
      DataArray < int > perm(num());
      remove(nums, n, perm.get_ptr());
   }

   /// removes all #LPCol%s associated to the \p n #Key%s \p keys, and stores the index permutation in array \p perm.
   void remove(Key keys[], int n, int* perm);

   /// removes #LPCol%s with numbers \p nums, where \p n is the length of the array \p nums, and stores the index permutation in array \p perm.
   void remove(int nums[], int n, int* perm);

   /// removes all #LPCol%s from the set.
   void clear();
   //@}


   /**@name Memory Management
      See #SVSet for a description of the memory management methods.
   */
   //@{
   /// reallocates memory to be able to store \p newmax #LPCol%s.
   void reMax(int newmax = 0)
   {
      SVSet::reMax (newmax);
      up.reSize (max());
      low.reSize (max());
      object.reSize(max());
   }

   /// returns used nonzero memory.
   int memSize() const
   {
      return SVSet::memSize();
   }

   /// returns length of nonzero memory.
   int memMax() const
   {
      return SVSet::memMax();
   }

   /// resets length of nonzero memory.
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
   int isConsistent() const;
   //@}

   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   /** The user can specify the initial maximum number of columns \p max
       and the initial maximum number of nonzero entries \p memmax. If these
       parameters are omitted, a default size is used. However, one can add
       an arbitrary number of columns to the #LPColSet, which may result in
       automated memory realllocation.
   */
   LPColSet(int pmax = -1, int pmemmax = -1)
      : SVSet(pmax, pmemmax), low(0), up(0), object(0)
   { }

   /// assignment operator.
   LPColSet& operator=(const LPColSet& rs)
   {
      if (this != &rs)
      {
         SVSet::operator=(rs);
         low = rs.low;
         up = rs.up;
         object = rs.object;
      }
      return *this;
   }
   //@}

};


} // namespace soplex
#endif // _LPCOLSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
