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
#pragma ident "@(#) $Id: lpcolset.h,v 1.3 2001/11/07 17:31:18 bzfbleya Exp $"

#ifndef _LPCOLSET_H_
#define _LPCOLSET_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */
#include "lpcol.h"
#include "dvector.h"
#include "svset.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** Set of LP columns.
    Class #LPColSet# implements a set of #LPCol#s. Unless for memory
    limitations, any number of #LPCol#s may be #add#ed to an #LPColSet#. Single
    or multiple #LPCol#s may be #add#ed to an #LPColSet#, where each method
    #add()# comes with two different signatures. One with an one without a
    parameter, used for returning the #Key#s assigned to the new #LPCol#s
    by the set. See \Ref{DataSet::Key} for a more detailed description of the
    concept of #Key#s. For the concept of renumbering #LPCol#s within an
    #LPColSet# after removal of some #LPCol#s see \Ref{DataSet}.
 
    @see        DataSet, DataSet::Key
 */
class LPColSet : protected SVSet
{
private:
   DVector low, up, object;

public:
   ///
   typedef SVSet::Key Key;

   /**@name Inquiry */
   //@{
   /// Return maximum number of #LPCol#s currently in #LPColSet#.
   int num() const
   {
      return SVSet::num();
   }
   /// Return maximum number of #LPCol#s currently fitting into #LPColSet#.
   int max() const
   {
      return SVSet::max();
   }

   ///
   const Vector& obj() const
   {
      return object;
   }
   /// return vector of #obj# values.
   Vector& obj()
   {
      return object;
   }

   ///
   double obj(int i) const
   {
      return object[i];
   }
   /// return #obj# of #i#-th #LPCol# in #LPColSet#.
   double& obj(int i)
   {
      return object[i];
   }

   ///
   double obj(const Key& k) const
   {
      return object[number(k)];
   }
   /// return #obj# of #k#-th #LPCol# in #LPColSet#.
   double& obj(const Key& k)
   {
      return object[number(k)];
   }

   ///
   const Vector& lower() const
   {
      return low;
   }
   /// return vector of #lower# values.
   Vector& lower()
   {
      return low;
   }

   ///
   double lower(int i) const
   {
      return low[i];
   }
   /// return #lower# of #i#-th #LPCol# in #LPColSet#.
   double& lower(int i)
   {
      return low[i];
   }

   ///
   double lower(const Key& k) const
   {
      return low[number(k)];
   }
   /// return #lower# of #k#-th #LPCol# in #LPColSet#.
   double& lower(const Key& k)
   {
      return low[number(k)];
   }

   ///
   const Vector& upper() const
   {
      return up;
   }
   /// return vector of #upper# values.
   Vector& upper()
   {
      return up;
   }

   ///
   double upper(int i) const
   {
      return up[i];
   }
   /// return #upper# of #i#-th #LPCol# in #LPColSet#.
   double& upper(int i)
   {
      return up[i];
   }

   ///
   double upper(const Key& k) const
   {
      return up[number(k)];
   }
   /// return #upper# of #k#-th #LPCol# in #LPColSet#.
   double& upper(const Key& k)
   {
      return up[number(k)];
   }

   ///
   SVector& colVector(int i)
   {
      return operator[](i);
   }
   /// return #colVector# of #i#-th #LPCol# in #LPColSet#.
   const SVector& colVector(int i) const
   {
      return operator[](i);
   }

   ///
   const SVector& colVector(const Key& k) const
   {
      return operator[](k);
   }
   /// return #colVector# of #k#-th #LPCol# in #LPColSet#.
   SVector& colVector(const Key& k)
   {
      return operator[](k);
   }

   /// return #Key# of #i#-th #LPCol# in #LPColSet#.
   Key key(int i) const
   {
      return SVSet::key(i);
   }
   /// return number of #k#-th #LPCol# in #LPColSet#.
   int number(const Key& k) const
   {
      return SVSet::number(k);
   }
   /// does #Key k# belong to #LPColSet#?.
   int has(const Key& k) const
   {
      return SVSet::has(k);
   }
   //@}


   /**@name Extension
       All extension methods come with two signatures, one of which providing a
       parameter to return the assigned #Key#(s). See \Ref{DataSet} for a more
       detailed description. All extension methods are designed to
       automatically realloc memory if required.
    */
   //@{
   ///
   void add(const LPCol& col)
   {
      Key k;
      add(k, col);
   }
   /// add #col# to #LPColSet#.
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
   /** add #LPCol# consisting of #lower#, #colVector# and #upper# to
       #LPColSet#
    */
   void add (Key& key,
             double obj,
             double lower,
             const SVector& colVector,
             double upper);

   ///
   void add(const LPColSet& set);
   /// add all #LPCol#s of #set# to #LPColSet#.
   void add(Key key[], const LPColSet& set);

   ///
   void add2(const Key& k, int n, int idx[], double val[])
   {
      SVSet::add2(colVector(k), n, idx, val);
   }
   /// add #n# nonzero (#idx#, #val#) to #i#-th #colVector#..
   void add2(int i, int n, int idx[], double val[])
   {
      SVSet::add2(colVector(i), n, idx, val);
   }

   ///
   SVector& create(int pnonzeros = 0, double pobj = 1, double plw = 0, double pupp = 1)
   {
      Key k;
      return create(k, pnonzeros, pobj, plw, pupp);
   }
   /** Create new #LPCol# with specified arguments and return a reference
       to its column vector.
    */
   SVector& create(Key& nkey, int nonzeros = 0, double obj = 1, double low = 0, double up = 1);
   //@}


   /**@name Shrinking
       See \Ref{DataSet} for a description of the renumbering of the remaining
       #LPCol#s in a #LPColSet# after the call of a removal method.
    */
   //@{
   /// remove #i#-th #LPCol#.
   void remove(int i);
   /// remove #k#-th #LPCol#.
   void remove(const Key& k)
   {
      remove(number(k));
   }


   /// remove multiple elements.
   void remove(int perm[]);

   ///
   void remove(Key keys[], int n)
   {
      DataArray < int > perm(num());
      remove(keys, n, perm.get_ptr());
   }

   ///
   void remove(int nums[], int n)
   {
      DataArray < int > perm(num());
      remove(nums, n, perm.get_ptr());
   }

   ///
   void remove(Key keys[], int n, int* perm);

   /// remove #n# #LPCol#s from set.
   void remove(int nums[], int n, int* perm);

   /// remove all #LPCol#s.
   void clear();
   //@}


   /**@name Memory Management
    *  See \Ref{SVSet} for a description of the memory management methods.
    */
   //@{
   ///
   void reMax(int newmax = 0)
   {
      SVSet::reMax (newmax);
      up.reSize (max());
      low.reSize (max());
      object.reSize(max());
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


   ///
   LPColSet& operator=(const LPColSet& rs)
   {
      SVSet::operator=(rs);
      low = rs.low;
      up = rs.up;
      object = rs.object;
      return *this;
   }

   ///
   LPColSet(int pmax = -1, int pmemmax = -1)
      : SVSet(pmax, pmemmax), low(0), up(0), object(0)
   { }
   //@}

   ///

   int isConsistent() const;
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
