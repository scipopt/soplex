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
#pragma ident "@(#) $Id: svset.h,v 1.3 2001/11/07 17:31:25 bzfbleya Exp $"


/*      \Section{Imports}
 */
#ifndef _SVSET_H_
#define _SVSET_H_


#include <assert.h>


#include "svector.h"
#include "dataset.h"
#include "dataarray.h"
#include "idlist.h"

namespace soplex
{

/*      \Section{Class Declarataion}
 */

/*      \SubSection{Double Linked SVector --- implementation class}
    This class is used for implementation purposes only. It should hence be
    better made a local class of #SVSet#. However, CRI T3D's C++ compiler does
    not support templete instantiation with subclasses. This is why we define it
    here and add a typedef in #SVSet#'s class scope.
 */
class SVSet_DLPSV : public SVector
{
   SVSet_DLPSV *thenext;
   SVSet_DLPSV *theprev;
public:
   SVSet_DLPSV*& next()
   {
      return thenext;
   }
   SVSet_DLPSV*const& next() const
   {
      return thenext;
   }
   SVSet_DLPSV*const& prev() const
   {
      return theprev;
   }
   SVSet_DLPSV*& prev()
   {
      return theprev;
   }
   SVector& svector()
   {
      return static_cast<SVector&>(*this);
   }
   SVSet_DLPSV()
   {}

   SVSet_DLPSV(const SVSet_DLPSV& copy) : SVector(copy)
   {}

};

typedef DataSet < SVSet_DLPSV > ::Key SVSet_Key;
typedef DataArray < SVector_Element > SVSet_Base;

/** sparse vector set.
Class #SVSet# provides a set of sparse vectors #SVector#. All #SVector#s in a
#SVSet# share one big memory block for their nonzeros. This memory is reffered
to as the {\em nonzero memory}. The #SVector#s themselfs are saved in another
memory block reffered to as the {\em vector memory}. Both memory blocks will
grow automatically if required, when adding more #SVector#s to the set or
enlarging #SVector#s within the set. For controlling memory consumption, methods
are provided to inquire and reset the size of the memory blocks used for vectors
and nonzeros.
 
#SVector#s in an #SVSet# are numbered from 0 thru #num()#-1. They can be
accessed using the index #operator[]#. When removing #SVector#s of a #SVSet#
the remaining ones will be renumbered. However, all #SVector# with a smaller
number than the lowest number of the removed #SVector#s will remain unchanged.
 
For providing a uniform access to #SVector#s in a set even if others are removed
or added, #SVSet# assigns a #Key# to each #SVector# in the set. Such a #Key#
remains unchanged as long as the corresponding #SVector# is in the #SVSet#, no
matter what other #SVector#s are added to or removed from the #SVSet#. Methods
are provided for getting the #Key# to a #SVector# or its number and vice versa.
Further, each #add()# method for enlarging an #SVSet# is provided with two
signatures. One of them returns the #Key#s assigned to the #SVector#s added to
the #SVSet#.
 */
class SVSet : protected SVSet_Base
{
   /*  The management of Keys is left for #DataSet#
    */
   typedef SVSet_DLPSV DLPSV;
   DataSet < SVSet_DLPSV > set;

   /*  while the management of nonzeros is done in this class. It requires a
       double linked #list# of #DLPSV#s, where the #SVector#s are kept in the
       order their indices occurr in the inherited #DataArray#. It keeps the
       #SVector#s without holes: If one is removed or moved to the end, the
       one preceeding it optains all the nonzeros that previously belonged to
       the (re-)moved one.  However, the nonzeros in use are uneffected by
       this.
    */
   IdList < SVSet_DLPSV > list;

   /*  Make sure to provide enough vector memory for #n# more #SVector#s.
    */
   void ensurePSVec(int n)
   {
      if (num() + n > max())
      {
         assert(factor > 1);
         reMax(int(factor*max()) + 8 + n);
      }
   }

   /*  Make sure to provide enough nonzero memory for #n# more #Elements#s.
    */
   void ensureMem(int n);

public:
   // typedef DataSet<DLPSV>::Key      Key;
   typedef SVSet_Key Key;

   /**@name Control Parameters */
   //@{
   /** Sparse vector memory enlargment factor.
       If the #SVSet# runs out of vector memory, it is enlareged by
       #factor#.
    */
   double factor;

   /** Nonzero element memory enlargment factor.
       If the #SVSet# runs out of nonzero memory it is enlareged by a
       #memFactor#.
    */
   double& memFactor;
   //@}


   /**@name Extension */
   //@{
   /** Add #SVector svec# to the set.
    *  This includes copying its nonzeros to the sets nonzero memory and
    *  creating an additional #SVector# entry in vector memory. If
    *  neccessary, the memory blocks are enlarged appropriately.
    */
   void add(const SVector& svec)
   {
      ensurePSVec(1);
      SVector* new_svec = create(svec.size());
      *new_svec = svec;
   }

   /** Add #svec# to #SVSet#.
    *  Adds #SVector svec# to the set. This includes copying its nonzeros
    *  to the sets nonzero memory and creating an additional #SVector#
    *  entry in vector memory. If neccessary, the memory blocks are
    *  enlarged appropriately. Upon return #nkey# contains the #Key#, that
    *  the #SVSet# has assosicated to the new #SVector#.
    */
   void add(Key& nkey, const SVector& svec)
   {
      ensurePSVec(1);
      SVector* new_svec = create(nkey, svec.size());
      *new_svec = svec;
   }

   /// Add all #n# #SVector#s in the array #svec# to the set.
   void add(const SVector svec[], int n);

   /** Add #n# #SVector#s to #SVSet#.
    *  Adds all #n# #SVector#s in the array #svec# to the set.  Upon return
    *  #nkey# contains the #Key#s, that the #SVSet# has assosicated to the
    *  new #SVector#s. Hence, array #nkey# must be large enough to fit #n#
    *  #Key#s.
    */
   void add(Key nkey[], const SVector svec[], int n);

   /// Add all #SVector#s in #set# to an #SVSet#.
   void add(const SVSet& set);

   /** Add all #SVector#s of #set# to #SVSet#.
    *  Adds all #n# #SVector#s in the #set# to an #SVSet#. Upon return
    *  #nkey# contains the #Key#s, that the #SVSet# has assosicated to the
    *  new #SVector#s. Hence, array #nkey# must be large enough to fit
    *  #set.num()# #Key#s.
    */
   void add(Key nkey[], const SVSet& set);

   /** Creates new #SVector# in set.
    *  The new #SVector# will be ready to fit at least #idxmax# nonzeros.
    */
   SVector* create(int idxmax = -1);

   /** Creates new #SVector# in set.
    *  The new #SVector# will be ready to fit at least #idxmax# nonzeros.
    *  Upon return #nkey# contains the #Key# associated to the new
    *  #SVector#.
    */
   SVector* create(Key& nkey, int idxmax = -1);

   /** Extend #svec# to fit #newmax# nonzeros.
       It is an error, if #svec# is not an #SVector# of the #SVSet#.
    */
   void xtend(SVector& svec, int newmax);

   /** Add nonzero (#idx#, #val#) to #svec# of this #SVSet#.
    *  Adds one nonzero (#idx#, #val#) to #SVector svec# in the #SVSet#. It
    *  is an error, if #svec# is not an #SVector# of the #SVSet#. If #svec#
    *  is not large enough to hold the additional nonzero, it will be
    *  automatically enlarged within the set.
    */
   void add2(SVector &svec, int idx, double val);

   /** Add #n# nonzeros to #svec# of this #SVSet#.
    *  Adds #n# nonzero to #SVector svec# in the #SVSet#. It is an error,
    *  if #svec# is not an #SVector# of the #SVSet#. If #svec# is not large
    *  enough to hold the additional nonzeros, it will be automatically
    *  enlarged within the set.
    */

   void add2(SVector &svec, int n, const int idx[], const double val[]);
   //@}


   /**@name Shrinking */
   //@{
   ///
   void remove(Key removekey);

   ///
   void remove(int removenum)
   {
      remove(key(removenum));
   }

   /// remove one #SVector# from set.
   void remove(SVector *svec)
   {
      remove(key(svec));
   }

   /** remove multiple elements.
       The following method removes all #Svector#s for the #SVSet# with an
       index #i# such that #perm[i] < 0#. Upon completion, #perm[i] >= 0#
       indicates the new index where the #i#-th #SVector# has been moved to
       due to this removal. Note, that #perm# must point to an array of at
       least #num()# #int#s.
    */
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

   /** Remove #n# #SVector#s from set.
    *  Removing multiple #SVector#s with one method invocation is available
    *  two flavours. An array #perm# can be passed as third argument or
    *  left not, which must be at least of size #num()# and returns the
    *  permutations resulting from this removal: #perm[i] < 0# indicates,
    *  that the element to index #i# has been removed. Otherwise, #perm[i]#
    *  is the new index of the element with index #i# before the removal.
    */
   void remove(int nums[], int n, int* perm);

   /// remove all #SVector#s from set.
   void clear()
   {
      DataArray < SVector::Element > ::clear();
      DataArray < SVector::Element > ::reMax(10000);
      set.clear();
      list.clear();
   }
   //@}


   /**@name Access */
   //@{
   ///
   SVector& operator[](int n)
   {
      return set[n];
   }

   ///
   const SVector& operator[](int n) const
   {
      return set[n];
   }

   ///
   SVector& operator[](const Key& k)
   {
      return set[k];
   }

   ///
   const SVector& operator[](const Key& k) const
   {
      return set[k];
   }
   //@}


   /**@name Inquiry */
   //@{
   /// current nr. of #SVector#s.
   int num() const
   {
      return set.num();
   }

   /// current maximum nr. of #SVector#s.
   int max() const
   {
      return set.max();
   }

   ///
   Key key(int n) const
   {
      return set.key(n);
   }

   ///
   Key key(const SVector* svec) const
   {
      return set.key(static_cast<const DLPSV*>(svec));
   }

   ///
   int number(const Key& k) const
   {
      return set.number(k);
   }

   ///
   int number(const SVector* svec) const
   {
      return set.number(static_cast<const DLPSV*>(svec));
   }


   ///
   int has(const Key& k) const
   {
      return set.has(k);
   }

   ///
   int has(int n) const
   {
      return set.has(n);
   }

   /// is an #SVector# in the set.
   int has(const SVector* svec) const
   {
      return set.has(static_cast<const DLPSV*>(svec));
   }
   //@}


   /**@name Memory Management */
   //@{
   /// used nonzero memory.
   int memSize() const
   {
      return DataArray < SVector::Element > ::size();
   }

   /// length of nonzero memory.
   int memMax() const
   {
      return DataArray < SVector::Element > ::max();
   }

   /// reset length of nonzero memory.
   void memRemax(int newmax);

   /// garbage collection in nonzero memory.
   void memPack();
   //@}


   /**@name Miscellaneous */
   //@{
   /// reset maximum number of #SVector#s.
   void reMax(int newmax = 0);

   /// consistency check.
   int isConsistent() const;

   /// assignment operator.
   SVSet& operator=(const SVSet& rhs);

   /// copy constructor.
   SVSet(const SVSet& old);

   /// default constructor.
   SVSet(int pmax = -1,
         int pmemmax = -1,
         double pfac = 1.1,
         double pmemFac = 1.2)
      : DataArray < SVector::Element >
         (0, (pmemmax > 0) ? pmemmax : 8 * ((pmax > 0) ? pmax : 8), pmemFac)
         , set ((pmax > 0) ? pmax : 8)
         , factor (pfac)
         , memFactor (DataArray < SVector::Element > ::memFactor)
   { }
   //@}
};

} // namespace soplex
#endif // _SVSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
