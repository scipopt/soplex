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
#pragma ident "@(#) $Id: idxset.h,v 1.3 2001/11/11 20:27:30 bzfkocht Exp $"


/*      \Section{Imports}
 */
#ifndef _IDXSET_H_
#define _IDXSET_H_

#include <assert.h>

namespace soplex
{


/*      \Section{Class Declaration}
 */
/** index set.
    Class #IdxSet# provides a set of indices. At construction it must be given
    an array of #int# where to store the indice and its length. The array will
    from then on be managed by the #IdxSet#.
 
    Indices are implicitely numbered from 0 thru #size()-1#. They can be
    accessed (and altered) via method #index()# with the desired index number as
    argument.  Range checking is performed in the debug version.
 
    Indices may be added or removed from the set, by calling #add()# or
    #remove()# methods, respectively. However, no #IdxSet# can hold more then
    #max()# indices, i.e. the number given at the constructor.
 
    When removing indices, the remaining ones are renumbered. However, all
    indices before the first removed index keep their number unchanged.
*/
class IdxSet
{
   /*      \SubSection{Datastructures}
       Die intern Darstellung eines #IdxSet#s besteht i.w.~aus einem Feld f\"ur
       die Indizes (#idx#).  In #len# wird dessen L\"ange gehalten, wovon die
       ersten #num# tats\"achlich benutzt werden. Der Feld wird von der Klasse
       nicht selbst alloziiert, sondern mu\ss\ Instanzen bei der Initialisierung
       \"ubergeben werden.
    */
protected:
   friend class Vector;

   int num;            // nr. of used indices
   int len;            // length of array #idx#
   int *idx;           // array of indices

public:
   /**@name Access */
   //@{
   ///
   int& index(int n)
   {
      assert(n >= 0 && n < size());
      return idx[n];
   }

   /// access #n#-th index.
   int index(int n) const
   {
      assert(n >= 0 && n < size());
      return idx[n];
   }
   //@}

   /**@name Inquiery */
   //@{
   /// Number of used indices.
   int size() const
   {
      return num;
   }

   /// Maximal number indices.
   int max() const
   {
      return len;
   }

   /// Maximal index.
   int dim() const;

   /** Number of index #i#.
       Return the number of the first index #i#. If no index #i# is
       available in the #IdxSet#, -1 is returned. Otherwise,
       #index(number(i)) == i# hods.
    */
   int number(int i) const;
   //@}

   /**@name Extension
       An #IdxSet# cannot be extended to fit more than #max()# elements. If
       neccessary, the user must explicitely provide the #IdxSet# with a
       suitable memory. Alternatively, one can use \Ref{DIdxSet}s which provide
       the required memory managemant.
    */
   //@{
   /// append #n# uninitialized indices.
   void add(int n)
   {
      assert(n >= 0 && n + size() <= max());
      num += n;
   }

   /// append all indices of #set#.
   void add(const IdxSet& set)
   {
      add(set.size(), set.idx);
   }

   /// Append #n# indices in #i#.
   void add(int n, const int i[]);

   /// Append index #i#.
   void addIdx(int i)
   {
      assert(size() < max());
      idx[num++] = i;
   }
   //@}

   /**@name Removal */
   //@{
   /// remove indices #n# through #m#.
   void remove(int n, int m);

   /// remove #n#-th index.
   void remove(int n)
   {
      if (n < size() && n >= 0)
         idx[n] = idx[--num];
   }

   /// remove all indices.
   void clear()
   {
      num = 0;
   }
   //@}

   /**@name Internals
       The use of the following functions is not encouraged since they
       interfere with the internal representation of #IdxSet#. Such consists of
       a pointer to the array of indices that has been passed to the
       constructor, along with an #int# for the maximal length of this array
       and the number of elements currently in use. These members may be
       changed with the following methods.
    */
   //@{
   /** reset the size of the index array.
       Handle with care to prevent memory leakages. It is not safe to
       enlarge #max()# without providing an extended index array.
    */
   void setMax(int mx)
   {
      len = mx;
   }

   /// reset the number of indices.
   void setSize(int sz)
   {
      num = sz;
      assert(size() <= max());
      assert(size() >= 0);
   }

   ///
   int*& indexMem()
   {
      return idx;
   }
   /// pointer to the index array.
   const int* indexMem() const
   {
      return idx;
   }

   ///
   operator int* ()
   {
      return idx;
   }
   /// pointer to the index array.
   operator const int* () const
   {
      return idx;
   }
   //@}

   /**@name Miscellaneous */
   //@{
   /// switch n'th and 0'th index.
   void toFront(int n);

   /** Assignment operator.
       The assignment operator copies all nonzeros of the right handside
       #IdxSet# to the left one. This implies, that the latter must have
       enough index memory.
    */
   IdxSet& operator=(const IdxSet& set);

   /// consistency check.
   int isConsistent() const;

   /** Constructor.
       The constructur receives the index memory #imem# to use for saving
       its indices. This must be large enough to fit #n# indices. #l# can
       be given to construct an #IdxSet# initialized to the #l# first
       indices in #imem#.
    */
   IdxSet(int n, int imem[], int l = 0)
      : num(l), len(n), idx(imem)
   {
      assert(isConsistent());
   }
   IdxSet()
      : num(0), len(0), idx(0)
   {
      assert(isConsistent());
   }

   /**@name Copy constructor
       The copy constructor creates an #IdxSet# that {\em shares} the same
       index array. This is fine for argument passing in function calls,
       but may be dangerous if one keeps a copy constructed #IdxSet# on
       data, that has been released.
    */
   //@}
};

} // namespace soplex
#endif // _IDXSET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
