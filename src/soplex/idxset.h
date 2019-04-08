/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  idxset.h
 * @brief Set of indices.
 */
#ifndef _IDXSET_H_
#define _IDXSET_H_

#include "soplex/spxdefines.h"
#include "soplex/spxalloc.h"
#include <assert.h>

namespace soplex
{
/**@brief   Set of indices.
   @ingroup Elementary

   Class IdxSet provides a set of indices. At construction it must be given
   an array of int where to store the indice and its length. The array will
   from then on be managed by the IdxSet.

   Indices are implicitely numbered from 0 thru size()-1. They can be
   accessed (and altered) via method index() with the desired index number as
   argument.  Range checking is performed in the debug version.

   Indices may be added or removed from the set, by calling add() or
   remove() methods, respectively. However, no IdxSet can hold more then
   max() indices, i.e. the number given at the constructor.

   When removing indices, the remaining ones are renumbered. However, all
   indices before the first removed index keep their number unchanged.

   The internal structure of an IdxSet consists of an array #idx storing the
   indices, its length len, and the actually used number of indices #num.
   The class IdxSet doesn't allocate memory for the #idx array. Instead, the
   user has to provide an adequate buffer to the constructor.

   An IdxSet cannot be extended to fit more than max() elements. If
   necessary, the user must explicitely provide the IdxSet with a
   suitable memory. Alternatively, one can use \ref IdxSet "IdxSets"
   which provide the required memory managemant.
*/
class IdxSet
{
protected:

   //---------------------------------------
   /**@name Data */
   //@{
  std::vector<int> idx;           ///< array of indices
   bool freeArray;     ///< true iff \ref soplex::IdxSet::idx "idx" should be freed inside of this object
   //@}

public:

   //---------------------------------------
   /**@name Construction / destruction */
   //@{
   /// constructor.
   /** The constructur receives the index memory \p imem to use for saving
       its indices. This must be large enough to fit \p n indices. \p l can
       be given to construct an #IdxSet initialized to the \p l first
       indices in \p imem.
    */
   IdxSet(int n, int imem[], int l = 0)
      : freeArray(false)
   {
      assert(isConsistent());
   }

   /// default constructor.
   /** The default constructor creates an index set with an empty index
       space. You cannot store any indices in an #IdxSet created with
       the default constructor.
   */
   // IdxSet()
   //    : freeArray(false)
   // {
   //    assert(isConsistent());
   // }

   /// destructor.
   virtual ~IdxSet()
   {
     ;
   }

   /// assignment operator.
   /** The assignment operator copies all nonzeros of the right handside
       #IdxSet to the left one. This implies, that the latter must have
       enough index memory.
    */
   IdxSet& operator=(const IdxSet& set);
   //@}

   //---------------------------------------
   /**@name Access */
   //@{
   /// access \p n 'th index.
   int index(int n) const
   {
     assert(n >= 0 && n < size() && !idx.empty());
      return idx[n];
   }
   /// returns the number of used indices.
   int size() const
   {
     return idx.size();
   }
   /// returns the maximal number of indices which can be stored in IdxSet.
   int max() const
   {
     return idx.capacity();
   }

   /// returns the maximal index.
   int dim() const;

   /// returns the position of index \p i.
   /** Finds the position of the first index \p i in the #IdxSet. If no index \p i is
       available in the #IdxSet, -1 is returned. Otherwise,
       index(pos(\p i)) == \p i holds.
    */
   int pos(int i) const;
   //@}

   //---------------------------------------
   /**@name Modification */
   //@{
   /// appends \p n uninitialized indices.
   void add(int n)
   {
      idx.resize(idx.size() + n);
   }

   /// appends all indices of \p set.
   void add(const IdxSet& set)
   {
     idx.insert(idx.end(), set.idx.begin(), set.idx.end());
   }

   /// appends \p n indices in \p i.
  void add(int n, const std::vector<int> i);

  // Add n elements from a array
  void add(int n, const int *i)
  {
    assert(n >= 0 && size() + n <= max());
    idx.insert(idx.end(), i, i+n);
  }

   /// appends index \p i.
   void addIdx(int i)
   {
      idx.push_back(i);
   }
   /// removes indices at position numbers \p n through \p m.
   void remove(int n, int m);

   /// removes \p n 'th index.
   void remove(int n)
   {
      assert(n >= 0 && n < size());
      // The value of the nth element is set of the last element
      idx[n] = idx.back();
      // The last element is removed. We do this, since the order of elements in
      // idx doesn't matter and this is more efficient (?)
      idx.erase(idx.end()-1);
   }

   /// removes all indices.
   void clear()
   {
     idx.clear();
   }
   //@}

   //---------------------------------------
   /**@name Consistency check */
   //@{
   /// consistency check.
   bool isConsistent() const;
   //@}

  // Functions from old DIdxSet

  void setMax(int newmax = 1)
  {
    if(newmax > idx.capacity())
      {
        idx.reserve(newmax);
      }
  }

  // copy constructor
  explicit IdxSet(const IdxSet& old)
  {
    idx.reserve(old.idx.size());
    // Calling the = operator of std::vector
    idx = old.idx;

    freeArray = true;
  }

  // constructor
  explicit IdxSet(int n = 8)
  {
    idx.reserve(n);
    freeArray = false;
  }

};

} // namespace soplex
#endif // _IDXSET_H_
