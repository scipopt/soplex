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
#pragma ident "@(#) $Id: dataset.h,v 1.13 2001/12/04 18:25:56 bzfkocht Exp $"

/**@file  dataset.h
 * @brief Set of data objects.
 */
#ifndef _DATASET_H_
#define _DATASET_H_


#include <stdlib.h>
// #include <memory.h>
#include <assert.h>

#include "dataarray.h"
#include "datakey.h"
#include "spxalloc.h"

namespace soplex
{
/**@brief   Set of data objects.
   @ingroup Elementary

   Class #DataSet manages of sets of \Ref DataObjects "Data Objects" of a
   template type #DATA. For constructing a #DataSet the maximum number 
   of entries must be given. The current maximum number may be inquired 
   with method #max().

   Adding more then #max() elements to a #DataSet will core dump. However,
   method #reMax() allows to reset #max() without loss of elements currently
   in the #DataSet. The current number of elements in a #DataSet is returned
   by method #num().
   
   Adding elements to a #DataSet is done via methods #add() or #create(),
   while #remove() removes elements from a #DataSet. When adding an element
   to a #DataSet the new element is assigned a #Key. #Key%s serve to
   access #DATA elements in a set via a version of the subscript
   #operator[](#Key).
   
   For convenience all elements in a #DataSet are implicitely numbered
   from 0 through #num()-1 and can be accessed with these numbers
   using a 2nd subscript #operator[](int). The reason for providing
   #Key%s to access elements of a #DataSet is that the #Key of an
   element remains unchanged as long as the element is a member of the
   #DataSet, while the numbers will change in an undefined way, if
   other elements are added to or removed from the #DataSet.

   The elements in a #DataSet and their #Key%s are stored in two arrays:
   - #theitem keeps the elements #data along with their number stored in #item.
   - #thekey  keeps the #Key::idx's of the elements in a #DataSet.

   Both arrays have size #themax.

   In #thekey only elements 0 thru #thenum-1 contain #Key::idx%'s of
   valid elements, i.e. elements currently in the #DataSet.
   The current number of elements in the #DataSet is counted in #thenum.
   
   In #theitem only elements 0 thru #thesize-1 are used, but only some of
   them actually contain real data elements of the #DataSet. They are
   recognized by having #info >= 0, which gives the number of that
   element. Otherwise #info < 0 indicates an unused element. Unused
   elements are linked in a single linked list: starting with element
   #-firstfree-1, the next free element is given by #-info-1. The last
   free element in the list is marked by #info == -themax-1. Finally all
   elements in #theitem with index >= #thesize are unused as well.  

   @todo We use here malloc/realloc for objects of type Item/DATA.
         I am not sure that this is a good idea. Maybe we should switch
         do new/delete here.
*/
template<class DATA>
class DataSet
{
private:
   /**@todo replace all "Key"s with "DataKey" and get rid of this typedef */
   typedef DataKey Key;

protected:
   struct Item
   {
      DATA data;       ///< data element
      int info;        ///< element number. info \f$\in\f$ [0,thesize-1] 
                       ///< iff element is used
   }* theitem;         ///< array of elements in the #DataSet

   Key* thekey;        ///< #Key::idx%s of elements

   int themax;         ///< length of arrays #theitem and #thekey
   int thesize;        ///< highest used element in #theitem
   int thenum;         ///< number of elements in #DataSet
   int firstfree;      ///< first unused element in #theitem

public:
   /**@name Extension
    *  Whenever a new element is added to a #DataSet, the latter assigns it a
    *  #Key. For this all methods that extend a #DataSet by one ore more
    *  elements are provided with two signatures, one of them having a
    *  parameter for returning the assigned #Key(s).
    */
   //@{
   /// adds an element.
   /**@return 0 on success and non-zero, if an error occured.
    */
   int add(Key& newkey, const DATA& item)
   {
      DATA* data = create(newkey);
      if (data == 0)
         return 1;
      memcpy(data, &item, sizeof(DATA));
      return 0;
   }
   /// adds element \p item.
   /**@return 0 on success and non-zero, if an error occured.
    */
   int add(const DATA& item)
   {
      DATA* data = create();
      if (data == 0)
         return 1;
      memcpy(data, &item, sizeof(DATA));
      return 0;
   }

   /// add several items.
   /**@return 0 on success and non-zero, if an error occured.
    *
    * @todo The add-methods for multiple items ignore the return value 
    *       of the standard add-methods; they should return a non-zero, 
    *       if one of the add's failed! 
    */
   int add(Key newkey[], const DATA* item, int n)
   {
      assert(n >= 0);
      if (num() + n > max())
         return 1;
      for (int i = 0; i < n; ++i)
         add(newkey[i], item[i]);
      return 0;
   }

   /// adds \p n elements from \p items.
   /**@return 0 on success and non-zero, if an error occured.
    */
   int add(const DATA* items, int n)
   {
      assert(n >= 0);
      if (num() + n > max())
         return 1;
      for (int i = 0; i < n; ++i)
         add(items[i]);
      return 0;
   }

   /// adds several new items.
   int add(Key newkey[], const DataSet < DATA > & set)
   {
      if (num() + set.num() > max())
         return 1;
      for (int i = 0; i < set.num(); ++i)
         add(newkey[i], set[i]);
      return 0;
   }

   /// adds all elements of \p set.
   /**@return 0 on success and non-zero, if an error occured.
    */
   int add(const DataSet < DATA > & set)
   {
      if (num() + set.num() > max())
         return 1;
      for (int i = 0; i < set.num(); ++i)
         add(set[i]);
      return 0;
   }

   /// creates new data element in #DataSet.
   /**@return Pointer to the newly created element.
    */
   DATA* create(Key& newkey)
   {
      if (num() >= max())
         return 0;

      if (firstfree != -themax - 1)
      {
         newkey.idx = -firstfree - 1;
         firstfree = theitem[newkey.idx].info;
      }
      else
         newkey.idx = thesize++;

      thekey[thenum] = newkey;
      theitem[newkey.idx].info = thenum;
      ++thenum;

      return &(theitem[newkey.idx].data);
   }
   /// creates new (uninitialized) data element in #DataSet.
   /**@return Pointer to the newly created element.
    */
   DATA* create()
   {
      Key tmp;
      return create(tmp);
   }
   //@}

   /**@name Shrinkage
    * When elements are removed from a #DataSet, the remaining ones are
    * renumbered from 0 through the new #size()-1. How this renumbering is
    * performed will not be revealed, since it might be target of future
    * changes. However, some methods provide a parameter int* #perm, which
    * returns the new order after the removal: If #perm[i] < 0, the element
    * numbered i prior to the removal operation has been removed from the
    * set. Otherwise, #perm[i] = j >= 0 means, that the element with number
    * i prior to the removal operation has been renumberd to j.
    * Removing a single elements from a #DataSet yields a simple
    * renumbering of the elements: The last element in the set (i.e.
    * element #num()-1) is moved to the index of the removed element.
    */
   //@{
   /// removes the \p removenum 'th element.
   void remove(int removenum)
   {
      if (has(removenum))
      {
         int idx = thekey[removenum].idx;

         theitem[idx].info = firstfree;
         firstfree = -idx - 1;

         while (-firstfree == thesize)
         {
            firstfree = theitem[ -firstfree - 1].info;
            --thesize;
         }

         --thenum;
         if (removenum != thenum)
         {
            thekey[removenum] = thekey[thenum];
            theitem[thekey[removenum].idx].info = removenum;
         }
      }
   }

   /// removes element with key \p removekey.
   void remove(const Key& removekey)
   {
      remove(number(removekey));
   }

   /// removes element \p item.
   void remove(DATA& item)
   {
      remove(number(&item));
   }

   /// remove multiple elements.
   /** This method removes all elements for the #DataSet with an
    *  index i such that #perm[i] < 0. Upon completion, #perm contains
    *  the new numbering of elements.
    */
   void remove(int perm[])
   {
      int k, j, first = -1;
      // setup permutation and remove items
      for (k = j = 0; k < num(); ++k)
      {
         if (perm[k] >= 0)      // #j# has not been removed ...
            perm[k] = j++;
         else
         {
            int idx = thekey[k].idx;
            theitem[idx].info = firstfree;
            firstfree = -idx - 1;
            if (first < 0)
               first = k;
         }
      }
      if (first >= 0)        // move remaining items
      {
         for (k = first, j = num(); k < j; ++k)
         {
            if (perm[k] >= 0)
            {
               thekey[perm[k]] = thekey[k];
               theitem[thekey[k].idx].info = perm[k];
               thekey[k].idx = -1;
            }
            else
               --thenum;
         }
      }
   }

   /// ???
   void remove(Key *keys, int n, int* perm)
   {
      assert(perm != 0);
      for (int i = num() - 1; i >= 0; --i)
         perm[i] = i;
      while (--n >= 0)
         perm[number(keys[n])] = -1;
      remove(perm);
   }
   /// remove \p n elements given by \p keys.
   void remove(Key *keys, int n)
   {
      DataArray<int> perm(num());
      remove(keys, n, perm.get_ptr());
   }

   /// ???
   void remove(int *nums, int n, int* perm)
   {
      assert(perm != 0);
      for (int i = num() - 1; i >= 0; --i)
         perm[i] = i;
      while (--n >= 0)
         perm[nums[n]] = -1;
      remove(perm);
   }
   /// remove \p n elements with numbers \p nums.
   void remove(int *nums, int n)
   {
      DataArray<int> perm(num());
      remove(nums, n, perm.get_ptr());
   }

   /// remove all elements.
   void clear()
   {
      thesize = 0;
      thenum = 0;
      firstfree = -themax - 1;
   }
   //@}

   /**@name Access   
    * When accessing elements from a #DataSet with one of the index
    * operators, it must be ensured, that the index is valid for the
    * #DataSet. If this is not known afore, it is the programmers
    * responsability to ensure this using the inquiry methods below.
    */
   //@{
   ///
   DATA& operator[](int n)
   {
      return theitem[thekey[n].idx].data;
   }
   /// returns element number \p n.
   const DATA& operator[](int n) const
   {
      return theitem[thekey[n].idx].data;
   }

   ///
   DATA& operator[](const Key& k)
   {
      return theitem[k.idx].data;
   }
   /// returns element with #Key \p k.
   const DATA& operator[](const Key& k) const
   {
      return theitem[k.idx].data;
   }
   //@}

   /**@name Inquiry */
   //@{
   /// returns maximum number of elements that would fit into #DataSet.
   int max() const
   {
      return themax;
   }

   /// returns number of elements currently in #DataSet.
   int num() const
   {
      return thenum;
   }

   /// returns the maximum #Key::idx currently in #DataSet.
   int size() const
   {
      return thesize;
   }

   /// returns #Key of \p n 'th element in #DataSet.
   Key key(int n) const
   {
      assert(n >= 0 && n < num());
      return thekey[n];
   }

   /// returns #Key of element \p item in #DataSet.
   Key key(const DATA* item) const
   {
      assert(number(item) >= 0);
      return thekey[number(item)];
   }

   /// returns the number of the element with #Key \p k in #DataSet or -1, 
   /// if it doesn't exist.
   int number(const Key& k) const
   {
      return (k.idx < 0 || k.idx >= size()) ? -1
          : theitem[k.idx].info;
   }

   /**@todo Please check, whether this is correctly implemented! */
   /// returns the number of element \p item in #DataSet or -1, 
   /// if it doesn't exist.
   int number(const DATA* item) const
   {      
      ptrdiff_t idx = reinterpret_cast<const struct Item*>(item) - theitem;

      if( idx < 0 || idx >= size())
         return -1;

      /* ??? old code:
         if ((reinterpret_cast<unsigned long>(item) < 
              reinterpret_cast<unsigned long>(theitem) )
             || 
             ( reinterpret_cast<unsigned long>(item) >= 
               reinterpret_cast<unsigned long>(&(theitem[size()]))))
            return -1;
         long idx = ((reinterpret_cast<long>(item)) 
                  - (reinterpret_cast<long>(theitem))) 
                    / sizeof(Item);
      */
      return theitem[idx].info;
   }

   /// Is \p k a valid #Key of an element in #DataSet?
   int has(const Key& k) const
   {
      return theitem[k.idx].info >= 0;
   }

   /// Is \p n a valid number of an element in #DataSet?
   int has(int n) const
   {
      return (n >= 0 && n < num());
   }

   /// Does \p item belong to #DataSet?
   int has(const DATA* item) const
   {
      return number(item) >= 0;
   }
   //@}

   /**@name Miscellaneous */
   //@{
   /// resets #max() to \p newmax.
   /** This method will not succeed if \p newmax < #size(), in which case
    *  \p newmax == #size() will be taken. As generally this method involves
    *  copying the #DataSet%s elements in memory, #reMax() returns the
    *  number of bytes the addresses of elements in the #DataSet have been
    *  moved. Note, that this is identical for all elements in the
    *  #DataSet.
    */
   ptrdiff_t reMax(int newmax = 0)
   {
      struct Item * old_theitem = theitem;
      newmax = (newmax < size()) ? size() : newmax;

      int* lastfree = &firstfree;
      while (*lastfree != -themax - 1)
         lastfree = &(theitem[ -1 - *lastfree].info);
      *lastfree = -newmax - 1;

      themax = newmax;

      spx_realloc(theitem, themax);
      spx_realloc(thekey,  themax);

      return reinterpret_cast<char*>(theitem) 
         - reinterpret_cast<char*>(old_theitem);
   }

   /// consistencty check.
   int isConsistent() const
   {
      if (theitem == 0 || thekey == 0)
      {
         std::cerr << "Inconsistency detected at " 
                   << __FILE__ << "(" << __LINE__ << ")" << std::endl;
         return 0;
      }

      if (thesize > themax || thenum > themax || thenum > thesize)
      {
         std::cerr << "Inconsistency detected at " 
                   << __FILE__ << "(" << __LINE__ << ")" << std::endl;
         return 0;
      }

      if (thesize == thenum && firstfree != -themax - 1)
      {
         std::cerr << "Inconsistency detected at " 
                   << __FILE__ << "(" << __LINE__ << ")" << std::endl;
         return 0;
      }

      if (thesize != thenum && firstfree == -themax - 1)
      {
         std::cerr << "Inconsistency detected at " 
                   << __FILE__ << "(" << __LINE__ << ")" << std::endl;
         return 0;
      }

      for (int i = 0; i < thenum; ++i)
      {
         if (theitem[thekey[i].idx].info != i)
         {
            std::cerr << "Inconsistency detected at " 
                      << __FILE__ << "(" << __LINE__ << ")" << std::endl;
            return 0;
         }
      }
      return 1;
   }
   //@}


   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   DataSet(int pmax = 8)
   {
      themax = (pmax < 1) ? 8 : pmax;
      thesize = thenum = 0;
      firstfree = -themax - 1;

      spx_alloc(theitem, themax);
      spx_alloc(thekey, themax);
   }

   /// copy constructor.
   DataSet(const DataSet& old)
   {
      themax = old.themax;
      thesize = old.thesize;
      thenum = old.thenum;
      if (old.firstfree == -old.themax - 1)
         firstfree = -themax - 1;
      else
         firstfree = old.firstfree;

      spx_alloc(theitem, themax);
      spx_alloc(thekey, themax);

      memcpy(theitem, old.theitem, themax * sizeof(Item));
      memcpy(thekey, old.thekey, themax * sizeof(Key));
   }

   /// assignment operator.
   /** The assignment operator involves #reMax()%ing the lvalue #DataSet
    *  to the size needed for copying all elements of the rvalue. After the
    *  assignment all #Key%s from the lvalue are valid for the rvalue as
    *  well. They refer to a copy of the corresponding data elements.
    */
   DataSet < DATA > & operator=(const DataSet < DATA > & rhs)
   {
      int i;
      if (rhs.size() > max())
         reMax(rhs.size());
      clear();
      for (i = 0; i < rhs.size(); ++i)
         memcpy(&(theitem[i]), &(rhs.theitem[i]), sizeof(Item));
      for (i = 0; i < rhs.num(); ++i)
         thekey[i] = rhs.thekey[i];

      if (rhs.firstfree == -rhs.themax - 1)
         firstfree = -themax - 1;
      else
      {
         firstfree = rhs.firstfree;
         i = rhs.firstfree;

         while (rhs.theitem[ -i - 1].info != -rhs.themax - 1)
            i = rhs.theitem[ -i - 1].info;
         theitem[ -i - 1].info = -themax - 1;
      }
      thenum = rhs.thenum;
      thesize = rhs.thesize;

      return *this;
   }

   /// destructor.
   ~DataSet()
   {
      spx_free(theitem);
      spx_free(thekey);
   }
   //@}
};

} // namespace soplex
#endif // _DATASET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
