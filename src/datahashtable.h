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
#pragma ident "@(#) $Id: datahashtable.h,v 1.10 2002/01/04 17:31:38 bzfkocht Exp $"

/**@file  datahashtable.h
 * @brief Generic hash table for data objects.
 */
#ifndef _DATAHAHSTABLE_H_
#define _DATAHAHSTABLE_H_

#include <iostream>
#include <assert.h>

namespace soplex
{
/**@brief   Generic hash table for data objects.
   @ingroup Elementary

   Class DataHashTable provides a generic hash table for 
   \ref DataObjects "Data Objects",
   i.e. a map that maps arguments called #HashItem%s to values called #Info%s.
   #HashItem and #Info types are passed as template arguments. #HashItem%s
   must provide a comparision #operator==().  Further both, the HashItem and
   Info must be data objects in the sense, that the assignment operator is
   equivalent to a #memcpy() of the structure and no destructor is required.
   
   The construction of a #DataHashTable requires a \em hash \em function that
   assigns every #HashItem to an integer value.  Provided this, pairs of a
   #HashItem and a #Info can be added to the #DataHashTable. No more
   than one #Info to the same #HashItem is possible at a time. The #Info
   to a #HashItem can be accessed through the subscript #operator[]() with
   #Info as subscript.
   
   A #DataHashTable can hold up to #max() entries. This #max() can be
   specified upon construction and may be reset with #reMax() later on.
   Further, a value #hashSize is required. This value must be #< max() and must
   not have a common dominator with #max(). If not specified explicitely, it
   is set automatically to a reasonable value. It may be reset with method
   #reHash(). Note that both, #reMax()%ing and #reHash()%ing renders any
   reference to entries of a #DataHashTable invalid. This also happens, if the
   #DataHashTable is #reMax()%ed automatically, if more than #max() entries
   are added to it.

   The implementation relies on an array of #DataHashTable::Element%s, from
   now on referred to as elements. Upon construction, all elements are
   marked #FREE in their member #status. When an entry is added
   to the #DataHashTable, the hash value is computed by calling #hashval
   for its #HashItem. If this array element is unused, it is
   taken right away. Otherwise, the array index is incremented by
   #hashsize (modulo the element array #size()%) until an unused element
   is found.
   
   Removing elements is simply done by marking it as #RELEASED. Hence,
   when searching for an element, the search loop may not stop, when a
   #RELEASED element is encountered. However, such an element may be
   reused when adding a new element to the #DataHashTable.
   
   Further, memory management with resizing of the element array is
   straight forward.
*/
template < class HashItem, class Info >
class DataHashTable
{
private:

   /// template class for elements stored in the hash table
   template < class ElemHashItem, class ElemInfo >
   class Element
   {
   public:
      ElemHashItem item;
      ElemInfo     info;
      enum
      {
         FREE,            ///< element has never been used
         RELEASED,        ///< element had been used, but released
         USED             ///< element is in use
      } status;
   };
   /// stores all elements of the hash table
   DataArray < Element < HashItem, Info > > element;   

   int hashsize;           ///< increment added to hash index, if allready used
   int thenum;             ///< current number of entries in the hash table

   /// pointer to hash function (mapping: #HashItem -> int)
   int (*hashval) (const HashItem*);  

   /// memory is #reMax()%ed by this factor, if a new element does't fit
   double factor;  

   mutable int theCurrent; ///< index for iterator
   
public:
   /**@name Inquiry Methods */
   //@{
   /// number of elements that would fit
   int max () const
   {
      return element.size();
   }

   /// number of hashed elements
   int num () const
   {
      return thenum;
   }

   /// returns #hashsize, i.e. the increment when searching for elements
   int hashSize () const
   {
      return hashsize;
   }

   /// Is item \p h present in #DataHashTable ?
   int has (const HashItem& h) const
   {
      return index(h) >= 0 ? 1 : 0;
   }
   //@}

   /**@name Access Methods */
   //@{
   /// returns pointer to #Info of #HashItem \p h or 0, if \p h is not found.
   /** Returns a pointer to #Info component of hash element \p h or a zero
       pointer if element \p h is not in the table.
    */
   Info* get (const HashItem& h)
   {
      int i = index(h);
      return i >= 0 ? &element[i].info : 0;
   }
   /// returns const pointer to #Info of #HashItem \p h or 0, 
   /// if item is not found.
   /** Returns a pointer to #Info component of hash element \p h or a zero
       pointer if element \p h is not in the table.
    */
   const Info* get (const HashItem& h) const
   {
      int i = index(h);
      return i >= 0 ? &element[i].info : 0;
   }

   /// references #Info of #HashItem \p h.
   /** Index operator for accessing the #Info associated to
       #HashItem \p h. It is required, that \p h belongs to the
       #DataHashTable, otherwise it core dumps. Methods #has() or
       #get() can be used for inquiring wheater \p h belongs to the
       #DataHashTable or not.
   */
   Info& operator[](const HashItem& h)
   {
      return element[index(h)].info;
   }
   /// references #Info of #HashItem \p h.
   /** Index operator for accessing the #Info associated to
       #HashItem \p h. It is required, that \p h belongs to the
       #DataHashTable, otherwise it core dumps. Methods #has() or
       #get() can be used for inquiring wheater \p h belongs to the
       #DataHashTable or not.
   */
   const Info& operator[](const HashItem& h) const
   {
      return element[index(h)].info;
   }
   //@}

   /**@name Iteration
    *  Often it is desired to loop through all elements in a #DataHashTable.
    *  This is provided by means of the following 5 methods. They imply an
    *  arbitray order to all elements currently in the #DataHashTable. This
    *  order may change after any non const member function invocation. When
    *  calling one of these methods, a marker is set that serves as reference
    *  point for the next call.
    *  All iteration methods return 0, if the marker points to a non existing
    *  element (for example, calling #next() returns 0, if there doesn't exist
    *  another element).
    */
   //@{
   /// returns first #HashItem in hash table and sets marker to it.
   const HashItem* first() const
   {
      theCurrent = -1;
      return next();
   }
   /// returns last #Item in hash table and sets marker to it.
   const HashItem* last() const
   {
      theCurrent = element.size();
      return prev();
   }
   /// returns #HashItem following current marker and increasing the marker.
   const HashItem* next() const
   {
      if (theCurrent < 0)
         theCurrent = -1;
      while (++theCurrent < element.size())
      {
         if (element[theCurrent].status
              == Element < HashItem, Info > ::USED)
            return &element[theCurrent].item;
      }
      theCurrent = -1;
      return 0;
   }
   /// returns #HashItem referenced by current marker.
   const HashItem* current() const
   {
      return (theCurrent < 0) ? 0 : &element[theCurrent].item;
   }
   /// returns #HashItem preceding current marker thereby decreasing marker.
   const HashItem* prev() const
   {
      if (theCurrent > element.size())
         theCurrent = element.size();
      while (--theCurrent >= 0)
      {
         if (element[theCurrent].status
              == Element < HashItem, Info > ::USED)
            return &element[theCurrent].item;
      }
      return 0;
   }
   //@}

   /**@name Manipulation Methods */
   //@{
   /// adds a new entry to the hash table.
   /** Adds a new entry consisting of #HashItem \p h and #Info \p x to the
    *  #DataHashTable. No entry with #HashItem \p h must yet be in the
    *  #DataHashTable. After completion, \p x may be accessed via #get() or
    *  #operator[]() with \p h as parameter. The #DataHashTable is #reMax()%ed
    *  if it becomes neccessary.
    */
   void add (const HashItem& h, const Info& x)
   {
      assert(!has(h));
      int i;

      if (thenum >= element.size())
         reMax(int(factor * thenum) + 1);

      assert(element.size() > 0);

      for(
         i = (*hashval)(&h) % element.size();
         element[i].status == Element < HashItem, Info > ::USED;
         i = (i + hashsize) % element.size())
        ;
      element[i].status = Element < HashItem, Info > ::USED;
      memcpy(&(element[i].info), &x, sizeof(Info));
      memcpy(&(element[i].item), &h, sizeof(HashItem));
      ++thenum;
   }

   /// remove #HashItem \p h from the #DataHashTable.
   void remove (const HashItem& h)
   {
      assert(has(h));
      element[index(h)].status
      = Element < HashItem, Info > ::RELEASED;
      --thenum;
   }

   /// remove all entries from #DataHashTable.
   void clear ()
   {
      for (int i = element.size() - 1; i >= 0; --i)
         element[i].status
         = Element < HashItem, Info > ::FREE;
      thenum = 0;
   }

   /// reset #max() and #hashSize().
   /** Reset the #max() of a #DataHashTable to \p nel. However, if
    *  \p nel < #num(), it is resized to #num() only. If \p hashsze < 1, a
    *  new hash size is computed automatically. Otherwise, the specified
    *  value will be taken.
    */
   void reMax (int nel = -1, int hashsze = 0)
   {
      DataArray < Element < HashItem, Info > > cpy(element);
      element.reSize(nel < num() ? num() : nel);
      clear();
      if (hashsze < 1)
         this->hashsize = autoHashSize();
      else
         this->hashsize = hashsze;
      for (int i = cpy.size() - 1; i >= 0; --i)
         if (cpy[i].status == Element < HashItem, Info > ::USED)
            add(cpy[i].item, cpy[i].info);
   }
   //@}

   /**@name Miscellaneous */
   //@{
   /// checks, whether #DataHashTable is consistent
   int isConsistent () const
   {
      int i, tot;

      for (i = element.size() - 1, tot = 0; i >= 0; --i)
         if (element[i].status
              == Element < HashItem, Info > ::USED)
         {
            ++tot;
            if (!has(element[i].item))
            {
               std::cout << "Inconsistency detected in class DataHashTable\n";
               return 0;
            }
         }

      if (tot != thenum)
      {
         std::cout << "Inconsistency detected in class DataHashTable\n";
         return 0;
      }
      return element.isConsistent();
   }

#ifdef DEFINE_OUTPUT_OPERATOR
   /// Output operator. Displays all elements contained in hash table.
   /**@todo Is there any reason not to define this operator? */
   friend std::ostream& operator<<(std::ostream& out,
      const DataHashTable < HashItem, Info > & h)
   {
      const HashItem* item;
      for (item = h.first(); item; item = h.next())
         out << "    " << *item << "\t\t" << h[*item] << std::endl;
      return out;
   }
#endif
   //@}

private:
   /// automatically computes a good #hashsize.
   /** Computes a good #hashsize as the product of all prime numbers 
    *  not divisors of #size() that are <= the maximum divisor of #size().
    *  @return good value for #hashsize
    */
   int autoHashSize() const
   {
      int i, j;
      int hashsze = 1;
      int size = element.size();
      DataArray < char > prime(size);

      for (i = 2; i < size; ++i)
         prime[i] = 1;

      for (i = 2; i < size; ++i)
      {
         if (prime[i])
         {
            for (j = i; j < size; j += i)
               prime[j] = 0;
            if (size % i != 0)
            {
               hashsze *= i;
               if (hashsze > size)
               {
                  hashsze /= i;
                  break;
               }
            }
         }
      }

      return hashsze;
   }

   /// returns hash index of #HashItem \p h or -1, if \p h is not present.
   /** Using the hash function #hashval, the hash value of \p h is calculated.
    *  Starting with this hash index, every #hashsize%-th #element is
    *  compared with \p h until \p h is found or all #element%s are checked.
    *
    *  @param  h  #HashItem, for which the hash index should be calculated
    *  @return hash index of \p h or -1, 
    *          if \p h is not a member of the hash table
    */
   int index(const HashItem& h) const
   {
      int i, j;

      if (thenum == 0)
         return -1;

      assert(element.size() > 0);

      for(i = j = (*hashval)(&h) % element.size();
          element[i].status != Element < HashItem, Info > ::FREE;)
      {
         if (element[i].item == h)
            return i;
         i = (i + hashsize) % element.size();
         if (i == j)
            break;
      }
      return -1;
   }

public:
   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   /** Allocates a #DataHashTable for \p nel entries using \p f as hash
    *  function. If \p hashsze > 0, #hashSize() is set to the specified
    *  value, otherwise a suitable #hashSize() is computed automatically.
    *  Parameter \p incr is used for memory management: If more than
    *  \p nel entries are added to the #DataHashTable, it will
    *  automatically be #reMax()%ed by a factor of \p incr.
    *
    *  @param f            pointer to hash function.
    *  @param nel          number of hash elements.
    *  @param hashsze      hash size.
    *  @param incr         factor for increasing data block.
    */
   DataHashTable
   (int (*f)(const HashItem*),
     int nel = 256 ,
     int hashsze = 0 ,
     double incr = 2.0
  )
      : element(nel)
         , hashval(f)
         , factor (incr)
   {
      clear();
      if (hashsze < 1)
         this->hashsize = autoHashSize();
      else
         this->hashsize = hashsze;
      assert(factor > 1);
   }
   //@}
};
} // namespace soplex
#endif   // _DATAHAHSTABLE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------




