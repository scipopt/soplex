/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*  Copyright (c) 1996-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SoPlex; see the file LICENSE. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  islist.h
 * @brief Generic single linked list.
 */
#ifndef _ISLIST_H_
#define _ISLIST_H_

#include <assert.h>
#include <stddef.h>
#include <iostream>


namespace soplex
{

//---------------------------------------------------------------------
//  class IsElement<T>
//---------------------------------------------------------------------

/**@brief   Elements for \ref soplex::IsList "IsList"s.
   @ingroup Elementary

   Class IsElement allows to easily construct list elements for an intrusive
   single linked list IsList out of a template class T. It adds a next
   pointer to each element. An instance of IdElement<T> a can be used just
   like an instance of T itself, except that method next() has been
   added (thereby overriding any method next() defined in T).
 */
template < class T >
class IsElement : public T
{
protected:

   //--------------------------
   /**@name Data */
   ///@{
   IsElement<T>* the_next;       ///< pointer to next element in the IsList.
   ///@}

public:

   //---------------------------------
   /**@name Successor */
   ///@{
   ///
   IsElement<T>*& next()
   {
      return the_next;
   }
   /// returns the next element in the IsList.
   IsElement<T>* next() const
   {
      return the_next;
   }
   ///@}

   //------------------------------------
   /**@name Constructors / destructors */
   ///@{

   /// default constructor.
   IsElement()
   {}

   ///
   explicit
   IsElement(const T& old)
      : T(old)
      , the_next(nullptr)
   {}

   /// copy constructor.
   /** Only the element itself is copied, while the link to the next list
       element is set to a zero pointer.
   */
   IsElement(const IsElement<T>& old)
      : T(old)
      , the_next(nullptr)
   {}
};

//---------------------------------------------------------------------
// class IsList<T>
//---------------------------------------------------------------------

/**@brief   Generic single linked list.
   @ingroup Elementary

   Class IsList implements an intrusive single linked list of elements of a
   template class T. As an \em intrusive list, the objects of type T
   must provide methods next() for setting and inquiring a pointer to the
   next element in a list. The user is responsible for not modifying the
   next() pointer of elements currently residing in a list, which may destroy
   the lists integrity. For this, class IsList provides enough methods for
   modifying a list in a save way. See the method list for a description.
 */
template < class T >
class IsList
{
protected:
   T* the_first;   ///< the first element in the IsList.
   T* the_last;    ///< the last element in the IsList.
   bool destroyElements;
   ///< should the destructor be called for each element when the list is destroyed?

public:
   typedef IsElement<T> Element;

   //--------------------------
   /**@name Extension */
   ///@{
   /// appends \p elem to IsList.
   void append(T* elem)
   {
      if(the_last)
         the_last->next() = elem;
      else
         the_first = elem;

      the_last = elem;
   }

   /// prepends \p elem to IsList.
   void prepend(T* elem)
   {
      if(the_first)
         elem->next() = the_first;
      else
         the_last = elem;

      the_first = elem;
   }

   /// inserts \p elem to IsList after its element \p after.
   void insert(T* elem, T* after)
   {
      assert(find(after));

      if(after == the_last)
         append(elem);
      else
      {
         elem->next() = after->next();
         after->next() = elem;
      }
   }

   /// appends all elements of \p list to IsList.
   /** Appending one list to another keeps the appended \p list. Instead,
       \p list remains an own IsList which is then part of the
       concatenated list. This means that modifying \p list will modify the
       concateneted list as well and vice versa. The programmer is
       responsible for such changes not to yield inconsistent lists.
    */
   void append(IsList<T>& list)
   {
      if(list.the_first)
      {
         append(list.the_first);
         the_last = list.the_last;
      }
   }

   /// prepends all elements of \p list to IsList.
   /** Appending one list to another keeps the appended \p list.  Instead,
       \p list remains an own IsList which is then part of the
       concatenated list. This means that modifying \p list will modify the
       concateneted list as well and vice versa. The programmer is
       responsible for such changes not to yield inconsistent lists.
   */
   void prepend(IsList<T>& list)
   {
      if(list.the_first)
      {
         prepend(list.the_last);
         the_first = list.the_first;
      }
   }

   /// inserts all elements of \p list after element \p after of an IsList.
   /** Inserting one list into another keeps the appended \p list. Instead,
       \p list remains an own IsList which is then part of the
       concatenated list. This means that modifying \p list will modify the
       concateneted list as well and vice versa. The programmer is
       responsible for such changes not to yield inconsistent lists.
   */
   void insert(IsList<T>& list, T* after)
   {
      assert(find(after));

      if(list.the_first)
      {
         list.the_last->next() = after->next();
         after->next() = list.first();

         if(after == last())
            the_last = list.last();
      }
   }
   ///@}

   //--------------------------
   /**@name Removal */
   ///@{
   /// removes the successor of \p after from an IsList.
   void remove_next(T* after)
   {
      assert(find(after));

      if(after->next())
      {
         if(after->next() == last())
            the_last = after;

         after->next() = after->next()->next();
      }
   }

   /// removes element \p elem from an IsList.
   void remove(const T* elem)
   {
      if(the_first)
      {
         if(elem == the_first)
         {
            the_first = next(elem);

            if(the_first == nullptr)
               the_last = nullptr;
         }
         else
         {
            T* after = the_first;

            for(; after != the_last; after = after->next())
               if(after->next() == elem)
               {
                  remove_next(after);
                  return;
               }
         }
      }
   }

   /// removes all elements of \p list from an IsList.
   /** Removing \p list from an IsList requires \p list to be part of the
       IsList. Such a situation can be achieved by previously adding
       (i.e., #append%ing, #insert%ing or #prepend%ing) a list or
       explicitely constructing a sublist with method sublist().
   */
   void remove(IsList<T>& list)
   {
      if(the_first != nullptr && list.the_first != nullptr)
      {
         assert(find(list.first()));
         assert(find(list.last()));

         if(the_first == list.the_first)
         {
            if(the_last == list.the_last)
               the_first = the_last = nullptr;
            else
               the_first = list.the_last->next();
         }
         else
         {
            T* after = the_first;

            for(; after->next() != list.the_first; after = after->next())
               ;

            if(the_last == list.the_last)
               the_last = after;
            else
               after->next() = list.the_last->next();
         }
      }
   }

   /// removes all elements from an IsList.
   void clear(bool pDestroyElements = false)
   {
      if(pDestroyElements)
      {
         T* nextElement;

         for(T* it = the_first; it; it = nextElement)
         {
            nextElement = next(it);
            it->~T();
            spx_free(it);
         }
      }

      the_first = the_last = nullptr;
   }
   ///@}

   //--------------------------
   /**@name Access */
   ///@{
   /// returns the IsList's first element.
   T* first() const
   {
      return the_first;
   }

   /// returns the IsList's last element.
   T* last() const
   {
      return the_last;
   }

   /// returns successor of \p elem in an IsList.
   /** The successor of \p elem in a list generally corresponds to the
       element returned by elem->next(). However, if \p elem is the last
       element in an IsList, this method will return 0, whereas
       elem->next() may yield an arbitrary value. For example, if the
       current list is actually a sublist of another, larger IsList,
       elem->next() returns the successor of \p elem in this larger
       IsList.
    */
   T* next(const T* elem) const
   {
      return (elem == the_last) ? nullptr : elem->next();
   }

   /// returns the number of elements in IsList.
   int length() const
   {
      int num;

      if(the_first)
      {
         T* test = the_first;

         for(num = 1; test != the_last; test = test->next())
            ++num;

         return num;
      }

      return 0;
   }

   /// returns the position of element \p elem within IsList.
   int find(const T* elem) const
   {
      const T* test;

      assert(elem != nullptr);

      for(test = the_first; test != nullptr; test = next(test))
         if(test == elem)
            return 1;

      return 0;
   }

   /// constructs sublist of an IsList.
   /** Returns a new IsList containing a sublist of an IsList starting
       with element \p start and reaching up to element \p end. Both must be
       members of the IsList or 0, in which case the first and last
       element are used, respectively.
    */
   IsList<T>sublist(const T* start = nullptr, const T* end = nullptr) const
   {
      IsList<T>part = *this;

      if(start)
      {
         assert(find(start));
         part.the_first = const_cast<T*>(start);
      }

      if(end)
      {
         assert(part.find(end));
         part.the_last = const_cast<T*>(end);
      }

      return part;
   }
   ///@}

   //--------------------------
   /**@name Miscellaneous */
   ///@{
   /// adjusts list pointers to a new memory address.
   /** This method is of a rather technical nature. If all list elements
       are taken form one array of elements, in certain circumstances the
       user may be forced to realloc this array. As a consequence all
       next() pointers of the list elements would become invalid.
       However, all addresses will be changed by a constant offset \p delta.
       Then move( \p delta ) may be called, which adjusts the next()
       pointers of all elements in the list.
   */
   void move(ptrdiff_t delta)
   {
      if(the_first)
      {
         T* elem;
         the_last  = reinterpret_cast<T*>(reinterpret_cast<char*>(the_last) + delta);
         the_first = reinterpret_cast<T*>(reinterpret_cast<char*>(the_first) + delta);

         for(elem = first(); elem; elem = next(elem))
            if(elem != last())
               elem->next() = reinterpret_cast<T*>(reinterpret_cast<char*>(elem->next()) + delta);
      }
   }

   /// consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS

      if(first() != nullptr && last() == nullptr)
         return SPX_MSG_INCONSISTENT("IsList");

      if(first() == nullptr && last() != nullptr)
         return SPX_MSG_INCONSISTENT("IsList");

      if(first() && find(last()) == nullptr)
         return SPX_MSG_INCONSISTENT("IsList");

#endif

      return true;
   }
   ///@}

   //------------------------------------
   /**@name Constructors / Destructors */
   ///@{
   /// default constructor.
   /** The default constructor may be used to setup a (sub-)list, by
       specifying a \p first and \p last element. Then \p last must be a
       successor of \p first.
   */
   explicit
   IsList(T* pfirst = nullptr, T* plast = nullptr, bool pDestroyElements = false)
      : the_first(pfirst)
      , the_last(plast)
      , destroyElements(pDestroyElements)
   {
      if(pfirst)
      {
         assert(plast != nullptr);
         assert(find(plast));
      }

      assert(isConsistent());
   }

   /// Assignment operator and copy constructor should be deleted to avoid memory problems
   IsList(const IsList<T>&) = delete;
   IsList<T>& operator=(const IsList<T>& old) = delete;

   /// destructor
   ~IsList()
   {
      clear(destroyElements);
   }
   ///@}
};

} // namespace soplex
#endif // _ISLIST_H_
