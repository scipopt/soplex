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
#pragma ident "@(#) $Id: dataset.h,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

#ifndef _DATASET_H_
#define _DATASET_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <stdlib.h>
// #include <memory.h>
#include <assert.h>


/*  and class header files
 */
#include "dataarray.h"

namespace soplex
{



//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
    The following declaration of #Key# is made this complicated in order to
    get this code thru AT\&T's cfront compiler. However, #Key# should better
    be declared as a local class of #DataSet# only.
 */
struct DataSet_Key
{
signed int info:
   8;
signed int idx:
   (8*sizeof(int) - 8);

   int isValid() const
   {
      return idx >= 0;
   }

   void inValidate()
   {
      idx = -1;
      info = 0;
   }

   DataSet_Key()
      : info(0)
         , idx(-1)
   {}

};

/** Set of data objects.
    Class #DataSet# manages of sets of \Ref{Data Objects} of a template type
    #DATA#. For constructing a #DataSet# the maximum number of entries must
    be given. The current maximum number may be inquired with method #max()#.
    Adding more then #max()# elements to a #DataSet# will core dump. However,
    method #reMax()# allows to reset #max()# without loss of elements currently
    in the #DataSet#. The current number of elements in a #DataSet# is returned
    by method #num()#.
 
    Adding elements to a #DataSet# is done via methods #add()# or #create()#,
    while #remove()# removes elements from a #DataSet#. When adding an element
    to a #DataSet# the new element is assigned a \Ref{Key}. #Key#s serve to
    access #DATA# elements in a set via a version of the subscript
    #operator[](Key)#.
 
    For convenience all elements in a #DataSet# are implicitely numbered from 0
    through #num()-1# and can be accessed with these numbers using a 2nd subscript
    #operator[](int)#. The reason for providing #Key#s to access elements of a
    #DataSet# is that the #Key# of an element remains unchanged as long as the
    element is a member of the #DataSet#, while the numbers will change in an
    undefined way, if other elements are added to or removed from the
    #DataSet#.
 */
template<class DATA>
class DataSet
{
public:
#ifdef  FOR_DOCXX
   /** Entry Identifiers.
       Every item in a #DataSet# is assigned a #Key# by which it can be
       accessed (using #DataSet::operator[]#). A #Key# consists of an integer
       member #idx#, which is a positive number for any valid #Key#. No
       #Key::idx# of an element in a #DataSet# may exceed the sets #max()#.
       This property may be used to build arrays with additional information to
       the elements of a #DataSet#.

       In addition, #Key#s provides member #info#, where the programmer is
       free to store other information.

       Each #Key# is unique for one #DataSet# but different #DataSet#s may (and
       generally will) manage the same #Key#s. When an element is removed from
       a #DataSet# its #Key# may (and generally will) be reused for other
       elements added to the #DataSet# later on.
    */
   struct Key
   {
      ///
signed int info:
      8;
      ///
signed int idx:
      (8*sizeof(int) - 8);
      ///
      int isValid() const;
      ///
      void inValidate();
   };
#endif  // FOR_DOCXX
   typedef DataSet_Key Key;

protected:
   /*
       \SubSection{Data Structures and Layout}
       The elements in a #DataSet# and their #Key#s are stored in two arrays:
       \begin{itemize}
       \item   \strut #theitem# keeps the elements #data# along with
               their number stored in #item#.
       \item   \strut #thekey# keeps the #Key::idx#'s of the
               elements in a #DataSet#.
       \end{itemize}
       Both arrays have size #themax#.

       In #thekey# only elements 0 thru #thenum-1# contain #Key::idx#'s of
       valid elements, i.e. elements currently in the #DataSet# --- #thenum# is
       the current number of elements in the #DataSet#.

       In #theitem# only elements 0 thru #thesize-1# are used, but only some of
       them actually contain real data elements of the #DataSet#. They are
       recognized by having #info >= 0#, which gives the number of that
       element. Otherwise #info < 0# indicates an unused element. Unused
       elements are linked in a single linked list: starting with element
       #-firstfree-1#, the next free element is given by #-info-1#. The last
       free element in the list is marked by #info == -themax-1#. Finally all
       elements in #theitem# with index #>= thesize# are unused as well.
    */
   struct Item
   {
      DATA data;           // data element
      int info;           //
   }* theitem;       // array of elements in the #DataSet#

   Key* thekey;         // #Key::idx#s of elements

   int themax;         // length of arrays #theitem# and #thekey#
   int thesize;        // highest used element in #theitem#
   int thenum;         // number of elements in #DataSet#
   int firstfree;      // first unused element in #theitem#

public:
   /**@name Extension
       Whenever a new element is added to a #DataSet#, the latter assigns it a
       \Ref{Key}. For this all methods that extend a #DataSet# by one ore more
       elements are provided with two signatures, one of them having a
       parameter for returning the assigned #Key#(s).
    */
   //@{
   ///
   int add(Key& newkey, const DATA& item)
   {
      DATA* data = create(newkey);
      if (data == 0)
         return 1;
      memcpy(data, &item, sizeof(DATA));
      return 0;
   }

   /// add element #item#.
   int add(const DATA& item)
   {
      DATA* data = create();
      if (data == 0)
         return 1;
      memcpy(data, &item, sizeof(DATA));
      return 0;
   }

   ///
   int add(Key newkey[], const DATA* item, int n)
   {
      assert(n >= 0);
      if (num() + n > max())
         return 1;
      for (int i = 0; i < n; ++i)
         add(newkey[i], item[i]);
      return 0;
   }

   /// add #n# elements from #items#.
   int add(const DATA* items, int n)
   {
      assert(n >= 0);
      if (num() + n > max())
         return 1;
      for (int i = 0; i < n; ++i)
         add(items[i]);
      return 0;
   }

   ///
   int add(Key newkey[], const DataSet < DATA > & set)
   {
      if (num() + set.num() > max())
         return 1;
      for (int i = 0; i < set.num(); ++i)
         add(newkey[i], set[i]);
      return 0;
   }

   /// add all elements of #set#.
   int add(const DataSet < DATA > & set)
   {
      if (num() + set.num() > max())
         return 1;
      for (int i = 0; i < set.num(); ++i)
         add(set[i]);
      return 0;
   }

   ///
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

   /// create new (uninitialized) data element in #DataSet#.
   DATA* create()
   {
      Key tmp;
      return create(tmp);
   }
   //@}

   /**@name Shrinkage
       When elements are removed form a #DataSet#, the remaining ones are
       renumbered from 0 through the new #size()#. How this renumbering is
       performed will not be revield, since it might be target of future
       changes. However, some methods provide a parameter #int* perm#, which
       returns the new order after the removal: If #perm[i] < 0#, the element
       numbered $i$ prior to the removal operation has been removed from the
       set. Otherwise, #perm[i] = j >= 0# means, that the element with number
       #i# prior to the removal operation has been renumberd to #j#.
    */
   //@{
   /**@name \  */
   /// remove #removenum#'th element.
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

   /// Remove element to #removekey#.
   void remove(const Key& removekey)
   {
      remove(number(removekey));
   }

   /** Remove element #item#.
       Removing a single elements from a #DataSet# yields a simple
       renumbering of the elements: The last element in the set (i.e.
       element #num()-1# is moved to the index of the removed element.
    */
   void remove(DATA& item)
   {
      remove(number(&item));
   }

   /** remove multiple elements.
       This method removes all elements for the #DataSet# with an
       index #i# such that #perm[i] < 0#. Upon completion, #perm# contains
       the new numbering of elements.
    */
   void remove(int perm[])
   {
      int i, j, first = -1;
      // setup permutation and remove items
      for (i = j = 0; i < num(); ++i)
      {
         if (perm[i] >= 0)      // #j# has not been removed ...
            perm[i] = j++;
         else
         {
            int idx = thekey[i].idx;
            theitem[idx].info = firstfree;
            firstfree = -idx - 1;
            if (first < 0)
               first = i;
         }
      }
      if (first >= 0)        // move remaining items
      {
         for (i = first, j = num(); i < j; ++i)
         {
            if (perm[i] >= 0)
            {
               thekey[perm[i]] = thekey[i];
               theitem[thekey[i].idx].info = perm[i];
               thekey[i].idx = -1;
            }
            else
               --thenum;
         }
      }
   }

   ///
   void remove(Key *keys, int n, int* perm)
   {
      assert(perm);
      for (int i = num() - 1; i >= 0; --i)
         perm[i] = i;
      while (--n >= 0)
         perm[number(keys[n])] = -1;
      remove(perm);
   }

   /// remove #n# elements in #keys#.
   void remove(Key *keys, int n)
   {
      DataArray<int>perm(num());
      remove(keys, n, perm.get_ptr());
   }

   ///
   void remove(int *nums, int n, int* perm)
   {
      assert(perm);
      for (int i = num() - 1; i >= 0; --i)
         perm[i] = i;
      while (--n >= 0)
         perm[nums[n]] = -1;
      remove(perm);
   }

   /// remove #n# elements with numbers #nums#.
   void remove(int *nums, int n)
   {
      DataArray<int>perm(num());
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
       When accessing elements from a #DataSet# with one of the index
       operators, it must be ensured, that the index is valid for the
       #DataSet#. If this is not known afore, it is the programmers
       responsability to ensure this using the inquiry methods below.
    */
   //@{
   ///
   DATA& operator[](int n)
   {
      return theitem[thekey[n].idx].data;
   }

   /// return element nubmer #n#.
   const DATA& operator[](int n) const
   {
      return theitem[thekey[n].idx].data;
   }

   ///
   DATA& operator[](const Key& k)
   {
      return theitem[k.idx].data;
   }

   /// return element to key #k#.
   const DATA& operator[](const Key& k) const
   {
      return theitem[k.idx].data;
   }
   //@}


   /**@name Inquiry */
   //@{
   /// maximum number of elements that would fit into #DataSet#.
   int max() const
   {
      return themax;
   }

   /// number of elements currently in #DataSet#.
   int num() const
   {
      return thenum;
   }

   /// maximum #Key::idx# currently in #DataSet#.
   int size() const
   {
      return thesize;
   }

   /// return #Key# of #n#'th element in #DataSet#.
   Key key(int n) const
   {
      assert(n >= 0 && n < num());
      return thekey[n];
   }

   /// return #Key# to element #item# in #DataSet#.
   Key key(DATA* item) const
   {
      assert(number(item) >= 0);
      return thekey[number(item)];
   }

   /// return number of element with #Key k# in #DataSet#.
   int number(const Key& k) const
   {
      return (k.idx < 0 || k.idx >= size()) ? -1
          : theitem[k.idx].info;
   }

   /// return number of element #item# in #DataSet#.
   int number(DATA* item) const
   {
      if ((unsigned long)item < (unsigned long)theitem
          || (unsigned long)item >= (unsigned long)&(theitem[size()]))
         return -1;
      long idx = (((long)item) - ((long)theitem)) / sizeof(Item);
      return theitem[idx].info;
   }

   /// is #k# valid #Key# of an element in #DataSet#?.
   int has(const Key& k) const
   {
      return theitem[k.idx].info >= 0;
   }

   /// is #n# valid number of an element in #DataSet#?.
   int has(int n) const
   {
      return (n >= 0 && n < num());
   }

   /// does #item# belong to #DataSet#?.
   int has(DATA* item) const
   {
      return number(item) >= 0;
   }
   //@}

   /**@name Miscellaneous */
   //@{
   /** Reset #max()# to #newmax#.
       This method will not succeed if #newmax < size()#, in which case
       #newmax == size# will be taken. As generally this method involves
       copying the #DataSet#s elements in memory, #reMax()# returns the
       number of bytes the addresses of elements in the #DataSet# have been
       moved. Note, that this is identical for all elements in the
       #DataSet#.
    */
   long reMax(int newmax = 0)
   {
      long delta = long(theitem);
      newmax = (newmax < size()) ? size() : newmax;

      int* lastfree = &firstfree;
      while (*lastfree != -themax - 1)
         lastfree = &(theitem[ -1 - *lastfree].info);
      *lastfree = -newmax - 1;

      themax = newmax;

      theitem = (Item*)realloc(theitem, themax * sizeof(Item));
      thekey = (Key *)realloc(thekey, themax * sizeof(Key));
      if (theitem == 0 || thekey == 0)
      {
         std::cerr << "ERROR: DataSet could not reallocate memory\n";
         exit(-1);
      }

      return long(theitem) - delta;
   }

   /** Assignment operator.
       The assignment operator involves #reMax()#ing the lvalue #DataSet#
       to the size needed for copying all elements of the rvalue. After the
       assignment all #Key#s from the lvalue are valid for the rvalue as
       well. They refer to a copy of the corresponding data elements.
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
      theitem = (Item*)malloc(themax * sizeof(Item));
      thekey = (Key *)malloc(themax * sizeof(Key));
      if (theitem == 0 || thekey == 0)
      {
         std::cerr << "ERROR: DataSet could not allocate memory\n";
         exit(-1);
      }
      assert(theitem);
      assert(thekey);
      memcpy(theitem, old.theitem, themax * sizeof(Item));
      memcpy(thekey, old.thekey, themax * sizeof(Key));
   }

   /// default constructor.
   DataSet(int max = 8)
   {
      themax = (max < 1) ? 8 : max;
      thesize = thenum = 0;
      firstfree = -themax - 1;
      theitem = (Item*)malloc(themax * sizeof(Item));
      thekey = (Key *)malloc(themax * sizeof(Key));
      if (theitem == 0 || thekey == 0)
      {
         std::cerr << "ERROR: DataSet could not allocate memory\n";
         exit(-1);
      }
      assert(theitem);
      assert(thekey);
   }

   /// destructor.
   ~DataSet()
   {
      free(theitem);
      free(thekey);
   }

   /// consistencty check.
   int isConsistent() const
   {
      if (theitem == 0 || thekey == 0)
      {
         std::cerr << "Inconsistency detected in class DataSet\n";
         return 0;
      }

      if (thesize > themax || thenum > themax || thenum > thesize)
      {
         std::cerr << "Inconsistency detected in class DataSet\n";
         return 0;
      }

      if (thesize == thenum && firstfree != -themax - 1)
      {
         std::cerr << "Inconsistency detected in class DataSet\n";
         return 0;
      }

      if (thesize != thenum && firstfree == -themax - 1)
      {
         std::cerr << "Inconsistency detected in class DataSet\n";
         return 0;
      }

      for (int i = 0; i < thenum; ++i)
      {
         if (theitem[thekey[i].idx].info != i)
         {
            std::cerr << "Inconsistency detected in class DataSet\n";
            return 0;
         }
      }
      return 1;
   }
   //@}
};

} // namespace soplex
#endif // _DATASET_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
