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
#pragma ident "@(#) $Id: dataarray.h,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

#ifndef _DATAARRAY_H_
#define _DATAARRAY_H_


//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>
#include <stdlib.h>
#include <memory.h>
#include <iostream>

namespace soplex
{

//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/**@name Data Objects
    {\em Data Objects} refer to C++ objects that do not allocate any resources,
    particularly that do not allocate any memory.  This makes them behave just
    like ordinary C structures, in that both, the copy constructor and
    assignment operator are equivalent to a memcopy of the memory taken by the
    object. Examples for data objects are all builtin types such as #int# or
    #double# or ``simple'' classes such as #complex#.
 
    We distinguish {\em data objects} from general C++ objects that may include
    some allocation of resources. (Note, that for general C++ objects that do
    allocate resources, this must be respected by providing apropriate copy
    constructor and assignment operators.) An example for a general C++ class is
    class \Ref{DataArray}.
 
    The distingtion between data and general C++ objects becomes relevant when
    using such objects in container classes such as \Ref{DataArray} or
    \Ref{Array}.
 */

/** Safe arrays of data objects.
    Class #DataArray# provides safe arrays of \Ref{Data Objects}. For general
    C++ objects (i.e. no data objects) class \Ref{Array} is provided which
    manages memory in a C++ compliant way.
 
    The elements of an instance of #DataArray# can be accessed just like
    ordinary C++ array elements by means of the index #operator[]#. Safety is
    provided by
    \begin{itemize}
    \item       automatic memory management in constructor and destructor
                preventing memory leaks
    \item       checking of array bounds when accessing elements with the
                indexing #operator[]# (only when compiled without #-DNDEBUG#).
    \end{itemize}
 
    Moreover, #DataArray#s may easily be extended by #insert#ing or #append#ing
    elements to the #DataArray# or shrunken by #remov#ing elements. Method
    #reSize(int n)# resets the #DataArray#s length to #n# thereby possibly
    appending elements or truncating the #DataArray# to the required size.
 
    A #DataArray#s may be used as arguments for standard C functions requiring
    pointers through a cast operator.
 
    Internally, a #DataArray# objects allocates a block of memory that fits up
    to #max()# elements, only #size()# of them are used. This makes extension
    and shrinking methods perform better.
    @see        Array, Data Objects
 */
template < class T >
class DataArray
{
   /*  For the implementation we simply store:
    */
   int thesize;                // nr. of used elements in array #data# and
   int themax;                 // the length of array #data# and
   T *data;                  // the array of elements

public:
   /// reference #n#'th element.
   T& operator[](int n)
   {
      assert(n >= 0);
      assert(n < thesize);
      return data[n];
   }
   /// reference #n#'th const element.
   const T& operator[](int n) const
   {
      assert(n >= 0);
      assert(n < thesize);
      return data[n];
   }

   /// reference last element.
   T& last()
   {
      assert(thesize > 0);
      return data[thesize -1];
   }
   /// reference last const element.
   const T& last() const
   {
      assert(thesize > 0);
      return data[thesize -1];
   }

   /// cast to C array.
   T* get_ptr()
   {
      return data;
   }
   const T* get_const_ptr() const
   {
      return data;
   }

   //           operator T* ()                          { return data; }
   /// cast to const C array for const object.
   //           operator const T* () const              { return data; }
   /// cast to const C array
   //                   operator const T* ()                    { return data; }

   /// append element #t#.
   long append(const T& t)
   {
      return insert(thesize, 1, &t);
   }
   /// append #n# elements from #t#.
   long append(int n, const T t[])
   {
      return insert(thesize, n, t);
   }
   /// append all elements from #t#.
   long append(const DataArray<T>& t)
   {
      return insert(thesize, t);
   }

   /// insert #n# uninitialized elements before #i#-th element.
   long insert(int i, int n)
   {
      long j = size();
      long delta = reSize(thesize + n);
      while (i < j--)
         data[j + n] = data[j];
      return delta;
   }

   /// insert #n# elements from #t# before #i#-the element.
   long insert(int i, int n, const T t[])
   {
      if (n > 0)
      {
         long delta = insert(i, n);
         memcpy(&(data[i]), t, n*sizeof(T));
         return delta;
      }
      return 0;
   }

   /// insert all elements from #t# before #i#-the element.
   long insert(int i, const DataArray<T>& t)
   {
      if (t.size())
      {
         long delta = insert(i, t.size());
         memcpy(&(data[i]), t.data, t.size()*sizeof(T));
         return delta;
      }
      return 0;
   }

   /// remove #m# elements starting at #n#.
   void remove(int n = 0, int m = 1)
   {
      assert(n < size() && n >= 0);
      if (n + m < size())
         memcpy(&(data[n]), &(data[n + m]), (size() - (n + m)) * sizeof(T));
      else
         m = size() - n;
      thesize -= m;
   }
   /// remove #m# last elements.
   void removeLast(int m = 1)
   {
      assert(m <= size() && m >= 0);
      thesize -= m;
   }
   /// remove all elements.
   void clear()
   {
      thesize = 0;
   }

   /// return nr. of elements.
   int size() const
   {
      return thesize;
   }

   /** reset size to #newsize#.
       Resizing a #DataArray# to less than the previous size, involves
       discarding its last elements. Resizing to a larger value involves
       adding uninitialized elements (similar to #append#). If neccessary,
       also memory will be reallocated.
    */
   long reSize(int newsize)
   {
      assert(memFactor >= 1);
      if (newsize > themax)
         return reMax(int(memFactor * newsize), newsize);
      else if (newsize < 0)
         thesize = 0;
      else
         thesize = newsize;
      return 0;
   }

   /** return maximum nr. of elements.
       Even though the #DataArray# currently holds no more than #size()#
       elements, up to #max()# elements could be added without need to
       reallocated free store.
    */
   int max() const
   {
      return themax;
   }

   /** reset maximum nr. of elements.
       The value of #max()# is reset to #newMax# thereby setting #size()#
       to #newSize#. However, if #newSize# has a value #< 0# (as the
       default argument does) #size()# remains unchanged and #max()# is set
       to #MIN(size(), newMax)#. Hence, calling #reMax()# without the
       default arguments, will reduce the memory consumption to a minimum.
       In no instance #max()# will be set to a value less than 1 (even if
       specified).
    */
   long reMax(int newMax = 1, int newSize = -1)
   {
      if (newSize >= 0)
         thesize = newSize;
      if (newMax < newSize)
         newMax = newSize;
      if (newMax < 1)
         newMax = 1;
      if (newMax == themax)
         return 0;
      themax = newMax;
      long olddata = long(data);
      if (thesize <= 0)
      {
         free(data);
         data = (T*)malloc(themax * sizeof(T));
      }
      else
         data = (T*)realloc(data, themax * sizeof(T));
      if (data == 0)
      {
         std::cerr << "ERROR: DataArray could not reallocate memory\n";
         exit(-1);
      }
      assert(data);
      return long(data) - olddata;
   }
   /** memory extension factor.
       When a #DataArray# is #reSize()#d to more than #max()# elements, the
       new value for #max()# is not just set to the new size but rather to
       #memFactor * size#. This makes #reSize#ing perform better in codes
       where a #DataArray# is extended often by a small number of elements
       only.
    */
   double memFactor;

   ///
   DataArray& operator=(const DataArray& rhs)
   {
      reSize(rhs.size());
      memcpy(data, rhs.data, size() * sizeof(T));
      return *this;
   }

   ///
   int isConsistent() const
   {
      if (data == 0
           || themax < 1
           || themax < thesize)
      {
         std::cout << "Inconsistency detected in class DataArray\n";
         return 0;
      }
      return 1;
   }


   ///
   DataArray(const DataArray& old)
      : thesize(old.thesize)
         , themax (old.themax)
         , memFactor (old.memFactor)
   {
      data = (T*)malloc(max() * sizeof(T));
      if (data == 0)
      {
         std::cerr << "ERROR: DataArray could not allocate memory\n";
         exit(-1);
      }
      if (thesize)
         memcpy(data, old.data, thesize * sizeof(T));
      assert(isConsistent());
   }

   /** Default constructor.
       The constructor allocates an #Array# containing #size# uninitialized
       elements. The internal array is allocated to have #max# nonzeros,
       and the memory extension factor is set to #fac#.
    */
   DataArray(int size = 0, int max = 0, double fac = 1.2)
      : memFactor (fac)
   {
      thesize = (size < 0) ? 0 : size;
      if (max > thesize)
         themax = max;
      else
         themax = (thesize == 0) ? 1 : thesize;
      data = (T*)malloc(themax * sizeof(T));
      if (data == 0)
      {
         std::cerr << "ERROR: DataArray could not allocate memory\n";
         exit(-1);
      }
      assert(isConsistent());
   }

   ///
   ~DataArray()
   {
      free(data);
   }
};

} // namespace soplex
#endif

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
