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
#pragma ident "@(#) $Id: array.h,v 1.2 2001/11/06 23:31:00 bzfkocht Exp $"

#ifndef _ARRAY_H_
#define _ARRAY_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <stdlib.h>
#include <assert.h>

namespace soplex
{

//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/**
    Class #Array# provides safe arrays of arbitrary type. #Array# elements are
    accessed just like ordinary C++ array elements by means of the index
    #operator[]#. Safety is provided by
    \begin{itemize}
    \item       automatic memory management in constructor and destructure
                preventing memory leaks
    \item       checking of array bound when accessing elements with the
                indexing #operator[]# (only when compiled without #-DNDEBUG#).
    \end{itemize}
 
    Moreover, #Array#s may easily be extended by #insert#ing or #append#ing elements
    to the #Array# or shrunken by #remov#ing elements. Method #reSize(int n)# resets
    the #Array#s length to #n# thereby appending elements or truncating the
    #Array# to the required size.
 
    Finally, #Array#s may be used as arguments of standard C functions requiring
    pointers thru a cast operator.
 
    An #Array# is implemented in a C++ complient way with respect to how memory
    is managemed: Only operators #new# and #delete# are used for allocating
    memory. This involves some overhead for all methods effecting the length of
    an #Array#, i.e. all methods #insert#, #append#, #remove# and #reSize#. This
    involves allocating a new C++ array of the new size and copying all elements
    with the template parameters #operator=#.
 
    For this reason, it is not convenient to use class #Array#, if its elements
    are \Ref{Data Objects}. In this case use class \Ref{DataArray} instead.
 
    @see        DataArray, Data Objects
 */
template < class T >
class Array
{
protected:
   /*  For the implementation we simply store:
    */
   int num;                    // the length of array #data# and
   T *data;                  // the array of elements

public:
   /// reference #n#'th element.
   T& operator[](int n)
   {
      assert(n >= 0 && n < size());
      return data[n];
   }
   /// reference #n#'th element.
   const T& operator[](int n) const
   {
      assert(n >= 0 && n < size());
      return data[n];
   }

   /// Cast to C arrey.
   operator T* ()
   {
      return data;
   }
   /** Cast to constant C array.
       This method allows it to use an #Array<T># as argument for all
       functions expecting #T*# or #const T*#. However, when used, bound
       checking can no longer be performed in the function using #T*#.
    */
   operator const T* () const
   {
      return data;
   }

   /// append #n# uninitialized elements.
   void append(int n)
   {
      insert(size(), n);
   }
   /// append #n# elements from #tarray#.
   void append(int n, const T* tarray)
   {
      insert(size(), n, tarray);
   }
   /// append all elements from #tarray#.
   void append(const Array<T>& tarray)
   {
      insert(size(), tarray);
   }

   /// insert #n# uninitialized elements before #i#-th element.
   void insert(int i, int n)
   {
      assert(i <= size());
      if (n > 0)
      {
         int k;
         T *olddata = data;
         data = new T[size() + n];
         assert(data);
         if (size() > 0)
         {
            for (k = 0; k < i; ++k)
               data[k] = olddata[k];
            for (; k < size(); ++k)
               data[k + n] = olddata[k];
            delete[] olddata;
         }
         num += n;
      }
   }

   /// insert #n# elements from #tarray# before #i#-th element.
   void insert(int i, int n, const T* tarray)
   {
      insert(i, n);
      for (n--; n >= 0; --n)
         data[n + i] = tarray[n];
   }

   /// insert all elements from #tarray# before #i#-th element.
   void insert(int i, const Array<T>& tarray)
   {
      int n = tarray.size();
      insert(i, n);
      for (n--; n >= 0; --n)
         data[n + i] = tarray.data[n];
   }

   /// remove #m# elements starting at #n#.
   void remove(int n = 0, int m = 1)
   {
      assert(n >= 0 && m >= 0);
      if (m > 0 && n < size())
      {
         T *olddata = data;
         m -= (n + m <= size()) ? 0 : n + m - size();
         num -= m;
         if (num > 0)
         {
            int i;
            data = new T[num];
            for (i = 0; i < n; ++i)
               data[i] = olddata[i];
            for (; i < num; ++i)
               data[i] = olddata[i + m];
         }
         delete[] olddata;
      }
   }

   /// remove all elements.
   void clear()
   {
      if (num > 0)
      {
         num = 0;
         delete[] data;
      }
   }


   /// return nr. of elements.
   int size() const
   {
      return num;
   }

   /// reset nr. of elements.
   void reSize(int newsize)
   {
      if (newsize < size())
         remove(newsize, size() - newsize);
      else if (newsize > size())
         append(newsize - size());
   }

   /** Assignment operator.
       Assigning an rvalue #Array# to an lvalue #Array# involves resizing
       the lvalue to the rvalues #size()# and copying all elements via
       the #Array# element's assignment #operator=#.
    */
   Array<T>& operator=(const Array<T>& rhs)
   {
      reSize(rhs.size());
      for (int i = 0; i < size(); ++i)
         data[i] = rhs.data[i];
      return *this;
   }

   /** Default Constructor.
       The constructor allocates an #Array# of #n# uninitialized elements.
    */
   Array(int n = 0)
      : data(0)
   {
      assert(n >= 0);
      num = n;
      if (num > 0)
      {
         data = new T[num];
         assert(data);
      }
   }

   ///
   Array(const Array<T>& old)
      : num(old.num)
   {
      if (num > 0)
      {
         data = new T[num];
         assert(data);
         *this = old;
      }
   }

   ///
   ~Array()
   {
      if (num > 0)
         delete[] data;
   }

   ///
   int isConsistent() const
   {
      if (num < 0 || (num > 0 && data == 0))
      {
         std::cerr << "Inconsistency detected in class array\n";
         return 0;
      }
      return 1;
   }
};

} // namespace soplex
#endif

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
