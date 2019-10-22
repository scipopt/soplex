/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  dataarray.h
 * @brief Save arrays of data objects.
 */
#ifndef _DATAARRAY_H_
#define _DATAARRAY_H_

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <iostream>

#include "soplex/spxdefines.h"
#include "soplex/spxalloc.h"

#include <boost/container/vector.hpp>

namespace soplex
{
/**@brief   Safe arrays of data objects.
   @ingroup Elementary

   Class DataArray provides safe arrays of \ref DataObjects. For general
   C++ objects (in contrast to data objects) class Array is provided which
   manages memory in a C++ compliant way.

   The elements of an instance of DataArray can be accessed just like
   ordinary C++ array elements by means of the index operator[](). Safety is
   provided by

    - automatic memory management in constructor and destructor
      preventing memory leaks
    - checking of array bounds when accessing elements with the
      indexing operator[]() (only when compiled without \c -DNDEBUG).

   Moreover, #DataArray%s may easily be extended by #insert%ing or #append%ing
   elements to the DataArray or shrunken by \ref remove() "removing" elements.
   Method reSize(int n) resets the DataArray%s length to \p n thereby possibly
   appending elements or truncating the DataArray to the required size.

   A DataArray may be used as arguments for standard C functions requiring
   pointers through the use of get_ptr() and get_const_ptr().

   Internally, a DataArray object allocates a block of memory that fits up
   to max() elements, only size() of them are used. This makes extension
   and shrinking methods perform better.

   @see Array, \ref DataObjects "Data Objects"
*/
template <typename T>
class DataArray
{
private:
   // A boost container vector is used instead of std::vector. This is because
   // std::vector<bool> is not consistent with the rest of the std::vectors<T>.
   // Namely, the current code needs the operator[] for DataArray<bool>.
   boost::container::vector<T>  data;              ///< the array of elements

public:

   /// reference \p n 'th element.
   T& operator[](int n)
   {
      assert(n >= 0);
      return data[n];
   }

   /// reference \p n 'th const element.
   const T& operator[](int n) const
   {
      assert(n >= 0);
      return data[n];
   }

   /// reference last element.
   T& last()
   {
      return data.back();
   }
   /// reference last const element.
   const T& last() const
   {
      return data.back();
   }

   /// get a C pointer to the data.
   T* get_ptr()
   {
      return data.data();
   }
   /// get a const C pointer to the data.
   const T* get_const_ptr() const
   {
      return data.data();
   }

   /// append element \p t.
   void append(const T& t)
   {
      data.push_back(t);
   }
   /// append \p n elements with value \p t.
   void append(int n, const T& t)
   {
      data.insert(data.end(), n, t);
   }
   /// append \p n elements from \p t.
   void append(int n, const T t[])
   {
      data.insert(data.end(), t, t + n);
   }
   /// append all elements from \p t.
   void append(const DataArray<T>& t)
   {
      data.insert(data.end(), t.data.begin(), t.data.end());
   }

   /// insert \p n elements with value \p t before \p i 'the element.
   void insert(int i, int n, const T& t)
   {
      if(n > 0)
      {
         data.insert(data.begin() + i - 1, n, t);
      }
   }

   /// insert \p n elements from \p t before \p i 'the element.
   void insert(int i, int n, const T t[])
   {
      if(n > 0)
      {
         // Inserts the elements of t (using legacy iterators, i.e., pointers)
         // before data's i th position.
         data.insert(data.begin() + i - 1, t, t + n);
      }
   }

   /// insert all elements from \p t before \p i 'th element.
   void insert(int i, const DataArray<T>& t)
   {
      if(t.size())
      {
         data.insert(data.begin() + i - 1, t.data.begin(), t.data.end());
      }
   }

   /// remove \p m elements starting at \p n.
   void remove(int n = 0, int m = 1)
   {
      assert(n < size() && n >= 0);

      if(n + m < size())
        {
          data.erase(data.begin() + n, data.begin() + n + m);
        }
      else
        {
          data.erase(data.begin() + n, data.end());
        }
   }
   /// remove \p m last elements.
   void removeLast(int m = 1)
   {
      assert(m <= size() && m >= 0);
      // Erase the last m elements
      data.erase(data.end() - m, data.end());
   }
   /// remove all elements.
   void clear()
   {
      data.clear();
   }

   /// return nr. of elements.
   int size() const
   {
      return int(data.size());
   }

   /// reset size to \p newsize.
   /** Resizing a DataArray to less than the previous size, involves
       discarding its last elements. Resizing to a larger value involves
       adding uninitialized elements (similar to append()). If neccessary,
       also memory will be reallocated.
       @param newsize the new number of elements the array can hold.
    */
   void reSize(int newsize)
   {
      if(newsize > int(data.capacity()))
         reMax(newsize, newsize);
      else if(newsize < 0)
      {
         data.clear();
      }
      else
      {
         data.resize(newsize);
      }
   }

   /// return maximum number of elements.
   /** Even though the DataArray currently holds no more than size()
       elements, up to max() elements could be added without need to
       reallocated free store.
    */
   int max() const
   {
      return int(data.capacity());
   }

   /// reset maximum number of elements.
   /** The value of max() is reset to \p newMax thereby setting size()
       to \p newSize. However, if \p newSize has a value \c < \c 0 (as the
       default argument does) size() remains unchanged and max() is set
       to MIN(size(), newMax). Hence, calling reMax() without the
       default arguments, will reduce the memory consumption to a minimum.
       In no instance max() will be set to a value less than 1 (even if
       specified).

    */
   void reMax(int newMax = 1, int newSize = -1)
   {
      if(newSize >= 0)
      {
         data.resize(newSize);
      }

      if(newMax < newSize)
         newMax = newSize;

      if(newMax < 1)
         newMax = 1;

      if(newMax == int(data.capacity()))
         return;

      int themax = newMax;

      if(newSize <= 0)
      {
         /* no data needs to be copied so do a clean free and alloc */
         data.clear();
         data.resize(themax);
      }
      else
      {
         data.resize(themax);
      }
   }
   /// assignment operator
   DataArray& operator=(const DataArray& rhs)
   {
      if(this != &rhs)
      {
         reSize(rhs.size());
         data = rhs.data;
      }

      return *this;
   }

  // Move assignment for Dataarray
  DataArray& operator=(const DataArray&& rhs)
  {
    data = std::move(rhs.data);
    return *this;
  }

   /// copy constructor
   DataArray(const DataArray& old)
   {
      data.reserve(max());
      data = old.data;
   }

   /// default constructor.
   /** The constructor allocates an Array containing \p size uninitialized
       elements. The internal array is allocated to have \p max nonzeros,
       and the memory extension factor is set to \p fac.

       @param p_size number of uninitialised elements.
       @param p_max  maximum number of elements the array can hold.
    */
   explicit DataArray(int p_size = 0, int p_max = 0)
   {
      // p_size number of uninitialized elements.
      //
      // That is, the operator[] is valid on for 0...(p_size -1) whereas, even
      // though the internal array can handle 0...(p_max -1) elements, the
      // operator[] on beyond (p_size -1) would give an exception.
      data.reserve(p_max);
      data.resize(p_size);
      // The underlying array will have at least p_max number of elements,
      // though.
   }

  // The move constructors
  DataArray(DataArray&& other) noexcept: data(std::move(other.data))
  {
  }

   /// destructor
   ~DataArray()
   {
      ;
   }

   void push_back(const T& val)
   {
      data.push_back(val);
   }

   void push_back(T&& val)
   {
      data.push_back(val);
   }
};

} // namespace soplex
#endif // _DATAARRAY_H_
