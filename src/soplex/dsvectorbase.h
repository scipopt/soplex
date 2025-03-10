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

/**@file  dsvectorbase.h
 * @brief Dynamic sparse vectors.
 */
#ifndef _DSVECTORBASE_H_
#define _DSVECTORBASE_H_

#include <assert.h>

#include "soplex/svectorbase.h"

namespace soplex
{
template < class R > class VectorBase;
template < class S > class SSVectorBase;
template < class R > class SLinSolver;

/**@brief   Dynamic sparse vectors.
 * @ingroup Algebra
 *
 *  Class DSVectorBase implements dynamic sparse vectors, i.e. #SVectorBase%s with an automatic memory management. This
 *  allows the user to freely add() as many nonzeros to a DSVectorBase as desired, without any precautions.  For saving
 *  memory method setMax() allows to reduce memory consumption to the amount really required.
 *
 *  @todo Both DSVectorBase and SVectorBase have a member variable that points to allocated memory. This does not seem to
 *        make too much sense.  Why doesn't DSVectorBase use the element of its base class?
 */
template < class R >
class DSVectorBase : public SVectorBase<R>
{
   friend class SLinSolver<R>;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   ///@{

   /// Memory.
   Nonzero<R>* theelem;

   ///@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Private helpers */
   ///@{

   /// Allocate memory for \p n nonzeros.
   void allocMem(int n)
   {
      spx_alloc(theelem, n);

      for(int i = 0; i < n; i++)
         new(&(theelem[i])) Nonzero<R>();

      SVectorBase<R>::setMem(n, theelem);
   }

   /// Ensure there is room for \p n new nonzeros.
   void makeMem(int n)
   {
      assert(n >= 0);

      if(SVectorBase<R>::max() - SVectorBase<R>::size() < n)
      {
         assert(SVectorBase<R>::size() + n > 0);
         setMax(SVectorBase<R>::size() + n);
      }
   }

   ///@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction, assignment, and destruction */
   ///@{

   /// Default constructor.
   /** Creates a DSVectorBase ready to hold \p n nonzeros. However, the memory is automatically enlarged, if more
    *  nonzeros are added to the DSVectorBase.
    */
   explicit DSVectorBase(int n = 8)
      : theelem(nullptr)
   {
      allocMem((n < 1) ? 2 : n);

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   explicit DSVectorBase(const SVectorBase<S>& old)
      : theelem(nullptr)
   {
      allocMem(old.size());
      SVectorBase<R>::operator=(old);

      assert(isConsistent());
   }

   /// Copy constructor.
   /** The redundancy with the copy constructor below is necessary since otherwise the compiler doesn't realize that it
    *  could use the more general one with S = R and generates a shallow copy constructor.
    */
   DSVectorBase(const DSVectorBase<R>& old)
      : SVectorBase<R>()
      , theelem(nullptr)
   {
      allocMem(old.size());
      SVectorBase<R>::operator=(old);

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   DSVectorBase(const DSVectorBase<S>& old)
      : SVectorBase<R>()
      , theelem(nullptr)
   {
      allocMem(old.size());
      SVectorBase<R>::operator=(old);

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   explicit DSVectorBase(const VectorBase<S>& vec);

   /// Copy constructor.
   template < class S >
   explicit DSVectorBase(const SSVectorBase<S>& old);

   /// Assignment operator.
   template < class S >
   DSVectorBase<R>& operator=(const SVectorBase<S>& vec)
   {
      if(this != &vec)
      {
         SVectorBase<R>::clear();
         makeMem(vec.size());
         SVectorBase<R>::operator=(vec);
      }

      return *this;
   }

   /// Assignment operator.
   DSVectorBase<R>& operator=(const DSVectorBase<R>& vec)
   {
      if(this != &vec)
      {
         SVectorBase<R>::clear();
         makeMem(vec.size());
         SVectorBase<R>::operator=(vec);
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DSVectorBase<R>& operator=(const DSVectorBase<S>& vec)
   {
      if(this != (DSVectorBase<R>*)(&vec))
      {
         SVectorBase<R>::clear();
         makeMem(vec.size());
         SVectorBase<R>::operator=(vec);
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DSVectorBase<R>& operator=(const VectorBase<S>& vec);

   /// Assignment operator.
   template < class S >
   DSVectorBase<R>& operator=(const SSVectorBase<S>& vec);

   /// Destructor.
   virtual ~DSVectorBase()
   {
      if(theelem)
      {
         for(int i = SVectorBase<R>::max() - 1; i >= 0; i--)
            theelem[i].~Nonzero<R>();

         spx_free(theelem);
      }
   }

   ///@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Modification */
   ///@{

   /// Append nonzeros of \p sv.
   template < class S >
   void add(const SVectorBase<S>& vec)
   {
      SVectorBase<R>::clear();
      makeMem(vec.size());
      SVectorBase<S>::add(vec);
   }

   /// Append one nonzero \p (i,v).
   void add(int i, const R& v)
   {
      makeMem(1);
      SVectorBase<R>::add(i, v);
   }

   /// Append one uninitialized nonzero.
   void add(int i)
   {
      makeMem(1);
      SVectorBase<R>::add(i);
   }

   /// Append \p n nonzeros.
   void add(int n, const int i[], const R v[])
   {
      makeMem(n);
      SVectorBase<R>::add(n, i, v);
   }

   /// Reset nonzero memory to >= \p newmax.
   /** This methods resets the memory consumption to \p newmax. However, if \p newmax < size(), it is
    *  reset to size() only.
    */
   void setMax(int newmax = 1)
   {
      int siz = SVectorBase<R>::size();
      int len = (newmax < siz) ? siz : newmax;

      if(len == SVectorBase<R>::max())
         return;

      Nonzero<R>* newmem = nullptr;

      /* allocate new memory */
      spx_alloc(newmem, len);

      /* call copy constructor for first elements */
      int i;

      for(i = 0; i < siz; i++)
         new((&newmem[i])) Nonzero<R>(theelem[i]);

      /* call default constructor for remaining elements */
      for(; i < len; i++)
         new((&newmem[i])) Nonzero<R>();

      /* free old memory */
      for(i = SVectorBase<R>::max() - 1; i >= 0; i--)
         theelem[i].~Nonzero<R>();

      if(theelem != nullptr)
         spx_free(theelem);

      /* assign new memory */
      theelem = newmem;
      SVectorBase<R>::setMem(len, theelem);
      SVectorBase<R>::set_size(siz);
   }

   ///@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Utilities */
   ///@{

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS

      if(theelem != nullptr && SVectorBase<R>::mem() != theelem)
         return SPX_MSG_INCONSISTENT("DSVectorBase");

#endif

      return true;
   }

   ///@}
};



/// Allocate memory for \p n nonzeros (specialization for Real).
template<>
inline
void DSVectorBase<Real>::allocMem(int n)
{
   spx_alloc(theelem, n);
   SVectorBase<Real>::setMem(n, theelem);
}



/// Destructor (specialization for Real).
template<>
inline
DSVectorBase<Real>::~DSVectorBase()
{
   if(theelem)
      spx_free(theelem);
}



/// Reset nonzero memory to >= \p newmax.
/** This methods resets the memory consumption to \p newmax. However, if \p newmax < size(), it is
 *  reset to size() only (specialization for Real).
 */
template<>
inline
void DSVectorBase<Real>::setMax(int newmax)
{
   int siz = size();
   int len = (newmax < siz) ? siz : newmax;

   spx_realloc(theelem, len);
   setMem(len, theelem);
   // reset 'size' to old size since the above call to setMem() sets 'size' to 0
   set_size(siz);
}
} // namespace soplex
#endif // _DSVECTORBASE_H_
