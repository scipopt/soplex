/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  dvectorbase.h
 * @brief Dynamic dense vectors.
 */
#ifndef _DVECTORBASE_H_
#define _DVECTORBASE_H_

#include <iostream>
#include <assert.h>

#include "spxdefines.h"
#include "spxalloc.h"
#include "vectorbase.h"

namespace soplex
{
template < class R > class SVectorBase;

/**@brief   Dynamic dense vectors.
 * @ingroup Algebra
 *
 *  Class DVectorBase is a derived class of VectorBase adding automatic memory management to such objects.  This allows
 *  to implement maths operations operator+() and operator-().  Further, it is possible to reset the dimension of a
 *  DVectorBase via method reDim().  However, this may render all references to values of a #reDim()%ed DVectorBase
 *  invalid.
 *
 *  For vectors that are often subject to reDim() it may be unconvenient to reallocate the required memory every time.
 *  Instead, an array of values of length memSize() is kept, where only the first dim() elements are used.  Initially,
 *  memSize() == dim().  However, if the dimension is increased, memSize() will be increased at least by a factor of 1.2
 *  to be prepared for future (small) #reDim()%s.  Finally, one can explicitly set memSize() with method reSize(), but
 *  not lower than dim().
 */
template < class R >
class DVectorBase : public VectorBase<R>
{
   template < class S > friend class DVectorBase;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   /// Length of array of values \ref soplex::DVectorBase::mem "mem"
   int memsize;

   /// Array of values.
   R* mem;

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction, destruction, and assignment */
   //@{

   /// Default constructor. \p dim is the initial dimension.
   explicit DVectorBase<R>(int dim = 0)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      memsize = (dim > 0) ? dim : 4;

      spx_alloc(mem, memsize);
      mem = new (mem) R[memsize]();

      VectorBase<R>::val = mem;
      VectorBase<R>::dimen = dim;

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   explicit DVectorBase<R>(const VectorBase<S>& old)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      VectorBase<R>::dimen = old.dim();
      memsize = VectorBase<R>::dimen;

      spx_alloc(mem, memsize);
      mem = new (mem) R[memsize]();

      VectorBase<R>::val = mem;
      *this = old;

      assert(isConsistent());
   }

   /// Copy constructor.
   DVectorBase<R>(const DVectorBase<R>& old)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      VectorBase<R>::dimen = old.dim();
      memsize = old.memsize;

      spx_alloc(mem, memsize);
      mem = new (mem) R[memsize]();

      VectorBase<R>::val = mem;
      *this = old;

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   DVectorBase<R>(const DVectorBase<S>& old)
      : VectorBase<R>(0, 0)
      , mem(0)
   {
      VectorBase<R>::dimen = old.dim();
      memsize = old.memsize;

      spx_alloc(mem, memsize);
      mem = new (mem) R[memsize]();

      VectorBase<R>::val = mem;
      *this = old;

      assert(isConsistent());
   }

   /// Assignment operator.
   DVectorBase<R>& operator=(const VectorBase<R>& vec)
   {
      if( (VectorBase<R>*)this != &vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DVectorBase<R>& operator=(const VectorBase<S>& vec)
   {
      if( (VectorBase<S>*)this != &vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   DVectorBase<R>& operator=(const DVectorBase<R>& vec)
   {
      if( (void*)this != (void*)&vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DVectorBase<R>& operator=(const DVectorBase<S>& vec)
   {
      if( this != (DVectorBase<R>*)&vec )
      {
         if( vec.dim() != VectorBase<R>::dim() )
            reDim(vec.dim());

         VectorBase<R>::operator=(vec);

         assert(isConsistent());
      }

      return *this;
   }

   /// Assignment operator.
   template < class S >
   DVectorBase<R>& operator=(const SVectorBase<S>& vec);

   /// Destructor.
   ~DVectorBase<R>()
   {
      if( mem != 0 )
      {
         for( int i = memsize-1; i >= 0; i-- )
            mem[i].~R();

         spx_free(mem);
      }
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Access and modification */
   //@{

   /// Returns \ref soplex::DVectorBase "DVectorBase"'s memory size.
   int memSize() const
   {
      return memsize;
   }

   /// Resets \ref soplex::DVectorBase "DVectorBase"'s dimension to \p newdim.
   void reDim(int newdim)
   {
      assert(memsize >= 0);

      if( newdim > memsize )
         reSize(int(newdim + 0.2 * memsize));

      for( int i = VectorBase<R>::dimen; i < newdim; i++ )
         mem[i] = 0;

      VectorBase<R>::dimen = newdim;
   }

   /// Resets \ref soplex::DVectorBase "DVectorBase"'s memory size to \p newsize.
   void reSize(int newsize)
   {
      assert(newsize > VectorBase<R>::dim());

      R* newmem = 0;

      /* allocate and initialize new memory */
      spx_alloc(newmem, newsize);
      newmem = new (newmem) R[newsize]();

      for( int i = 0; i < VectorBase<R>::dim(); i++ )
         newmem[i] = mem[i];


      /* free old memory */
      for( int i = memsize-1; i >= 0; i-- )
         mem[i].~R();

      spx_free(mem);

      /* assign new memory */
      mem = newmem;
      memsize = newsize;
      VectorBase<R>::val = mem;
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Utilities */
   //@{

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( VectorBase<R>::val != mem || VectorBase<R>::dimen > memsize || VectorBase<R>::dimen < 0 )
         return MSGinconsistent("DVectorBase");

      return VectorBase<R>::isConsistent();
#endif

      return true;
   }

   //@}
};

} // namespace soplex
#endif // _DVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
