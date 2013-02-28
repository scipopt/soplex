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

/**@file  dsvectorbase.h
 * @brief Dynamic sparse vectors.
 */
#ifndef _DSVECTORBASE_H_
#define _DSVECTORBASE_H_

#include <assert.h>

#include "svectorbase.h"

namespace soplex
{
template < class R > class VectorBase;
template < class S > class SSVectorBase;

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
   friend class SLinSolver;

private:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Data */
   //@{

   /// Memory.
   Nonzero<R>* theelem;

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Private helpers */
   //@{

   /// Allocate memory for \p n nonzeros.
   void allocMem(int n)
   {
      spx_alloc(theelem, n);
      theelem = new (theelem) Nonzero<R>[n]();
      setMem(n, theelem);
   }

   /// Ensure there is room for \p n new nonzeros.
   void makeMem(int n)
   {
      if( SVectorBase<R>::max() - SVectorBase<R>::size() < n )
         setMax(SVectorBase<R>::size() + n);
   }

   //@}

public:

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Construction, assignment, and destruction */
   //@{

   /// Default constructor.
   /** Creates a DSVectorBase ready to hold \p n nonzeros. However, the memory is automatically enlarged, if more
    *  nonzeros are added to the DSVectorBase.
    */
   explicit DSVectorBase<R>(int n = 8)
      : theelem(0)
   {
      allocMem((n < 1) ? 2 : n);

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   explicit DSVectorBase<R>(const SVectorBase<S>& old)
      : theelem(0)
   {
      allocMem(old.size());
      SVectorBase<R>::operator=(old);

      assert(isConsistent());
   }

   /// Copy constructor.
   DSVectorBase<R>(const DSVectorBase<R>& old)
      : SVectorBase<R>()
      , theelem(0)
   {
      allocMem(old.size());
      SVectorBase<R>::operator=(old);

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   DSVectorBase<R>(const DSVectorBase<S>& old)
      : SVectorBase<R>()
      , theelem(0)
   {
      allocMem(old.size());
      SVectorBase<R>::operator=(old);

      assert(isConsistent());
   }

   /// Copy constructor.
   template < class S >
   explicit DSVectorBase<R>(const VectorBase<S>& vec);

   /// Copy constructor.
   template < class S >
   explicit DSVectorBase<R>(const SSVectorBase<S>& old);

   /// Assignment operator.
   template < class S >
   DSVectorBase<R>& operator=(const SVectorBase<S>& vec)
   {
      if( this != &vec )
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
      if( this != &vec )
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
      if( this != (DSVectorBase<R>*)(&vec) )
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
   ~DSVectorBase<R>()
   {
      std::cout << "~DSVectorBase\n";

      if( theelem )
      {
         for( int i = SVectorBase<R>::max() - 1; i >= 0; i-- )
            theelem[i].~Nonzero<R>();

         spx_free(theelem);
      }
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Modification */
   //@{

   /// Append nonzeros of \p sv.
   template < class S >
   void add(const SVectorBase<S>& vec)
   {
      SVectorBase<R>::clear();
      makeMem(vec.size());
      SVectorBase<S>::add(vec);
   }

   /// Append one nonzero \p (i,v).
   void add(int i, R v)
   {
      makeMem(1);
      SVectorBase<R>::add(i, v);
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

      if( len == SVectorBase<R>::max() )
         return;

      Nonzero<R>* newmem = 0;

      /* allocate and initialize new memory */
      spx_alloc(newmem, len);
      newmem = new (newmem) Nonzero<R>[len]();

      for( int i = 0; i < siz; i++ )
         newmem[i] = theelem[i];

      /* free old memory */
      for( int i = SVectorBase<R>::max()-1; i >= 0; i-- )
         theelem[i].~Nonzero<R>();

      spx_free(theelem);

      /* assign new memory */
      theelem = newmem;
      setMem(len, theelem);
      SVectorBase<R>::set_size(siz);
   }

   //@}

   // ------------------------------------------------------------------------------------------------------------------
   /**@name Utilities */
   //@{

   /// Consistency check.
   bool isConsistent() const
   {
#ifdef ENABLE_CONSISTENCY_CHECKS
      if( theelem != 0 && SVectorBase<R>::mem() != theelem )
         return MSGinconsistent("DSVectorBase");
#endif

      return true;
   }

   //@}
};

} // namespace soplex
#endif // _DSVECTORBASE_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
