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
#pragma ident "@(#) $Id: subsvector.h,v 1.8 2002/01/19 18:59:18 bzfkocht Exp $"


/**@file  subsvector.h
 * @brief Part of an #SVector.
 */
#ifndef _SUBSVECTOR_H_
#define _SUBSVECTOR_H_

#include <assert.h>

#include "real.h"
#include "svector.h"

namespace soplex
{
/**@brief   Part of an #SVector.
   @ingroup Algebra

   Class #SubSVector provides a means to references a subset of nonzeros of an
   existing #SVector and have it appear similar to an #SVector of its own
   right. However, the user is responsible for avoiding any problems. Most
   notably, if an #SVector changes, this may corrupt any #SubSVector
   referencing to it.
*/
class SubSVector
{
private:
   /// output operator.
   friend std::ostream& operator<<(std::ostream& os, const SubSVector& v);

   friend Vector& Vector::multAdd(Real x, const SubSVector& vec);

   const SVector::Element* elem;   ///< element array.
   int                     num;    ///< number of nonzero elements.
#ifndef NDEBUG
   const SVector*          svec;   ///< pointer to underlying #SVector.
#endif

public:
   ///
   int size() const
   {
      return num;
   }
   /// returns the maximal index.
   int dim() const;

   /// returns number of index \p i.
   /** Returns the number of the first index \p i. 
       If no index \p i is available
       in the #IdxSet, -1 is returned. Otherwise, #index(number(i)) == i
       holds.
    */
   int number(int i) const;

   /// gets value to index \p i.
   Real operator[](int i) const
   {
      int n = number(i);
      if (n >= 0)
         return elem[n].val;
      return 0;
   }
   /// returns the \p n 'th nonzero index/value-pair.
   const SVector::Element& element(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n];
   }
   /// returns the index of the \p n 'th nonzero.
   int index(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n].idx;
   }
   /// returns the value of the \p n 'th nonzero.
   Real value(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n].val;
   }

   /// returns eucledian norm.
   Real length() const
   {
      return sqrt(length2());
   }

   /// returns squared eucledian norm.
   Real length2() const;

   /// returns inner product with \p w.
   Real operator*(const Vector& w) const;

   /// consistency check.
   bool isConsistent() const;

   /// direct assignment.
   void assign(const SVector* sv, int first , int len)
   {
      assert(sv          != 0);
      assert(first       >= 0);
      assert(first + len < sv->max());

      elem = &sv->element(first);
      num  = len;
#ifndef NDEBUG
      svec = sv;
#endif
      assert(isConsistent());      
   }
   /// default constructor.
   explicit SubSVector(const SVector* sv = 0, int first = 0, int len = 0)
      : elem((sv && first < sv->max()) ? &sv->element(first) : 0)
      , num (len)
#ifndef NDEBUG
      , svec(sv)
#endif
   { 
      assert(isConsistent()); 
   }

   /// copy constructor.
   SubSVector(const SubSVector& old)
      : elem(old.elem)
      , num (old.num)
#ifndef NDEBUG
      , svec(old.svec)
#endif
   {
      assert(isConsistent());
   }
};

inline Vector& Vector::multAdd(Real x, const SubSVector& vec)
{
   assert(vec.dim() <= dim());

   for(int i = 0; i < vec.size(); i++)
      val[vec.elem[i].idx] += x * vec.elem[i].val;

   return *this;
}

} // namespace soplex
#endif // _SUBSVECTOR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
