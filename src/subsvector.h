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
#pragma ident "@(#) $Id: subsvector.h,v 1.5 2001/12/28 14:55:13 bzfkocht Exp $"


/**@file  subsvector.h
 * @brief Part of an #SVector.
 */
#ifndef _SUBSVECTOR_H_
#define _SUBSVECTOR_H_

#include <assert.h>

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
   friend Vector& Vector::multAdd(double x, const SubSVector& vec);

   SVector::Element* elem;   ///< element array.
   int               num;    ///< number of nonzero elements.
#ifndef NDEBUG
   SVector*          svec;   ///< pointer to underlying #SVector.
#endif

public:
   /**@name Modification */
   //@{
   /// switches \p n 'th with 0'th nonzero.
   void toFront(int n);

   /// sorts nonzeros to increasing indices.
   void sort();
   //@}


   /**@name Inquiry */
   //@{
   ///
   int size() const
   {
      return num;
   }
   /// returns the maximal index.
   int dim() const;

   /// returns number of index \p i.
   /** Returns the number of the first index \p i. If no index \p i is available
       in the #IdxSet, -1 is returned. Otherwise, #index(number(i)) == i
       holds.
    */
   int number(int i) const;

   /// gets value to index \p i.
   double operator[](int i) const
   {
      int n = number(i);
      if (n >= 0)
         return elem[n].val;
      return 0;
   }

   ///
   SVector::Element& element(int n)
   {
      assert(n >= 0 && n < size());
      return elem[n];
   }
   /// returns the \p n 'th nonzero index/value-pair.
   SVector::Element element(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n];
   }

   ///
   int& index(int n)
   {
      assert(n >= 0 && n < size());
      return elem[n].idx;
   }
   /// returns the index of the \p n 'th nonzero.
   int index(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n].idx;
   }

   ///
   double& value(int n)
   {
      assert(n >= 0 && n < size());
      return elem[n].val;
   }
   /// returns the value of the \p n 'th nonzero.
   double value(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n].val;
   }
   //@}


   /**@name Mathematical Operations */
   //@{
   /// returns eucledian norm.
   double length() const
   {
      return sqrt(length2());
   }

   /// returns squared eucledian norm.
   double length2() const;

   /// scales vector with \p x.
   SubSVector& operator*=(double x);

   /// returns inner product with \p w.
   double operator*(const Vector& w) const;
   //@}


   /**@name Miscellaneous */
   //@{
   /// consistency check.
   int isConsistent() const;
   //@}

   
   /**@name Constructors / Destructors */
   //@{
   /// default constructor.
   SubSVector(SVector* sv = 0, int first = 0, int len = 0)
      : elem((sv && first < sv->max()) ? &sv->element(first) : 0)
      , num (len)
#ifndef NDEBUG
      , svec(sv)
#endif
   { assert(isConsistent()); }

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
   //@}

   /// output operator.
   friend std::ostream& operator<<(std::ostream& os, const SubSVector& v);
};

inline Vector& Vector::multAdd(double x, const SubSVector& vec)
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
