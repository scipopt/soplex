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
#pragma ident "@(#) $Id: subsvector.h,v 1.2 2001/11/06 23:31:06 bzfkocht Exp $"


#ifndef _SUBSVECTOR_H_
#define _SUBSVECTOR_H_

//@ ----------------------------------------------------------------------------
/*      \Section{Imports}
    Import required system include files
 */
#include <assert.h>


/*  and class header files
 */

#include "svector.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */

/** Part of an #SVector#.
    Class #SubSVector# provides a means to references a subset of nonzeros of an
    existing #SVector# and have it appear similar to an #SVector# of its own
    right. However, the user is responsible for avoiding any problems. Most
    notably, if an #SVector# changes, this may corrupt any #SubSVector#
    referencing to it.
 */
class SubSVector
{
   friend Vector& Vector::multAdd(double x, const SubSVector& vec);

   SVector::Element* elem;
   int               num;
#ifndef NDEBUG
   SVector*          svec;
#endif

public:
   /**@name Modification */
   //@{
   /// switch n'th with 0'th nonzero.
   void toFront(int n);

   /// sort nonzero to increasing indices.
   void sort();
   //@}


   /**@name Inquiery */
   //@{
   /// number of used indeces.
   int size() const
   {
      return num;
   }
   ///
   int& size()
   {
      return num;
   }

   /// maximal index.
   int dim() const;

   /** Number of index #i#.
       Return the number of the first index #i#. If no index #i# is available
       in the #IdxSet#, -1 is returned. Otherwise, #index(number(i)) == i#
       hods.
    */
   int number(int i) const;

   /// get value to index #i#.
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

   /// get #n#-th nonzero element.
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

   /// get index of #n#-th nonzero.
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

   /// get value of #n#-th nonzero.
   double value(int n) const
   {
      assert(n >= 0 && n < size());
      return elem[n].val;
   }
   //@}


   /**@name Mathematical Operations */
   //@{
   /// eucledian norm.
   double length() const
   {
      return sqrt(length2());
   }

   /// squared eucledian norm.
   double length2() const;

   /// scale with #x#.
   SubSVector& operator*=(double x);

   /// inner product.
   double operator*(const Vector& w) const;
   //@}


   /**@name Miscellaneous */
   //@{
   ///
   friend std::ostream& operator<<(std::ostream& os, const SubSVector& v);

   ///
   SubSVector(SVector* sv = 0, int first = 0, int len = 0)
      : elem((sv && first < sv->max()) ? &sv->element(first) : 0)
      , num (len)
#ifndef NDEBUG
      , svec(sv)
#endif
   { assert(isConsistent()); }
   ///
   SubSVector(const SubSVector& old)
      : elem(old.elem)
      , num (old.num)
#ifndef NDEBUG
      , svec(old.svec)
#endif
   {
      assert(isConsistent());
   }

   /// check consistency.
   int isConsistent() const;
   //@}

};

inline Vector& Vector::multAdd(double x, const SubSVector& vec)
{
   assert(vec.dim() <= dim());
   Vector_MultAddSVector(val, x, vec.size(), vec.elem);
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
