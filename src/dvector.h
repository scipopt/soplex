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
#pragma ident "@(#) $Id: dvector.h,v 1.2 2001/11/06 23:31:01 bzfkocht Exp $"

#ifndef _DVECTOR_H_
#define _DVECTOR_H_

/*      \Section{Imports}
    Import required system include files
 */
#include <iostream>
#include <stdlib.h>
#include <assert.h>

#include "vector.h"
#include "svector.h"

namespace soplex
{

//@ ----------------------------------------------------------------------------
/* \Section{Class Declaration}
 */
/** Dynamic vectors.
    Class #DVector# is a derived class of #Vector# adding automatic memory
    management to such objects. This allows to implement maths operations
    #operator+# and #operator-#. Further, it is possible to reset the dimension
    of a #DVector#s via method #reDim()#. However, this may render all
    references to values of a #reDim()#ed #DVector# invalid.
 
    For vectors that are often subject to #reDim()# it may be
    unconvenient to reallocate the required memory every
    time. Instead, an array of values of length #memSize()# is kept,
    where only the first #dim()# elements are used.  Initially,
    #memSize() == dim()#. However, if the dimension is increased,
    #memSize()# will be increased at least by a factor of 1.2 to be
    prepared for futur (small) #reDim()#s. Finally, one can explicitly
    set #memSize()# with method #reSize()#, but not lower than
    #dim()#.
*/
class DVector : public Vector
{
   int memsize;        // length of array of values #mem#
   double* mem;            // value array to be used

public:
   /**@name Maths */
   //@{
   ///
   friend DVector operator+(const Vector& v, const Vector& w);
   ///
   friend DVector operator+(const SVector& v, const Vector& w);
   /// adding vectors.
   friend DVector operator+(const Vector& v, const SVector& w);

   ///
   friend DVector operator-(const Vector& v, const Vector& w);
   ///
   friend DVector operator-(const SVector& v, const Vector& w);
   /// subtracting vectors.
   friend DVector operator-(const Vector& v, const SVector& w);

   ///
   friend DVector operator-(const Vector& vec);
   /// negative vectors.
   friend DVector operator-(const SVector& vec);

   ///
   friend DVector operator*(const Vector& v, double x);
   /// scaled vectors.
   friend DVector operator*(double x, const Vector& v);
   //@}


   /**@name Assignments */
   //@{
   ///
   DVector& operator+=(const Vector& vec)
   {
      Vector::operator+=(vec);
      return *this;
   }

   ///
   DVector& operator+=(const SVector& vec)
   {
      Vector::operator+=(vec);
      return *this;
   }

   ///
   DVector& operator-=(const Vector& vec)
   {
      Vector::operator-=(vec);
      return *this;
   }

   ///
   DVector& operator-=(const SVector& vec)
   {
      Vector::operator-=(vec);
      return *this;
   }

   ///
   DVector& operator*=(double x)
   {
      Vector::operator*=(x);
      return *this;
   }

   ///
   DVector& operator=(const Vector& vec)
   {
      if (vec.dim() != dim())
         reDim(vec.dim());
      Vector::operator=(vec);
      return *this;
   }

   /// assingment operator
   //lint -e1529 Test for self assignment is not neccessary here.
   DVector& operator=(const DVector& vec)
   {
      if (vec.dim() != dim())
         reDim(vec.dim());
      Vector::operator=(vec);

      return *this;
   }
   ///
   DVector& operator=(const SVector& vec)
   {
      if (vec.dim() != dim())  // ??? TK inserted here. I am not sure
         reDim(vec.dim());    // ??? TK if it is really needed

      Vector::operator=(vec);
      return *this;
   }

   //@}

   /**@name Miscellaneous */
   //@{
   /// reset #DVector#s dimension to #newdim#.
   void reDim(int newdim);
   /// reset #DVector#s memory size to #newsize#.
   void reSize(int newsize);
   ///  reset #DVector#s memory size to #newsize# and dimension to #newdim#
   void reSize(int newsize, int newdim);
   /// get #DVector#s memory size.
   int memSize() const
   {
      return memsize;
   }

   ///
   friend std::istream& operator>>(std::istream& s, DVector& vec);
   ///
   DVector(const Vector& old);
   ///
   DVector(const DVector& old);
   ///
   DVector(int dim = 0);
   ///
   ~DVector();
   ///
   int isConsistent() const;
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
