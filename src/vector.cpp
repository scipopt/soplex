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
#pragma ident "@(#) $Id: vector.cpp,v 1.2 2001/11/06 23:31:07 bzfkocht Exp $"


/*  \Section{Complex Methods}
 */

#include <stdlib.h>
#include <iostream>
#include "vector.h"


#include "ssvector.h"
#include "subsvector.h"
#include "svector.h"

namespace soplex
{




/* \SubSection{Assignment Operators}
 */
Vector& Vector::operator=(const Vector& vec)
{
   assert(dim() == vec.dim());
   memcpy(val, vec.val, dimen*sizeof(double));
   return *this;
}

Vector& Vector::operator=(const SVector& vec)
{
   clear();
   assign(vec);
   return *this;
}

Vector& Vector::assign(const SVector& psv)
{
   for (int i = psv.size(); i-- > 0;)
      val[psv.index(i)] = psv.value(i);
   return *this;
}

Vector& Vector::operator+=(const Vector& vec)
{
   assert(dim() == vec.dim());
   for (int i = 0; i < dim(); ++i)
      val[i] += vec[i];
   return *this;
}

Vector& Vector::operator+=(const SVector& vec)
{
   for (int i = vec.size(); i > 0; --i)
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] += vec.value(i);
   }
   return *this;
}

Vector& Vector::operator+=(const SubSVector& vec)
{
   for (int i = vec.size(); i > 0; --i)
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] += vec.value(i);
   }
   return *this;
}

Vector& Vector::operator-=(const Vector& vec)
{
   assert(dim() == vec.dim());
   for (int i = 0; i < dim(); ++i)
      val[i] -= vec[i];
   return *this;
}

Vector& Vector::operator-=(const SVector& vec)
{
   for (int i = vec.size(); i--;)
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }
   return *this;
}

Vector& Vector::operator-=(const SubSVector& vec)
{
   for (int i = vec.size(); i--;)
   {
      assert(vec.index(i) >= 0);
      assert(vec.index(i) < dim());
      val[vec.index(i)] -= vec.value(i);
   }
   return *this;
}

Vector& Vector::operator*=(double x)
{
   for (int i = 0; i < dim(); ++i)
      val[i] *= x;
   return *this;
}


/* \SubSection{Products}
 */
double Vector::length() const
{
   return sqrt(length2());
}

double Vector::length2() const
{
   return *this * *this;
}

double Vector::maxAbs() const
{
   double x = 0;
   int n = dim();
   double* v = val;
   while (n--)
   {
      x = (*v > x) ? *v : ((-*v > x) ? -*v : x);
      v++;
   }
   return x;
}


/* \SubSection{Miscellaneous}
 */
std::ostream& operator<<(std::ostream& s, const Vector& vec)
{
   int i;
   s << '(';
   for (i = 0; i < vec.dim() - 1; ++i)
      s << vec[i] << ", ";
   s << vec[i] << ')';
   return s;
}

double Vector::operator*(const SVector& v) const
{
   assert(dim() >= v.dim());
   int i;
   double x = 0;
   for (i = v.size(); i-- > 0;)
      x += val[v.index(i)] * v.value(i);
   return x;
}

double Vector::operator*(const SubSVector& v) const
{
   assert(dim() >= v.dim());
   int i;
   double x = 0;
   for (i = v.size(); i-- > 0;)
      x += val[v.index(i)] * v.value(i);
   return x;
}

int Vector::isConsistent() const
{
   if (dim() > 0 && val == 0)
   {
      std::cerr << "Inconsistency detected in class Vector\n";
      return 0;
   }
   return 1;
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
