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
#pragma ident "@(#) $Id: ssvector.cpp,v 1.8 2001/12/25 17:00:09 bzfkocht Exp $"

#include <assert.h>

#include "ssvector.h"
#include "svset.h"
#include "spxmessage.h"

/**@file ssvector.cpp
 * @todo There is a lot pointer arithmetic done here. It is not clear if
 *       this is an advantage at all. See all the function int() casts.
 */

namespace soplex
{

static const double shortProductFactor = 0.5;

void SSVector::reDim (int newdim)
{
   for (int i = DIdxSet::size() - 1; i >= 0; --i)
      if (index(i) >= newdim)
         remove(i);
   DVector::reDim(newdim);
   DIdxSet::setMax(DVector::memSize() + 1);
   assert(isConsistent());
}

void SSVector::reMem(int newsize)
{
   DVector::reSize(newsize);
   assert(isConsistent());
   DIdxSet::setMax(DVector::memSize() + 1);
}

void SSVector::clear ()
{
   if (isSetup())
   {
      int i = DIdxSet::size();
      int* iptr;
      for (iptr = idx; i--; iptr++)
         val[*iptr] = 0;
   }
   else
      Vector::clear();

   IdxSet::clear();
   setupStatus = 1;
   assert(isConsistent());
}

void SSVector::setValue(int i, double x)
{
   assert(i >= 0 && i < DVector::dim());

   if (isSetup())
   {
      int n = number(i);

      if (n < 0)
      {
         if (x > epsilon || x < -epsilon)
            IdxSet::add(1, &i);
      }
      else if (x == 0)
         clearNum(n);
   }
   val[i] = x;

   assert(isConsistent());
}

void SSVector::setup()
{
   if (!isSetup())
   {
      IdxSet::clear();

      // #define      TWO_LOOPS
#ifdef  TWO_LOOPS
      int i = 0;
      int n = 0;
      int* id = idx;
      double* v = val;
      const double* end = val + dim();

      while (v < end)
      {
         id[n] = i++;
         n += (*v++ != 0);
      }

      double x;
      int* ii = idx;
      int* last = idx + n;
      const double eps = epsilon;
      const double meps = -eps;
      v = val;

      for (; id < last; ++id)
      {
         x = v[*id];
         if (x > eps || x < meps)
            *ii++ = *id;
         else
            v[*id] = 0;
      }
      num = ii - idx;

#else

      if (dim() <= 1)
      {
         if (dim())
         {
            if (*val > epsilon || *val < -epsilon)
               IdxSet::add(0);
            else
               *val = 0;
         }
      }
      else
      {
         int* ii = idx;
         double* v = val;
         double* end = v + dim() - 1;
         const double eps = epsilon;
         const double meps = -eps;

         /* setze weissen Elefanten */
         double last = *end;
         *end = 1e-100;

         /* erstes element extra */
         if (*v > eps || *v < meps)
            *ii++ = 0;
         else
            *v = 0;

         for(;;)
         {
            while (!*++v);
            if (*v > eps || *v < meps)
            {
               *ii++ = int(v - val);
            }
            else
            {
               *v = 0;
               if (v == end)
                  break;
            }

         }

         /* fange weissen Elefanten wieder ein */
         if (last > eps || last < meps)
         {
            *v = last;
            *ii++ = dim() - 1;
         }
         else
            *v = 0;

         num = int(ii - idx);
      }

#endif

      setupStatus = 1;
      assert(isConsistent());
   }
}

SSVector& SSVector::operator+=(const Vector& vec)
{
   Vector::operator+=(vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator+=(const SVector& vec)
{
   Vector::operator+=(vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator+=(const SubSVector& vec)
{
   Vector::operator+=(vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator+=(const SSVector& vec)
{
   for (int i = vec.size() - 1; i >= 0; --i)
      val[vec.index(i)] += vec.value(i);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const Vector& vec)
{
   Vector::operator-=(vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const SVector& vec)
{
   Vector::operator-=(vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const SubSVector& vec)
{
   Vector::operator-=(vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator-=(const SSVector& vec)
{
   if (vec.isSetup())
   {
      for (int i = vec.size() - 1; i >= 0; --i)
         val[vec.index(i)] -= vec.value(i);
   }
   else
   {
      Vector::operator-=(Vector(vec));
   }

   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }

   return *this;
}

SSVector& SSVector::operator*=(double x)
{
   for (int i = size() - 1; i >= 0; --i)
      val[index(i)] *= x;
   assert(isConsistent());
   return *this;
}

double SSVector::maxAbs() const
{
   if (isSetup())
   {
      double x;
      int* i = idx;
      int* end = idx + num;
      double* v = val;
      double abs = 0;

      for (; i < end; ++i)
      {
         x = v[*i];
         if (x > abs)
            abs = x;
         else if (-x > abs)
            abs = -x;
      }
      return abs;
   }
   else
      return Vector::maxAbs();
}

double SSVector::length2() const
{
   if (isSetup())
   {
      int* i = idx;
      int* end = idx + num;
      double* v = val;
      double x = 0;

      for (; i < end; ++i)
         x += v[*i] * v[*i];
      return x;
   }
   else
      return Vector::length2();
}

double SSVector::length() const
{
   return sqrt(length2());
}


SSVector& SSVector::multAdd(double xx, const SSVector& svec)
{
   if (svec.isSetup())
   {
      if (isSetup())
      {
         int i, j;
         double x;
         const double eps = epsilon;

         for (i = svec.size() - 1; i >= 0; --i)
         {
            j = svec.index(i);
            if (val[j])
            {
               x = val[j] + xx * svec.value(i);
               if (x > eps || x < -eps)
                  val[j] = x;
               else
               {
                  val[j] = 0;
                  for (--i; i >= 0; --i)
                     val[svec.index(i)] += xx * svec.value(i);
                  unSetup();
                  break;
               }
            }
            else
            {
               x = xx * svec.value(i);
               if (x > eps || x < -eps)
               {
                  val[j] = x;
                  addIdx(j);
               }
            }
         }
      }
      else
         Vector::multAdd(xx, svec);
   }
   else
   {
      double y;
      int* ii = idx;
      double* v = val;
      double* rv = static_cast<double*>(svec.val);
      double* last = rv + svec.dim() - 1;
      double x = *last;
      const double eps = epsilon;
      const double meps = -eps;

      *last = 1e-100;
      for(;;)
      {
         while (!*rv)
         {
            ++rv;
            ++v;
         }
         y = *rv++ * xx;
         if (y < meps || y > eps)
         {
            *ii++ = int(v - val);
            *v++ = y;
         }
         else if (rv == last)
            break;
         else
            v++;
      }
      *rv = x;

      x *= xx;
      if (x < meps || x > eps)
      {
         *ii++ = int(v - val);
         *v = x;
      }
      num = int(ii - idx);

      setupStatus = 1;
   }

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::multAdd(double xx, const SVector& svec)
{
   if (isSetup())
   {
      int i, j;
      double x;
      double* v = val;
      const double eps = epsilon;
      const double meps = -eps;
      const double mark = 1e-100;
      int adjust = 0;

      for (i = svec.size() - 1; i >= 0; --i)
      {
         j = svec.index(i);
         if (v[j])
         {
            x = v[j] + xx * svec.value(i);
            if (x > eps || x < meps)
               v[j] = x;
            else
            {
               adjust = 1;
               v[j] = mark;
            }
         }
         else
         {
            x = xx * svec.value(i);
            if (x > eps || x < meps)
            {
               v[j] = x;
               addIdx(j);
            }
         }
      }

      if (adjust)
      {
         int* iptr = idx;
         int* iiptr = idx;
         int* endptr = idx + num;
         for (; iptr < endptr; ++iptr)
         {
            x = v[*iptr];
            if (x > eps || x < meps)
               *iiptr++ = *iptr;
            else
               v[*iptr] = 0;
         }
         num = iiptr - idx;
      }
   }
   else
      Vector::multAdd(xx, svec);

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::multAdd(double xx, const SubSVector& svec)
{
   if (isSetup())
   {
      int i, j;
      double x;
      double* v = val;
      const double eps = epsilon;
      const double meps = -epsilon;
      const double mark = 1e-100;
      int adjust = 0;

      for (i = svec.size() - 1; i >= 0; --i)
      {
         j = svec.index(i);
         if (v[j])
         {
            x = v[j] + xx * svec.value(i);
            if (x > eps || x < meps)
               v[j] = x;
            else
            {
               adjust = 1;
               v[j] = mark;
            }
         }
         else
         {
            x = xx * svec.value(i);
            if (x > eps || x < meps)
            {
               v[j] = x;
               addIdx(j);
            }
         }
      }

      if (adjust)
      {
         int* iptr = idx;
         int* iiptr = idx;
         int* endptr = idx + num;
         for (; iptr < endptr; ++iptr)
         {
            x = v[*iptr];
            if (x > eps || x < meps)
               *iiptr++ = *iptr;
            else
               v[*iptr] = 0;
         }
         num = int(iiptr - idx);
      }
   }
   else
      Vector::multAdd(xx, svec);

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::multAdd(double x, const Vector& vec)
{
   Vector::multAdd(x, vec);
   if (isSetup())
   {
      setupStatus = 0;
      setup();
   }
   return *this;
}

SSVector& SSVector::operator=(const SSVector& rhs)
{
   clear();

   DIdxSet::setMax(rhs.max());
   DIdxSet::operator=(rhs);
   DVector::reDim(rhs.dim());

   if (rhs.isSetup())
   {
      int i, j;
      for (i = size() - 1; i >= 0; --i)
      {
         j = index(i);
         val[j] = rhs.val[j];
      }
   }
   else
   {
      int* ii = idx;
      double* v = val;
      double* rv = static_cast<double*>(rhs.val);
      double* last = rv + rhs.dim() - 1;
      double x = *last;
      const double eps = epsilon;
      const double meps = -eps;

      *last = 1e-100;
      for(;;)
      {
         while (!*rv)
         {
            ++rv;
            ++v;
         }
         if (*rv < meps || *rv > eps)
         {
            *ii++ = int(v - val);
            *v++ = *rv++;
         }
         else if (rv == last)
            break;
         else
         {
            v++;
            rv++;
         }
      }
      *rv = x;

      if (x < meps || x > eps)
      {
         *ii++ = int(v - val);
         *v++ = x;
      }
      num = int(ii - idx);
   }
   setupStatus = 1;

   assert(isConsistent());
   return *this;
}

void SSVector::setup_and_assign(SSVector& rhs)
{
   clear();

   DIdxSet::setMax(rhs.max());
   DVector::reDim(rhs.dim());

   if (rhs.isSetup())
   {
      int i, j;
      DIdxSet::operator=(rhs);
      for (i = size() - 1; i >= 0; --i)
      {
         j = index(i);
         val[j] = rhs.val[j];
      }
   }
   else
   {
      int* ri = rhs.idx;
      int* ii = idx;
      double* rv = rhs.val;
      double* v = val;
      double* last = rv + rhs.dim() - 1;
      double x = *last;
      const double eps = rhs.epsilon;
      const double meps = -eps;

      *last = 1e-100;
      for(;;)
      {
         while (!*rv)
         {
            ++rv;
            ++v;
         }
         if (*rv < meps || *rv > eps)
         {
            *ri++ = *ii++ = int(v - val);
            *v++ = *rv++;
         }
         else if (rv == last)
            break;
         else
         {
            v++;
            *rv++ = 0;
         }
      }

      if (x < meps || x > eps)
      {
         *ri++ = *ii++ = int(v - val);
         *v++ = *rv = x;
      }
      else
         *rv = 0;
      num = rhs.num = int(ii - idx);
      rhs.setupStatus = 1;
   }
   setupStatus = 1;

   assert(isConsistent());
}

SSVector& SSVector::operator=(const SVector& rhs)
{
   clear();
   return assign(rhs);
}

SSVector& SSVector::assign(const SVector& rhs)
{
   assert(rhs.dim() <= Vector::dim());

   const SVector::Element* e = rhs.m_elem;
   int* p = idx;
   int i = rhs.size();

   while (i--)
   {
      val[*p = e->idx] = e->val;
      p += ((e++)->val != 0);
   }
   num = int(p - idx);
   setupStatus = 1;

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::assign2product1(const SVSet& A, const SSVector& x)
{
   assert(x.isSetup());

   const double* vl = x.val;
   const int* xi = x.idx;

   int* ii = idx;
   SVector* svec = const_cast<SVector*>( & A[*xi] );
   const SVector::Element* e = &svec->element(0);
   const SVector::Element* last = e + (num = svec->size());
   double* v = val;
   double y = vl[*xi];

   for (; e < last; ++e)
      v[ *ii++ = e->idx ] = y * e->val;

   return *this;
}

SSVector& SSVector::assign2productShort(const SVSet& A, const SSVector& x)
{
   assert(x.isSetup());

   int i, j;
   const double* vl = x.val;
   const int* xi = x.idx;

   double y;
   int* ii = idx;
   SVector* svec = const_cast<SVector*>( & A[*xi] );
   const SVector::Element* e = &svec->element(0);
   const SVector::Element* last = e + (num = svec->size());
   double* v = val;
   double xx = vl[*xi++];
   for (; e < last; ++e)
   {
      v[ *ii = e->idx ] = y = xx * e->val;
      ii += (y != 0);
   }

   int k;
   double mark = 1e-100;
   for (i = x.size(); --i > 0;)
   {
      xx = vl[*xi];
      svec = const_cast<SVector*>( & A[*xi++] );
      e = &svec->element(0);
      for (k = svec->size(); --k >= 0;)
      {
         *ii = j = e->idx;
         ii += ((y = v[j]) == 0);
         y += xx * e++->val;
         v[j] = y + (y == 0) * mark;
      }
   }

   const double eps = epsilon;
   const double meps = -eps;
   int* is = idx;
   int* it = idx;
   for (; is < ii; ++is)
   {
      y = v[*is];
      if (y > eps || y < meps)
         *it++ = *is;
      else
         v[*is] = 0;
   }
   num = int(it - idx);

   assert(isConsistent());
   return *this;
}

SSVector& SSVector::assign2productFull(const SVSet& A, const SSVector& x)
{
   assert(x.isSetup());

   int i;
   const double* vl = x.val;
   const int* xi = x.idx;

   SVector* svec;
   const SVector::Element* elem;
   const SVector::Element* last;
   double y;
   double* v = val;

   for (i = x.size(); i-- > 0; ++xi)
   {
      svec = const_cast<SVector*>( & A[*xi] );
      elem = &svec->element(0);
      last = elem + svec->size();
      y = vl[*xi];
      for (; elem < last; ++elem)
         v[elem->idx] += y * elem->val;
   }

   return *this;
}

SSVector& SSVector::assign2product4setup(const SVSet& A, const SSVector& x)
{
   assert(A.num() == x.dim());

   assert(x.isSetup());

   clear();

   if (x.size() == 1)
   {
      assign2product1(A, x);
      setupStatus = 1;
   }

   else if (double(x.size())*A.memSize() <= shortProductFactor*dim()*A.num()
             && isSetup())
   {
      assign2productShort(A, x);
      setupStatus = 1;
   }

   else
   {
      assign2productFull(A, x);
      setupStatus = 0;
   }

   return *this;
}

SSVector& SSVector::assign2product(const SSVector& x, const SVSet& A)
{
   assert(A.num() == dim());
   const double eps = epsilon;
   const double minuseps = -epsilon;
   double y;

   clear();
   for (int i = dim(); i-- > 0;)
   {
      y = A[i] * x;
      if (y > eps || y < minuseps)
      {
         val[i] = y;
         IdxSet::addIdx(i);
      }
   }

   return *this;
}

SSVector& SSVector::assign2productAndSetup(const SVSet& A, SSVector& x)
{
   if (x.isSetup())
      return assign2product4setup(A, x);

   SVector* svec;
   const SVector::Element* elem;
   const SVector::Element* last;
   double y;
   double* v = val;
   int* xi = x.idx;
   double* xv = x.val;
   double* end = xv + x.dim() - 1;
   const double eps = epsilon;
   const double meps = -epsilon;

   /* setze weissen Elefanten */
   double lastval = *end;
   *end = 1e-100;

   for(;;)
   {
      while (!*xv)
         ++xv;
      if (*xv > eps || *xv < meps)
      {
         y = *xv;
         svec = const_cast<SVector*>( & A[ *xi++ = int(xv - x.val) ] );
         elem = &svec->element(0);
         last = elem + svec->size();
         for (; elem < last; ++elem)
            v[elem->idx] += y * elem->val;
      }
      else
      {
         *xv = 0;
         if (xv == end)
            break;
      }
      xv++;
   }

   /* fange weissen Elefanten wieder ein */
   if (lastval > eps || lastval < meps)
   {
      y = *xv = lastval;
      svec = const_cast<SVector*>( & A[ *xi++ = int(xv - x.val) ] );
      elem = &svec->element(0);
      last = elem + svec->size();
      for (; elem < last; ++elem)
         v[elem->idx] += y * elem->val;
   }
   else
      *xv = 0;

   x.num = int(xi - x.idx);
   x.setupStatus = 1;
   setupStatus = 0;

   return *this;
}

int SSVector::isConsistent() const
{
   if (Vector::dim() > DIdxSet::max())
      return SPXinconsistent("SSVector");
   if (Vector::dim() < DIdxSet::dim())
      return SPXinconsistent("SSVector");

   if (isSetup())
   {
      for (int i = Vector::dim() - 1; i >= 0; --i)
         if (val[i] != 0 && number(i) < 0)
            return SPXinconsistent("SSVector");
   }

   return DVector::isConsistent() && DIdxSet::isConsistent();
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
