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
#pragma ident "@(#) $Id: svector.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "svector.h"


#include "ssvector.h"

namespace soplex
{






//@ ----------------------------------------------------------------------------
/*      \SubSection{Modification}
 */
void SVector::toFront(int n)
{
   assert(n < size() && n >= 0);
   Element dummy = elem[n];
   elem[n] = elem[0];
   elem[0] = dummy;
}

void SVector::add(int n, const int i[], const double v[])
{
   assert(n + size() <= max());
   Element* e = elem + size();
   size() += n;
   while (n--)
   {
      e->idx = *i++;
      e->val = *v++;
      e++;
   }
}

void SVector::add(int n, const Element e[])
{
   assert(n + size() <= max());
   Element* ee = elem + size();
   size() += n;
   while (n--)
      *ee++ = *e++;
}

void SVector::remove(int n, int m)
{
   assert(n <= m && m < size() && n >= 0);
   ++m;

   int cpy = m - n;
   cpy = (size() - m >= cpy) ? cpy : size() - m;

   Element* e = &elem[size() - 1];
   Element* r = &elem[n];
   size() -= cpy;
   do
   {
      *r++ = *e--;
   }
   while (--cpy);
}

int SVector::dim() const
{
   Element* e = elem;
   int d = -1;
   int n = size();
   while (n--)
   {
      d = (d > e->idx) ? d : e->idx;
      e++;
   }
   return d;
}

void SVector::sort()
{
   Element dummy;
   Element* w;
   Element* l;
   Element* s = &(elem[0]);
   Element* e = s + size();
   for (l = s, w = s + 1; w < e; l = w, ++w)
   {
      if (l->idx > w->idx)
      {
         dummy = *w;
         do
         {
            l[1] = *l;
            if (l-- == s)
               break;
         }
         while (l->idx > dummy.idx);
         l[1] = dummy;
      }
   }
}


//@ ----------------------------------------------------------------------------
/*      \SubSection{Maths}
 */
double SVector::length2() const
{
   double x = 0;
   int n = size();
   Element* e = elem;
   while (n--)
   {
      x += e->val * e->val;
      e++;
   }
   return x;
}

double SVector::maxAbs() const
{
   double x = 0;
   int n = size();
   Element* e = elem;
   while (n--)
   {
      x = (e->val > x) ? e->val : ((-e->val > x) ? -e->val : x);
      e++;
   }
   return x;
}

SVector& SVector::operator*=(double x)
{
   int n = size();
   Element* e = elem;
   while (n--)
   {
      e->val *= x;
      e++;
   }
   return *this;
}

//@ ----------------------------------------------------------------------------
/*      \SubSection{Miscellaneous}
 */
SVector& SVector::operator=(const SSVector& sv)
{
   assert(max() >= sv.size());
   size() = sv.size();
   int i = size();
   Element *e = elem;

   while (i--)
   {
      e->idx = sv.index(i);
      e->val = sv[e->idx];
      ++e;
   }

   return *this;
}

SVector& SVector::operator=(const Vector& vec)
{
   int n = 0;
   int i = vec.dim();
   Element *e = elem;
   clear();
   while (i--)
   {
      if (vec[i])
      {
         assert(n < max());
         e->idx = i;
         e->val = vec[i];
         ++e;
         ++n;
      }
   }
   size() = n;
   assert(isConsistent());
   return *this;
}

SVector& SVector::assign(const Vector& vec, double eps)
{
   int n = 0;
   int i = vec.dim();
   double x;
   Element* e = elem;
   clear();
   while (i--)
   {
      x = vec[i];
      if (x > eps || x < -eps)
      {
         assert(n < max());
         e->idx = i;
         e->val = x;
         ++e;
         ++n;
      }
   }
   size() = n;
   assert(isConsistent());
   return *this;
}

SVector& SVector::operator=(const SVector& sv)
{
   assert(max() >= sv.size());
   int i = size() = sv.size();
   Element *e = elem;
   Element *s = sv.elem;
   while (i--)
      *e++ = *s++;
   return *this;
}

std::ostream& operator<<(std::ostream& os, const SVector& v)
{
   int i, j;
   for (i = j = 0; i < v.size(); ++i)
   {
      if (j)
      {
         if (v.value(i) < 0)
            os << " - " << -v.value(i);
         else
            os << " + " << v.value(i);
      }
      else
         os << v.value(i);
      os << " x" << v.index(i);
      j = 1;
      if ((i + 1) % 4 == 0)
         os << "\n\t";
   }
   return os;
}


//@ ----------------------------------------------------------------------------
/*      \SubSection{Consistency}
 */
#define inconsistent                                                    \
{                                                                       \
std::cout << "ERROR: Inconsistency detected in class SVector\n"; \
return 0;                                                          \
}

int SVector::isConsistent() const
{
   if (elem)
   {
      if (size() > max())
         inconsistent;
      for (int i = 1; i < size(); ++i)
      {
         for (int j = 0; j < i; ++j)
         {
            if (elem[i].idx == elem[j].idx)
               inconsistent;
         }
      }
   }
   return 1;
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
