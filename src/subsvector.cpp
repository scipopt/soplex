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
#pragma ident "@(#) $Id: subsvector.cpp,v 1.3 2001/12/25 16:03:25 bzfkocht Exp $"

#include <assert.h>
#include <iostream>

#include "subsvector.h"
#include "spxmessage.h"

namespace soplex
{

void SubSVector::sort()
{
   SVector::Element dummy;
   SVector::Element* w;
   SVector::Element* l;
   SVector::Element* e = elem + size();
   for (l = elem, w = elem + 1; w < e; l = w, ++w)
   {
      if (l->idx > w->idx)
      {
         dummy = *w;
         do
         {
            l[1] = *l;
            if (l-- == elem)
               break;
         }
         while (l->idx > dummy.idx);
         l[1] = dummy;
      }
   }
}

int SubSVector::dim() const
{
   SVector::Element* e = elem;
   int d = -1;
   int n = size();
   while (n--)
   {
      d = (d > e->idx) ? d : e->idx;
      e++;
   }
   return d;
}

int SubSVector::number(int i) const
{
   int n = size();
   SVector::Element* e = &(elem[n]);
   while (n--)
   {
      --e;
      if (e->idx == i)
         return n;
   }
   return -1;
}

double SubSVector::length2() const
{
   double x = 0;
   int n = size();
   SVector::Element* e = elem;
   while (n--)
   {
      x += e->val * e->val;
      e++;
   }
   return x;
}

SubSVector& SubSVector::operator*=(double x)
{
   int n = size();
   SVector::Element* e = elem;
   while (n--)
   {
      e->val *= x;
      e++;
   }
   return *this;
}

double SubSVector::operator*(const Vector& w) const
{
   double x = 0;
   int n = size();
   SVector::Element* e = elem;
   while (n--)
   {
      x += e->val * w[e->idx];
      e++;
   }
   return x;
}

std::ostream& operator<<(std::ostream& os, const SubSVector& v)
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

int SubSVector::isConsistent() const
{
   if (elem)
   {
#ifndef NDEBUG
      if (elem < &svec->element(0))
         return SPXinconsistent("SubSVector");
      if (elem + num > (&svec->element(0)) + svec->size())
         return SPXinconsistent("SubSVector");
      return svec->isConsistent();
#endif
   }
   else if (num)
      return SPXinconsistent("SubSVector");
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
