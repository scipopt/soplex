/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1997-1999 Roland Wunderling                              */
/*                  1997-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: subsvector.cpp,v 1.8 2002/03/03 13:50:35 bzfkocht Exp $"

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "subsvector.h"
#include "message.h"

namespace soplex
{
int SubSVector::dim() const
{
   const SVector::Element* e = elem;
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
   const SVector::Element* e = &(elem[n]);
   while (n--)
   {
      --e;
      if (e->idx == i)
         return n;
   }
   return -1;
}

Real SubSVector::length2() const
{
   Real x = 0;
   int n = size();
   const SVector::Element* e = elem;
   while (n--)
   {
      x += e->val * e->val;
      e++;
   }
   return x;
}

Real SubSVector::operator*(const Vector& w) const
{
   Real x = 0;
   int n = size();
   const SVector::Element* e = elem;
   while (n--)
   {
      x += e->val * w[e->idx];
      e++;
   }
   return x;
}

bool SubSVector::isConsistent() const
{
   if (elem)
   {
#ifndef NDEBUG
      if (elem < &svec->element(0))
         return MSGinconsistent("SubSVector");
      if (elem + num > (&svec->element(0)) + svec->size())
         return MSGinconsistent("SubSVector");
      return svec->isConsistent();
#endif
   }
   else if (num != 0)
      return MSGinconsistent("SubSVector");
   return true;
}

std::ostream& operator<<(std::ostream& os, const SubSVector& v)
{
   bool first = true;

   for (int i = 0; i < v.size(); ++i)
   {
      if (first)
         os << v.value(i);
      else
      {
         if (v.value(i) < 0)
            os << " - " << -v.value(i);
         else
            os << " + " << v.value(i);
      }
      os << " x" << v.index(i);
      first = false;
      if ((i + 1) % 4 == 0)
         os << "\n\t";
   }
   return os;
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
