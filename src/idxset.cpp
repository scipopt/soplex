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
#pragma ident "@(#) $Id: idxset.cpp,v 1.2 2001/11/06 23:31:01 bzfkocht Exp $"


/*      \Section{Complex Members}
 */
#include <iostream>
#include "idxset.h"

namespace soplex
{

/*      \SubSection{Inquiry}
 */
int IdxSet::dim() const
{
   int dim = -1;
   for (int i = 0; i < size(); ++i)
      dim = (idx[i] > dim) ? idx[i] : dim;
   return dim -1;
}

int IdxSet::number(int i) const
{
   for (int n = size() - 1; n >= 0; --n)
   {
      if (idx[n] == i)
         return n;
   }
   return -1;
}

/*      \SubSection{Manipulation}
 */
void IdxSet::add(int n, const int i[])
{
   assert(n >= 0 && size() + n <= max());
   for (int j = 0; j < n; ++j)
      idx[size() + j] = i[j];
   add(n);
}

void IdxSet::toFront(int n)
{
   assert(n >= 0 && n < size());
   int idx = index(0);
   index(0) = index(n);
   index(n) = idx;
}

void IdxSet::remove(int n, int m)
{
   assert(n <= m && m < size() && n >= 0);
   ++m;

   int cpy = m - n;
   int newnum = num - cpy;
   cpy = (size() - m >= cpy) ? cpy : size() - m;

   do
   {
      --num;
      --cpy;
      idx[n + cpy] = idx[num];
   }
   while (cpy > 0);
   num = newnum;
}

/*      \SubSection{Miscellaneous}
 */
IdxSet& IdxSet::operator=(const IdxSet& set)
{
   assert(max() >= set.size());
   for (num = 0; num < set.size(); ++num)
      idx[num] = set.idx[num];
   return *this;
}

#define inconsistent                                            \
{                                                               \
std::cerr << "Inconsistency detected in class IdxSet\n";        \
return 0;                                                  \
}

int IdxSet::isConsistent() const
{
   int i, j;

   if (len > 0 && idx == 0)
      inconsistent;

   for (i = 0; i < size(); ++i)
   {
      if (index(i) < 0)
         inconsistent;
      for (j = 0; j < i; ++j)
         if (index(i) == index(j))
            inconsistent;
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
