/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "soplex/idxset.h"

namespace soplex
{

int IdxSet::dim() const
{
   int ddim = -1;

   for (int i = 0; i < size(); i++)
      if (ddim < idx[i])
         ddim = idx[i];

   return ddim;
}

int IdxSet::pos(int i) const
{
   for(int n = 0; n < size(); n++)
      if (idx[n] == i)
         return n;

   return -1;
}

  // Add n elements from an std::vector
  void IdxSet::add(int n, const std::vector<int> i)
{
   assert(n >= 0 && size() + n <= max());
   idx.insert(idx.end(), i.begin(), i.begin() + n);
   add(n);
}

void IdxSet::remove(int n, int m)
{
   assert(n <= m && m < size() && n >= 0);

   // erase(i, j) does an erase of the range [i, j)
   idx.erase(idx.begin() + n, idx.begin() + m + 1);
}

IdxSet& IdxSet::operator=(const IdxSet& rhs)
{
   if (this != &rhs)
   {
     if(idx.empty())
      {
         freeArray = true;
      }

     idx = rhs.idx;
   }

   return *this;
}

bool IdxSet::isConsistent() const
{
   return true;
}
} // namespace soplex
