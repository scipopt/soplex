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
#pragma ident "@(#) $Id: lpcolset.cpp,v 1.1 2001/11/06 16:18:31 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "lpcolset.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
/*      \SubSection{Extension}
 */
void LPColSet::add
(
   Key& key,
   double obj,
   double lower,
   const SVector& colVector,
   double upper
)
{
   SVSet::add(key, colVector);
   if (num() > low.dim())
   {
      low.reDim (num());
      up.reDim (num());
      object.reDim (num());
   }
   low [num() - 1] = lower;
   up [num() - 1] = upper;
   object[num() - 1] = obj;
}

void LPColSet::add(const LPColSet& set)
{
   int i = num();

   SVSet::add(set);
   if (num() > low.dim())
   {
      low.reDim (num());
      up.reDim (num());
      object.reDim (num());
   }

   for (int j = 0; i < num(); ++i, ++j)
   {
      low [i] = set.lower(j);
      up [i] = set.upper(j);
      object[i] = set.obj(j);
   }
}

void LPColSet::add(Key nkey[], const LPColSet& set)
{
   int i = num();
   add(set);

   for (int j = 0; i < num(); ++i, ++j)
      nkey[j] = key(i);
}

SVector& LPColSet::create(Key& nkey, int nonzeros, double obj, double lhs, double rhs)
{
   if (num() + 1 > low.dim())
   {
      low.reDim (num() + 1);
      up.reDim (num() + 1);
      object.reDim (num() + 1);
   }
   low [num()] = lhs;
   up [num()] = rhs;
   object[num()] = obj;

   return *SVSet::create(nkey, nonzeros);
}


/*      \SubSection{Shrinking}
 */
void LPColSet::remove(int i)
{
   SVSet::remove(i);
   low[i] = low[num()];
   up[i] = up[num()];
   object[i] = object[num()];
   low.reDim (num());
   up.reDim (num());
   object.reDim(num());
}

void LPColSet::remove(int perm[])
{
   int i;
   int j = num();
   SVSet::remove(perm);
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0 && perm[i] != i)
      {
         low[perm[i]] = low[i];
         up[perm[i]] = up[i];
         object[perm[i]] = object[i];
      }
   }
   low.reDim (num());
   up.reDim (num());
   object.reDim(num());
}

void LPColSet::remove(int nums[], int n, int* perm)
{
   SVSet::remove(nums, n, perm);
   int i;
   int j = num();
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0 && perm[i] != i)
      {
         low[perm[i]] = low[i];
         up[perm[i]] = up[i];
         object[perm[i]] = object[i];
      }
   }
   low.reDim (num());
   up.reDim (num());
   object.reDim(num());
}

void LPColSet::clear()
{
   SVSet::clear();
   low.reDim (num());
   up.reDim (num());
   object.reDim(num());
}


/*      \SubSection{Consistency}
 */
#define inconsistent                                                    \
{                                                                       \
std::cout << "ERROR: Inconsistency detected in class LPColSet\n";        \
return 0;                                                          \
}

int LPColSet::isConsistent() const
{
   if (low.dim() != object.dim())
      inconsistent;
   if (low.dim() != up.dim())
      inconsistent;
   if (low.dim() != num())
      inconsistent;

   return low.isConsistent() && up.isConsistent() && SVSet::isConsistent();
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
