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
#pragma ident "@(#) $Id: lprowset.cpp,v 1.2 2001/11/06 23:31:02 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "lprowset.h"

namespace soplex
{


//@ ----------------------------------------------------------------------------
/*      \SubSection{Extension}
 */
void LPRowSet::add(Key& key, double lhs, const SVector& vector, double rhs)
{
   SVSet::add(key, vector);
   if (num() > left.dim())
   {
      left.reDim(num());
      right.reDim(num());
   }
   left [num() - 1] = lhs;
   right[num() - 1] = rhs;
}

void LPRowSet::add(const LPRowSet& set)
{
   int i = num();

   SVSet::add(set);
   if (num() > left.dim())
   {
      left.reDim(num());
      right.reDim(num());
   }

   for (int j = 0; i < num(); ++i, ++j)
   {
      left [i] = set.lhs(j);
      right[i] = set.rhs(j);
   }
}

void LPRowSet::add(Key nkey[], const LPRowSet& set)
{
   int i = num();
   add(set);

   for (int j = 0; i < num(); ++i, ++j)
      nkey[j] = key(i);
}

SVector& LPRowSet::create(Key& nkey, int nonzeros, double lhs, double rhs)
{
   if (num() + 1 > left.dim())
   {
      left.reDim(num() + 1);
      right.reDim(num() + 1);
   }
   left [num()] = lhs;
   right[num()] = rhs;

   return *SVSet::create(nkey, nonzeros);
}


/*      \SubSection{Shrinking}
 */
void LPRowSet::remove(int i)
{
   SVSet::remove(i);
   left[i] = left[num()];
   right[i] = right[num()];
   left.reDim (num());
   right.reDim(num());
}

void LPRowSet::remove(int perm[])
{
   int i;
   int j = num();
   SVSet::remove(perm);
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0 && perm[i] != i)
      {
         left[perm[i]] = left[i];
         right[perm[i]] = right[i];
      }
   }
   left.reDim (num());
   right.reDim(num());
}

void LPRowSet::remove(int nums[], int n, int* perm)
{
   SVSet::remove(nums, n, perm);
   int i;
   int j = num();
   for (i = 0; i < j; ++i)
   {
      if (perm[i] >= 0 && perm[i] != i)
      {
         left[perm[i]] = left[i];
         right[perm[i]] = right[i];
      }
   }
   left.reDim (num());
   right.reDim(num());
}

void LPRowSet::clear()
{
   SVSet::clear();
   left.reDim (num());
   right.reDim(num());
}

// TK13OCT1998
void LPRowSet::setType(
   int i,
   LPRow::Type type)
{
   switch (type)
   {
   case LPRow::LESS_EQUAL:
      lhs(i) = -LPRow::infinity;
      break;
   case LPRow::EQUAL:
      if (lhs(i) > -LPRow::infinity)
         rhs(i) = lhs(i);
      else
         lhs(i) = rhs(i);
      break;
   case LPRow::GREATER_EQUAL:
      rhs(i) = LPRow::infinity;
      break;
   case LPRow::RANGE :
      std::cerr << __FILE__ << __LINE__
      << "RANGE not supported in LPRowSet::setType()";
      /*FALLTHROUGH*/
   default:
      assert(0);
   }
}


/*      \SubSection{Consistency}
 */
#define inconsistent                                                    \
{                                                                       \
std::cout << "ERROR: Inconsistency detected in class LPRowSet\n";        \
return 0;                                                          \
}

int LPRowSet::isConsistent() const
{
   if (left.dim() != right.dim())
      inconsistent;
   if (left.dim() != num())
      inconsistent;

   return left.isConsistent() && right.isConsistent() && SVSet::isConsistent();
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
