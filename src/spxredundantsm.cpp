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
#pragma ident "@(#) $Id: spxredundantsm.cpp,v 1.10 2002/01/29 15:38:49 bzfkocht Exp $"

#include <iostream>

#include "real.h"
#include "spxredundantsm.h"
#include "dataarray.h"

namespace soplex
{
static const Real eps = 1e-10;    ///< epslion for what is regarded equal.

int SPxRedundantSM::simplify()
{
   int    num;
   int    j;
   int    i;
   int    k;
   int    upcnt;
   int    locnt;
   Real up;
   Real lo;
   Real x;
   Real y;
   DataArray < int > rem(lp->nRows());

   num = 0;
   rem.reSize(lp->nCols());
   for (i = lp->nCols() - 1; i >= 0; --i)
   {
      const SVector& col = lp->colVector(i);
      rem[i] = 0;
      if (lp->upper(i) != lp->lower(i))
      {
         up = 0;
         lo = 0;
         for (j = col.size() - 1; j >= 0 && (up == 0 || lo == 0); --j)
         {
            x = col.value(j);
            k = col.index(j);
            if (x > 0)
            {
               up += (lp->rhs(k) < infinity);
               lo += (lp->lhs(k) > -infinity);
            }
            else if (x < 0)
            {
               lo += (lp->rhs(k) < infinity);
               up += (lp->lhs(k) > -infinity);
            }
         }
         x = lp->maxObj(i);
         if (lo == 0 && x < 0)
         {
            if (lp->lower(i) <= -infinity)
               return 1;           // LP is unbounded
            lp->changeUpper(i, lp->lower(i));
         }
         else if (up == 0 && x > 0)
         {
            if (lp->upper(i) >= infinity)
               return 1;           // LP is unbounded
            lp->changeLower(i, lp->upper(i));
         }
         else if (x == 0)
         {
            up += (lp->upper(i) < infinity);
            lo += (lp->lower(i) > -infinity);
            if (lo == 0)
            {
               lp->changeUpper(i, infinity);
               for (j = col.size() - 1; j >= 0; --j)
               {
                  x = col.value(j);
                  k = col.index(j);
                  if (x > 0)
                     lp->changeRhs(k, infinity);
                  else
                     lp->changeLhs(k, -infinity);
               }
            }
            if (up == 0)
            {
               lp->changeLower(i, -infinity);
               for (j = col.size() - 1; j >= 0; --j)
               {
                  x = col.value(j);
                  k = col.index(j);
                  if (x < 0)
                     lp->changeRhs(k, infinity);
                  else
                     lp->changeLhs(k, -infinity);
               }
            }
         }
      }
      if ((x = lp->upper(i)) == lp->lower(i))
      {
         rem[i] = -1;
         num++;
         if (x != 0)
         {
            for (j = col.size() - 1; j >= 0; --j)
            {
               k = col.index(j);
               if (lp->rhs(k) < infinity)
                  lp->changeRhs(k, lp->rhs(k) - x*col.value(j));
               if (lp->lhs(k) > -infinity)
                  lp->changeLhs(k, lp->lhs(k) - x*col.value(j));
            }
            delta += x * lp->obj(i);
         }
      }
   }
   if (num)
   {
      lp->removeCols(rem.get_ptr());
      std::cout << "SPxRedundantSM: removed " << num << " column(s)\n";
      assert(lp->isConsistent());
   }

   num = 0;
   for (i = lp->nRows() - 1; i >= 0; --i)
   {
      if (lp->rhs(i) < infinity || lp->lhs(i) > -infinity)
      {
         const SVector& row = lp->rowVector(i);

         rem[i] = 0;

         up = lo = 0;
         upcnt = locnt = 0;
         for (j = row.size() - 1; j >= 0; --j)
         {
            x = row.value(j);
            k = row.index(j);
            if (x > 0)
            {
               if (lp->upper(k) >= infinity)
                  upcnt++;
               else
                  up += lp->upper(k) * x;
               if (lp->lower(k) <= -infinity)
                  locnt++;
               else
                  lo += lp->lower(k) * x;
            }
            else if (x < 0)
            {
               if (lp->upper(k) >= infinity)
                  locnt++;
               else
                  lo += lp->upper(k) * x;
               if (lp->lower(k) <= -infinity)
                  upcnt++;
               else
                  up += lp->lower(k) * x;
            }
         }

         if (((lp->rhs(i) >= up - eps && upcnt <= 0)
               || lp->rhs(i) >= infinity)
              && ((lp->lhs(i) <= lo + eps && locnt <= 0)
                  || lp->lhs(i) <= -infinity))
         {
            rem[i] = -1;
            num++;
         }
         else if ((lp->rhs(i) < lo - eps && locnt <= 0)
                   || (lp->lhs(i) > up + eps && upcnt <= 0))
            return -1;
         else
         {
            /*
                if(lp->lhs(i) <= lo+eps && locnt <= 0)
                    lp->changeLhs(i, -infinity);
                else if(lp->rhs(i) >= up-eps && upcnt <= 0)
                    lp->changeRhs(i, infinity);
                else
             */
            if (upcnt < 2 || locnt < 2)
            {
               for (j = row.size() - 1; j >= 0; --j)
               {
                  x = row.value(j);
                  k = row.index(j);
                  if (x > 0)
                  {
                     if (lp->lhs(i) > -infinity
                          && lp->lower(k) > -infinity
                          && upcnt < 2)
                     {
                        y = -infinity;
                        if (lp->upper(k) < infinity && upcnt < 1)
                           y = lp->upper(k) + (lp->lhs(i) - up) / x;
                        else if (lp->upper(k) >= infinity)
                           y = lp->lhs(i) - up;
                        if (y >= lp->lower(k))
                        {
                           lp->changeLower(k, -infinity);
                           break;
                        }
                     }
                     if (lp->rhs(i) < infinity
                          && lp->upper(k) < infinity
                          && locnt < 2)
                     {
                        y = infinity;
                        if (lp->lower(k) > -infinity && locnt < 1)
                           y = lp->lower(k) + (lp->rhs(i) - lo) / x;
                        else if (lp->lower(k) <= -infinity)
                           y = lp->rhs(i) - lo;
                        if (y <= lp->upper(k))
                        {
                           lp->changeUpper(k, infinity);
                           break;
                        }
                     }
                  }
                  else if (x < 0)
                  {
                     if (lp->lhs(i) >= -infinity
                          && lp->upper(k) < infinity
                          && upcnt < 2)
                     {
                        y = infinity;
                        if (lp->lower(k) > -infinity && upcnt < 1)
                           y = lp->lower(k) + (lp->lhs(i) - up) / x;
                        else if (lp->lower(k) <= -infinity)
                           y = -(lp->lhs(i) - up);
                        if (y <= lp->upper(k))
                        {
                           lp->changeUpper(k, infinity);
                           break;
                        }
                     }
                     if (lp->rhs(i) <= infinity
                          && lp->lower(k) > -infinity
                          && locnt < 2)
                     {
                        y = -infinity;
                        if (lp->upper(k) < infinity && locnt < 1)
                           y = lp->upper(k) + (lp->rhs(i) - lo) / x;
                        else if (lp->upper(k) >= infinity)
                           y = -(lp->rhs(i) - lo);
                        if (y >= lp->lower(k))
                        {
                           lp->changeLower(k, -infinity);
                           break;
                        }
                     }
                  }
               }
            }
         }
      }
      else
         rem[i] = -1;
   }
   if (num)
   {
      lp->removeRows(rem.get_ptr());
      std::cout << "SPxRedundantSM:\tremoved " << num << " row(s)\n";
      assert(lp->isConsistent());
   }
   return 0;
}

void SPxRedundantSM::unsimplify()
{
   std::cerr << "SPxRedundantSM::unsimplify() not implemented\n";
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
