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
#pragma ident "@(#) $Id: spxredundantsm.cpp,v 1.20 2003/01/05 19:03:17 bzfkocht Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxredundantsm.h"
#include "dataarray.h"

namespace soplex
{

SPxSimplifier::Result SPxRedundantSM::treat_cols(SPxLP& lp)
{
   DataArray<int> rem(lp.nCols());
   DataArray<int> tmp(lp.nCols());
   int            num = 0;
   int            j;
   int            k;
   Real           x;

   for( int i = 0; i < lp.nCols(); ++i )
   {
      rem[i] = 0;

      const SVector& col = lp.colVector(i);

      if (lp.upper(i) != lp.lower(i))
      {
         // test if all coefficents are going in one direction
         int upcnt = 0;
         int locnt = 0;
         
         for( j = 0; j < col.size(); ++j )
         {            
            if (upcnt > 0 && locnt > 0)
               break;

            x = col.value(j);
            k = col.index(j);

            assert(isNotZero(x));

            if (x > 0.0)
            {
               upcnt += (lp.rhs(k) <  infinity) ? 1 : 0;
               locnt += (lp.lhs(k) > -infinity) ? 1 : 0;
            }
            else if (x < 0.0)
            {
               locnt += (lp.rhs(k) <  infinity) ? 1 : 0;
               upcnt += (lp.lhs(k) > -infinity) ? 1 : 0;
            }
         }
         x = lp.maxObj(i);

         // max -3 x
         // s.t. 5 x <= 8
         if (locnt == 0 && x < 0.0)
         {
            if (lp.lower(i) <= -infinity)
               return UNBOUNDED;           // LP is unbounded

            lp.changeUpper(i, lp.lower(i));
         }
         // max  3 x
         // s.t. 5 x >= 8
         else if (upcnt == 0 && x > 0.0)
         {
            if (lp.upper(i) >= infinity)
               return UNBOUNDED;           // LP is unbounded

            lp.changeLower(i, lp.upper(i));
         }
         else if (x == 0.0)
         {
            upcnt += (lp.upper(i) <  infinity) ? 1 : 0;
            locnt += (lp.lower(i) > -infinity) ? 1 : 0;

            if (locnt == 0)
            {
               // make variable free
               lp.changeUpper(i, infinity);

               for( j = 0; j < col.size(); ++j )
               {
                  if (col.value(j) < 0.0)
                     assert(lp.rhs(col.index(j)) >= infinity);
                  //                     lp.changeRhs(col.index(j), infinity);
                  else
                     assert(lp.lhs(col.index(j)) <= -infinity);
                  //                     lp.changeLhs(col.index(j), -infinity);
               }
            }
            if (upcnt == 0)
            {
               // make variable free
               lp.changeLower(i, -infinity);

               for( j = 0; j < col.size(); ++j )
               {
                  if (col.value(j) > 0.0)
                     assert(lp.rhs(col.index(j)) >= infinity);
                  //                     lp.changeRhs(col.index(j), infinity);
                  else
                     assert(lp.lhs(col.index(j)) <= -infinity);
                  //                     lp.changeLhs(col.index(j), -infinity);
               }
            }
         }
      }
      // remove fixed variables
      if (lp.upper(i) == lp.lower(i))
      {
         x      = lp.upper(i);
         rem[i] = -1;
         num++;

         if (x != 0.0)
         {
            for( j = 0; j < col.size(); ++j )
            {
               k = col.index(j);

               if (lp.rhs(k) <  infinity)
                  lp.changeRhs(k, lp.rhs(k) - x * col.value(j));
               if (lp.lhs(k) > -infinity)
                  lp.changeLhs(k, lp.lhs(k) - x * col.value(j));
            }
            // Zielfunktiondelta delta += x * lp.obj(i);
         }
      }
   }
   if (num > 0)
   {
      lp.removeCols(rem.get_ptr());
      assert(lp.isConsistent());

      VERBOSE1({ std::cout << "SPxRedundantSM: removed " << num
                           << " column(s)" << std::endl; });
   }
   return OKAY;
}


SPxSimplifier::Result SPxRedundantSM::treat_rows(SPxLP& lp)
{
   DataArray < int > rem(lp.nRows());
   int               num = 0;
   int               j;
   int               k;
   Real              x;

   for( int i = 0; i < lp.nRows(); ++i )
   {
      // unconstraint constraints and infeasible constraints are
      // also handled by rem1
      if (lp.rhs(i) >= infinity && lp.lhs(i) <= -infinity)
      {
         rem[i] = -1;
         num++;
         continue;
      }
      if (LT(lp.rhs(i), lp.lhs(i)))
      {
         rem[i] = -1;
         num++;
         continue;
      }
      rem[i] = 0;

      const SVector& row   = lp.rowVector(i);       
      Real           upbnd  = 0.0;
      Real           lobnd = 0.0;
      int            upcnt = 0;
      int            locnt = 0;

      for( j = 0; j < row.size(); ++j )
      {
         x = row.value(j);
         k = row.index(j);
         
         assert(isNotZero(x));

         if (x > 0)
         {
            if (lp.upper(k) >= infinity)
               upcnt++;
            else
               upbnd += lp.upper(k) * x;
            
            if (lp.lower(k) <= -infinity)
               locnt++;
            else
               lobnd += lp.lower(k) * x;
         }
         else if (x < 0)
         {
            if (lp.upper(k) >= infinity)
               locnt++;
            else
               lobnd += lp.upper(k) * x;
            
            if (lp.lower(k) <= -infinity)
               upcnt++;
            else
               upbnd += lp.lower(k) * x;
         }
      }
      // infeasible ?
      if ((LT(lp.rhs(i), lobnd) && locnt == 0) || (GT(lp.lhs(i), upbnd) && upcnt == 0))
         return INFEASIBLE;

      // redundant rhs ?
      if (lp.rhs(i) <=  infinity && upcnt == 0 && GE(lp.rhs(i), upbnd))
         lp.changeRhs(i, infinity);

      // redundant lhs ?
      if (lp.lhs(i) >= -infinity && locnt == 0 && LE(lp.lhs(i), lobnd))
         lp.changeLhs(i, -infinity);

      // redundant constraint ?
      if (lp.rhs(i) >= infinity && lp.lhs(i) <= -infinity)
      {
         rem[i] = -1;
         num++;
         continue;
      }
#if 0

      {

            Real              y;

            if (upcnt < 2 || locnt < 2)
            {
               for( j = 0; j < row.size(); ++j )
               {
                  x = row.value(j);
                  k = row.index(j);

                  if (x > 0.0)
                  {
                     if (lp.lhs(i) > -infinity && lp.lower(k) > -infinity && upcnt < 2)
                     {
                        y = -infinity;
                        if (lp.upper(k) < infinity && upcnt < 1)
                           y = lp.upper(k) + (lp.lhs(i) - up) / x;
                        else if (lp.upper(k) >= infinity)
                           y = lp.lhs(i) - up;

                        if (y >= lp.lower(k))
                        {
                           lp.changeLower(k, -infinity);
                           break;
                        }
                     }
                     if (lp.rhs(i) < infinity && lp.upper(k) < infinity && locnt < 2)
                     {
                        y = infinity;

                        if (lp.lower(k) > -infinity && locnt < 1)
                           y = lp.lower(k) + (lp.rhs(i) - lo) / x;
                        else if (lp.lower(k) <= -infinity)
                           y = lp.rhs(i) - lo;

                        if (y <= lp.upper(k))
                        {
                           lp.changeUpper(k, infinity);
                           break;
                        }
                     }
                  }
                  else if (x < 0.0)
                  {
                     if (lp.lhs(i) >= -infinity && lp.upper(k) < infinity && upcnt < 2)
                     {
                        y = infinity;

                        if (lp.lower(k) > -infinity && upcnt < 1)
                           y = lp.lower(k) + (lp.lhs(i) - up) / x;
                        else if (lp.lower(k) <= -infinity)
                           y = -(lp.lhs(i) - up);

                        if (y <= lp.upper(k))
                        {
                           lp.changeUpper(k, infinity);
                           break;
                        }
                     }
                     if (lp.rhs(i) <= infinity && lp.lower(k) > -infinity
                        && locnt < 2)
                     {
                        y = -infinity;

                        if (lp.upper(k) < infinity && locnt < 1)
                           y = lp.upper(k) + (lp.rhs(i) - lo) / x;
                        else if (lp.upper(k) >= infinity)
                           y = -(lp.rhs(i) - lo);

                        if (y >= lp.lower(k))
                        {
                           lp.changeLower(k, -infinity);
                           break;
                        }
                     }
                  }
               }
            }
         }
#endif
   }
   if (num > 0)
   {
      lp.removeRows(rem.get_ptr());
      assert(lp.isConsistent());

      VERBOSE1({ std::cout << "SPxRedundantSM:\tremoved " << num
                           << " row(s)" << std::endl; });
   }
   return OKAY;
}

SPxSimplifier::Result SPxRedundantSM::simplify(SPxLP& lp)
{
   Result ret;

   if( (ret = treat_cols(lp)) == OKAY)
      ret = treat_rows(lp);

   return ret;
}

const Vector& SPxRedundantSM::unsimplifiedPrimal(const Vector& x)
{
   std::cout << "SPxRedundantSM::unsimplifiedPrimal() not implemented\n";

   abort();

   return x;
}

const Vector& SPxRedundantSM::unsimplifiedDual(const Vector& pi)
{
   std::cout << "SPxRedundantSM::unsimplifiedDual() not implemented\n";

   abort();

   return pi;
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
