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
#pragma ident "@(#) $Id: spxrem1sm.cpp,v 1.15 2003/01/05 19:03:17 bzfkocht Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxrem1sm.h"
#include "dataarray.h"

namespace soplex
{
Real SPxRem1SM::epsilonSimplifier() const
{
   return 1e-14;
}

void SPxRem1SM::fixColumn(SPxLP& lp, int i)
{
   METHOD( "SPxRem1SM::fixColumn" );

   assert(EQ(lp.lower(i), lp.upper(i), epsilonSimplifier()));

   Real x = lp.upper(i);

   VERBOSE3({ std::cout << "Fixed column " << i 
                        << " lower= " << std::setprecision(16) << lp.lower(i)
                        << " upper= " << std::setprecision(16) << lp.upper(i)
                        << std::endl; });

   pval.add(cperm[i], x);
   
   if (isNotZero(x, epsilonSimplifier()))
   {
      const SVector& col = lp.colVector(i);

      for(int j = 0; j < col.size(); ++j)
      {
         int k = col.index(j);
         
         if (lp.rhs(k) < infinity)
         {
            Real y     = x * col.value(j);
            Real scale = fabs(lp.rhs(k)) > fabs(y) ? fabs(lp.rhs(k)) : fabs(y);
            Real rhs   = (lp.rhs(k) / scale) - (y / scale);
                  
            if (isZero(rhs, epsilonSimplifier()))
               rhs = 0.0;
            else
               rhs *= scale;

            VERBOSE3({ std::cout << "\trhs " << k 
                                 << " r= " << std::setprecision(16) << rhs 
                                 << " rhs= " << std::setprecision(16) << lp.rhs(k) 
                                 << " x= " << std::setprecision(16) << col.value(j) 
                                 << std::endl; });

            lp.changeRhs(k, rhs);
         }
         if (lp.lhs(k) > -infinity)
         {
            Real y     = x * col.value(j);
            Real scale = fabs(lp.lhs(k)) > fabs(y) ? fabs(lp.lhs(k)) : fabs(y);
            Real lhs   = (lp.lhs(k) / scale) - (y / scale);

            if (isZero(lhs, epsilonSimplifier()))
               lhs = 0.0;
            else
               lhs *= scale;
                  
            VERBOSE3({ std::cout << "\tlhs " << k 
                                 << " l= " << std::setprecision(16) << lhs 
                                 << " lhs= " << std::setprecision(16) << lp.lhs(k) 
                                 << " x= " << std::setprecision(16) << col.value(j) 
                                 << std::endl; });

            lp.changeLhs(k, lhs);
         }
      }
   }
}

SPxSimplifier::Result SPxRem1SM::redundantRows(SPxLP& lp, bool& again)
{
   DataArray<int> rem(lp.nRows());
   int            num = 0;

   for( int i = 0; i < lp.nRows(); ++i )
   {
      //      std::cout << "row " << i << std::endl;

      const SVector& row    = lp.rowVector(i);       
      Real           upbnd  = 0.0;
      Real           lobnd  = 0.0;
      int            upcnt  = 0;
      int            locnt  = 0;

      rem[i] = 0;

      for(int j = 0; j < row.size(); ++j )
      {
         Real x = row.value(j);
         int  k = row.index(j);
         
         // std::cout << "\t\tx= " << x << " lower= " << lp.lower(k) << " upper= " << lp.upper(k) << std::endl;

         assert(isNotZero(x, epsilonSimplifier()));

         if (x > 0.0)
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
         else if (x < 0.0)
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
      if (  (LT(lp.rhs(i), lobnd, epsilonSimplifier()) && locnt == 0) 
         || (GT(lp.lhs(i), upbnd, epsilonSimplifier()) && upcnt == 0))
      {
         VERBOSE3({ std::cout << "infeasible row " << i 
                              << " lo= " << lobnd
                              << " up= " << upbnd 
                              << " lhs= " << lp.lhs(i) 
                              << " rhs= " << lp.rhs(i)
                              << std::endl; });
         return INFEASIBLE;
      }
      // forcing equality constraint ?
      if (EQ(lp.lhs(i), lp.rhs(i), epsilonSimplifier()))
      {
         // all fixed on upper bound ?
         if (upcnt == 0 && EQ(lp.rhs(i), upbnd, epsilonSimplifier()))
         {
            VERBOSE3({ std::cout << "\trhs fixed on upbnd row " << i
                                 << " rhs= " << lp.rhs(i)
                                 << " up= " << upbnd 
                                 << std::endl; });

            for(int j = 0; j < row.size(); ++j )
            {
               Real x = row.value(j);
               int  k = row.index(j);

               assert(isNotZero(x, epsilonSimplifier()));
               assert(GE(lp.lower(k), 0.0) || LE(lp.upper(k), 0.0));

               if (x > 0.0 && GE(lp.lower(k), 0.0))
                  lp.changeLower(k, lp.upper(k));
               else
                  lp.changeUpper(k, lp.lower(k));
            }
            rem[i] = -1;
            num++;
            continue;
         }
         // all fixed on lower bound ?
         if (locnt == 0 && EQ(lp.lhs(i), lobnd, epsilonSimplifier()))
         {
            VERBOSE3({ std::cout << "\trhs fixed on lowbnd row " << i
                                 << " lhs= " << lp.lhs(i)
                                 << " lo= " << lobnd 
                                 << std::endl; });

            for(int j = 0; j < row.size(); ++j )
            {
               Real x = row.value(j);
               int  k = row.index(j);

               assert(isNotZero(x, epsilonSimplifier()));
               assert(GE(lp.lower(k), 0.0) || LE(lp.upper(k), 0.0));

               if (x > 0.0 && GE(lp.lower(k), 0.0))
                  lp.changeUpper(k, lp.lower(k));
               else
                  lp.changeLower(k, lp.upper(k));
            }
            rem[i] = -1;
            num++;
            continue;
         }
      }

      // redundant rhs ?
      if (lp.rhs(i) <=  infinity && upcnt == 0 && GE(lp.rhs(i), upbnd, epsilonSimplifier()))
      {
         VERBOSE3({ std::cout << "\tredundant rhs row " << i
                              << " rhs= " << lp.rhs(i)
                              << " up= " << upbnd 
                              << std::endl; });

         lp.changeRhs(i, infinity);
      }
      // redundant lhs ?
      if (lp.lhs(i) >= -infinity && locnt == 0 && LE(lp.lhs(i), lobnd, epsilonSimplifier()))
      {
         VERBOSE3({ std::cout << "\tredundant lhs row " << i
                              << " lhs= " << lp.lhs(i)
                              << " lo= " << lobnd 
                              << std::endl; });

         lp.changeLhs(i, -infinity);
      }
      // unconstraint constraints are also handled by simpleRows
      // but since they might come up here we do it again.
      if (lp.rhs(i) >= infinity && lp.lhs(i) <= -infinity)
      {
         VERBOSE3({ std::cout << "Unconstraint row " << i << " removed" << std::endl; });

         rem[i] = -1;
         num++;
         continue;
      }
#if 0
      if (upcnt <= 1 || locnt <= 1)
      {
         for(int j = 0; j < row.size(); ++j )
         {
            Real x = row.value(j);
            int  k = row.index(j);

            assert(isNotZero(x));

            if (x > 0.0)
            {
               if (lp.lhs(i) > -infinity && lp.lower(k) > -infinity && upcnt <= 1)
               {
                  Real y = -infinity;

                  if (lp.upper(k) < infinity && upcnt == 0)
                     y = lp.upper(k) + (lp.lhs(i) - upbnd) / x;
                  else if (lp.upper(k) >= infinity)
                     y = lp.lhs(i) - upbnd;

                  if (isZero(y, epsilonSimplifier()))
                     y = 0.0;

                  if (GE(y, lp.lower(k)))
                  {
                     VERBOSE3({ std::cout << "dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; });

                     lp.changeLower(k, -infinity);
                     break;
                  }
               }
               if (lp.rhs(i) < infinity && lp.upper(k) < infinity && locnt <= 1)
               {
                  Real y = infinity;
                  
                  if (lp.lower(k) > -infinity && locnt == 0)
                     y = lp.lower(k) + (lp.rhs(i) - lobnd) / x;
                  else if (lp.lower(k) <= -infinity)
                     y = lp.rhs(i) - lobnd;
                  
                  if (isZero(y, epsilonSimplifier()))
                     y = 0.0;

                  if (LE(y, lp.upper(k)))
                  {
                     VERBOSE3({ std::cout << "dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; });

                     lp.changeUpper(k, infinity);
                     break;
                  }
               }
            }
            else if (x < 0.0)
            {
               if (lp.lhs(i) >= -infinity && lp.upper(k) < infinity && upcnt <= 1)
               {
                  Real y = infinity;
                  
                  if (lp.lower(k) > -infinity && upcnt == 0)
                     y = lp.lower(k) + (lp.lhs(i) - upbnd) / x;
                  else if (lp.lower(k) <= -infinity)
                     y = -(lp.lhs(i) - upbnd);
                  
                  if (isZero(y, epsilonSimplifier()))
                     y = 0.0;

                  if (LE(y, lp.upper(k)))
                  {
                     VERBOSE3({ std::cout << "dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; });

                     lp.changeUpper(k, infinity);
                     break;
                  }
               }
               if (lp.rhs(i) <= infinity && lp.lower(k) > -infinity && locnt <= 1)
               {
                  Real y = -infinity;

                  if (lp.upper(k) < infinity && locnt == 0)
                     y = lp.upper(k) + (lp.rhs(i) - lobnd) / x;
                  else if (lp.upper(k) >= infinity)
                     y = -(lp.rhs(i) - lobnd);
                  
                  if (isZero(y, epsilonSimplifier()))
                     y = 0.0;

                  if (GE(y, lp.lower(k)))
                  {
                     VERBOSE3({ std::cout << "dominated bound row " << i
                                          << " col " << k
                                          << " removed y= " << y
                                          << " lower= " << lp.lower(k)
                                          << std::endl; });

                     lp.changeLower(k, -infinity);
                     break;
                  }
               }
            }
         }
      }
#endif
#if 0 // org

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
                        if (lp.upper(k) < infinity && upcnt == 0)
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

                        if (lp.lower(k) > -infinity && locnt == 0)
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
   again = removeRows(lp, rem, num, "redundantRows");

   return OKAY;
}

SPxSimplifier::Result SPxRem1SM::redundantCols(SPxLP& lp, bool& again)
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

      if (NE(lp.upper(i), lp.lower(i), epsilonSimplifier()))
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

            assert(isNotZero(x, epsilonSimplifier()));

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
         if (locnt == 0 && LT(x, 0.0, epsilonSimplifier()))
         {
            if (lp.lower(i) <= -infinity)
               return UNBOUNDED;           // LP is unbounded

            lp.changeUpper(i, lp.lower(i));
         }
         // max  3 x
         // s.t. 5 x >= 8
         else if (upcnt == 0 && GT(x, 0.0, epsilonSimplifier()))
         {
            if (lp.upper(i) >= infinity)
               return UNBOUNDED;           // LP is unbounded

            lp.changeLower(i, lp.upper(i));
         }
         else if (isZero(x, epsilonSimplifier()))
         {
            //TODO free variable with zero objective can be removed, or?
            //see simpleCols

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
      if (EQ(lp.lower(i), lp.upper(i), epsilonSimplifier()))
      {
         fixColumn(lp, i);
         rem[i] = -1;
         num++;
      }
   }
   again = removeCols(lp, rem, num, "redundantCols");

   return OKAY;
}

/* Handle rows -------------------------------------------------------
 */
SPxSimplifier::Result SPxRem1SM::simpleRows(SPxLP& lp, bool& again)
{
   DataArray<int> rem(lp.nRows());
   DataArray<int> tmp(lp.nRows());
   int            num = 0;
   int            i;

   for(i = 0; i < lp.nRows(); ++i)
   {
      const SVector& row = lp.rowVector(i);
      rem[i]             = 0;

      // infeasible range row
      if (LT(lp.rhs(i), lp.lhs(i), epsilonSimplifier()))
      {
         VERBOSE3({ std::cout << "Infeasible row " << i 
                              <<"  lhs= " << lp.lhs(i) 
                              << " rhs= " << lp.rhs(i) 
                              << std::endl; });
         return INFEASIBLE;
      }
      // empty row ?
      if (row.size() == 0)
      {
         VERBOSE3({ std::cout << "Empty row " << i; });

         if (LT(lp.rhs(i), 0.0, epsilonSimplifier()) || GT(lp.lhs(i), 0.0, epsilonSimplifier()))
         {
            VERBOSE3({ std::cout << " infeasible lhs= " << lp.lhs(i) 
                                 << " rhs= " << lp.rhs(i) << std::endl; });
            return INFEASIBLE;
         }         
         VERBOSE3({ std::cout << " removed" << std::endl; });

         rem[i] = -1;
         num++;
         continue;
      }
      // unconstraint constraint ?
      if (lp.rhs(i) >= infinity && lp.lhs(i) <= -infinity)
      {
         VERBOSE3({ std::cout << "Unconstraint row " << i << " removed" << std::endl; });

         rem[i] = -1;
         num++;
         continue;
      }
      // row singleton
      if (row.size() == 1)
      {
         Real x = row.value(0);
         int  j = row.index(0);
         Real up;
         Real lo;

         VERBOSE3({ std::cout << "Row singleton " << i 
                              << " x= " << x 
                              << " lhs= " << lp.lhs(i) 
                              << " rhs= " << lp.rhs(i); });

         if (GT(x, 0.0, epsilonSimplifier()))           // x > 0
         {
            up = (lp.rhs(i) >=  infinity) ?  infinity : (lp.rhs(i) / x);
            lo = (lp.lhs(i) <= -infinity) ? -infinity : (lp.lhs(i) / x);
         }
         else if (LT(x, 0.0, epsilonSimplifier()))      // x < 0
         {
            lo = (lp.rhs(i) >=  infinity) ? -infinity : (lp.rhs(i) / x);
            up = (lp.lhs(i) <= -infinity) ?  infinity : (lp.lhs(i) / x);
         }
         else if (LT(lp.rhs(i), 0.0, epsilonSimplifier()) || GT(lp.lhs(i), 0.0, epsilonSimplifier()))  
         {
            // x == 0 rhs/lhs != 0
            VERBOSE3({ std::cout << " infeasible" << std::endl; });

            return INFEASIBLE;
         }
         else                     // x == 0
         {
            lo = lp.lower(j);
            up = lp.upper(j);
         }
         rem[i] = -1;
         num++;
         
         if (isZero(lo, epsilonSimplifier()))
            lo = 0.0;
         if (isZero(up, epsilonSimplifier()))
            up = 0.0;
         
         assert(LE(lp.lower(j), lp.upper(j)));

         VERBOSE3({ std::cout << " removed lo= " << lo
                              << " up= " << up
                              << " lower= " << lp.lower(j)
                              << " upper= " << lp.upper(j)
                              << std::endl; });

         if (LT(up, lp.upper(j), epsilonSimplifier()))
            lp.changeUpper(j, up);
         if (GT(lo, lp.lower(j), epsilonSimplifier()))
            lp.changeLower(j, lo);

         assert(LE(lp.lower(j), lp.upper(j), epsilonSimplifier()));
      }
   }
   again = removeRows(lp, rem, num, "simpleRows");

   return OKAY;
}

/* Handle columns -----------------------------------------------------
 */
SPxSimplifier::Result SPxRem1SM::simpleCols(SPxLP& lp, bool& again)
{
   DataArray<int> rem(lp.nCols());
   DataArray<int> tmp(lp.nCols());
   int            num = 0;
   int            i;

   for(i = 0; i < lp.nCols(); ++i)
   {
      const SVector& col = lp.colVector(i);
      rem[i]             = 0;

      // Empty column ? 
      if (col.size() == 0)
      {
         VERBOSE3({ std::cout << "Empty column " << i 
                              << " maxObj= " << lp.maxObj(i)
                              << " lower= " << lp.lower(i)
                              << " upper= " << lp.upper(i); });

         if (GT(lp.maxObj(i), 0.0, epsilonSimplifier()))
         {
            if (lp.upper(i) >= infinity)
            {
               VERBOSE3({ std::cout << " unbounded" << std::endl; });

               return UNBOUNDED;
            }
            pval.add(cperm[i], lp.upper(i));
         }
         else if (LT(lp.maxObj(i), 0.0, epsilonSimplifier()))
         {
            if (lp.lower(i) <= -infinity)
            {
               VERBOSE3({ std::cout << " unbounded" << std::endl; });

               return UNBOUNDED;
            }
            pval.add(cperm[i], lp.lower(i));
         }
         else 
         {
            assert(isZero(lp.maxObj(i), epsilonSimplifier()));
            // any value within the bounds is ok
            if (lp.lower(i) > -infinity)
               pval.add(cperm[i], lp.lower(i));
            else if (lp.upper(i) < infinity)
               pval.add(cperm[i], lp.upper(i));
            else
               pval.add(cperm[i], 0.0);
         }
         VERBOSE3({ std::cout << " removed" << std::endl; });

         rem[i] = -1;
         num++;
         continue;
      }

      // infeasible bounds ?
      if (GT(lp.lower(i), lp.upper(i), epsilonSimplifier()))
      {
         VERBOSE3({ std::cout << "Infeasible bounds column " << i 
                              << " lower= " << lp.lower(i)
                              << " upper= " << lp.upper(i)
                              << std::endl; });
         return INFEASIBLE;
      }
      // Fixed column ?
      if (EQ(lp.lower(i), lp.upper(i), epsilonSimplifier()))
      {
         fixColumn(lp, i);
         rem[i] = -1;
         num++;
      }
      // a lot of this is also done in redundantRows.
      else if (col.size() == 1)
      {
         Real x = col.value(0);
         int  j = col.index(0);
         
         assert(isNotZero(x));

         if (x > 0.0)
         {
            // max -3 x
            // s.t. 5 x <= 8
            // l <= x
            if (lp.lhs(j) <= -infinity && lp.rhs(j) < infinity && LE(lp.maxObj(i), 0.0) && lp.lower(i) > -infinity)
            {
               lp.changeUpper(i, lp.lower(i));
               fixColumn(lp, i);
               rem[i] = -1;
               num++;
               continue;
            }
            // max  3 x
            // s.t. 5 x >= 8
            // x <= u
            if (lp.lhs(j) > -infinity && lp.rhs(j) >= infinity && GE(lp.maxObj(i), 0.0) && lp.upper(i) < infinity)
            {
               lp.changeLower(i, lp.upper(i));
               fixColumn(lp, i);
               rem[i] = -1;
               num++;
               continue;
            }
         }
      }
#if 0 // a lot of this is allready done in redundantRows.

      // Column singleton without objective
      else if (col.size() == 1 && isZero(lp.maxObj(i), epsilonSimplifier()))
      {
         Real x = col.value(0);
         int  j = col.index(0);
         
         Real           up;
         Real           lo;


         if (GT(x, 0.0, epsilonSimplifier()))
         {
            if (lp.lower(i) > -infinity)
               // TODO here should be checked if rhs < infinity
               up = lp.rhs(j) - lp.lower(i) * x;
            else
               up = infinity;
            
            if (lp.upper(i) < infinity)
               lo = lp.lhs(j) - lp.upper(i) * x;
            else
               lo = -infinity;
         }
         else if (LT(x, 0.0, epsilonSimplifier()))
         {
            if (lp.lower(i) > -infinity)
               lo = lp.lhs(j) - lp.lower(i) * x;
            else
               lo = -infinity;
            
            if (lp.upper(i) < infinity)
               up = lp.rhs(j) - lp.upper(i) * x;
            else
               up = infinity;
         }
         else
         {
            up = lp.rhs(j);
            lo = lp.lhs(j);
         }
         
         if (isZero(lo, epsilonSimplifier()))
            lo = 0.0;
         if (isZero(up, epsilonSimplifier()))
            up = 0.0;
         
         lp.changeRange(j, lo, up);
         rem[i] = -1;
         num++;
      }
#endif // 0
   }
   again = removeCols(lp, rem, num, "simpleCols");

   return OKAY;
}

bool SPxRem1SM::removeRows(SPxLP& lp, DataArray<int>& rem, int num, const char* msg) 
{
   if (num == 0)
      return false;

   int i;

   assert(rem.size() == lp.nRows());

   lp.removeRows(rem.get_ptr());

   DataArray<int> tmp(rperm);
      
   for(i = 0; i < rem.size(); ++i)
      if (rem[i] >= 0)
         rperm[rem[i]] = tmp[i];
   
   for(i = lp.nRows(); i < lp.nRows() + num; ++i)
      rperm[i] = -1;
   
   assert(lp.isConsistent());
   
   VERBOSE1({ std::cout << msg << " simplifier removed " 
                        << num << " row(s)" << std::endl; });   
   return true;
}

bool SPxRem1SM::removeCols(SPxLP& lp, DataArray<int>& rem, int num, const char* msg) 
{
   if (num == 0)
      return false;

   int i;

   assert(rem.size() == lp.nCols());

   lp.removeCols(rem.get_ptr());
      
   DataArray<int> tmp(cperm);
      
   for(i = 0; i < rem.size(); ++i)
      if (rem[i] >= 0)
         cperm[rem[i]] = tmp[i];
   
   for(i = lp.nCols(); i < lp.nCols() + num; ++i)
      cperm[i] = -1;
#if 0
     for(i = 0; i < cperm.size(); ++i)
      std::cout << i << " -> " << cperm[i] << " " << tmp[i] << " " << rem[i] << "\n";
#endif
   assert(lp.isConsistent());
   
   VERBOSE1({ std::cout << msg << " simplifier removed " 
                        << num << " column(s)" << std::endl; });
   return true;
}

SPxSimplifier::Result SPxRem1SM::simplify(SPxLP& lp)
{
   Result ret;
   bool   again;
   bool   rcagain;
   bool   rragain;
   int    i;

   prim.reDim(lp.nCols());
   dual.reDim(lp.nRows());

   cperm.reSize(lp.nCols());
   rperm.reSize(lp.nRows());

   for(i = 0; i < lp.nRows(); ++i)
      rperm[i] = i;

   for(i = 0; i < lp.nCols(); ++i)
      cperm[i] = i;

   //std::cout << lp << std::endl;

   do
   {
      again = false;
 
      if( (ret = simpleRows(lp, again)) != OKAY)
         return ret;

      // std::cout << lp << std::endl;

      if( (ret = simpleCols(lp, again)) != OKAY)
         return ret;

      //std::cout << lp << std::endl;
      //abort();
   }
   while(again);

   //std::cout << lp << std::endl;
#if 1
   if( (ret = redundantCols(lp, rcagain)) != OKAY)
      return ret;

   if( (ret = redundantRows(lp, rragain)) != OKAY)
      return ret;

   if (rcagain || rragain)
   {
      if( (ret = simpleRows(lp, again)) != OKAY)
         return ret;

      if( (ret = simpleCols(lp, again)) != OKAY)
         return ret;
   }
#endif
   return OKAY;
}

const Vector& SPxRem1SM::unsimplifiedPrimal(const Vector& x)
{
   assert(x.dim()               <= cperm.size());
   assert(x.dim() + pval.size() == prim.dim());

   int i;

   for(i = 0; i < x.dim(); ++i)
      prim[cperm[i]] = x[i];

   assert(cperm.size() == x.dim() || cperm[i] < 0);

   for(i = 0; i < pval.size(); ++i)
      prim[pval.index(i)] = pval.value(i);

   return prim;
}

const Vector& SPxRem1SM::unsimplifiedDual(const Vector& pi)
{
   std::cout << "SPxRem1SM::getDual() not implemented\n";

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
