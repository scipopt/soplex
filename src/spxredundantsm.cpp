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
#pragma ident "@(#) $Id: spxredundantsm.cpp,v 1.18 2002/05/15 13:38:44 bzfpfend Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxredundantsm.h"
#include "dataarray.h"

namespace soplex
{

int SPxRedundantSM::treat_cols()
{
   DataArray < int > rem(lp->nCols());
   int               num = 0;
   int               j;
   int               k;
   Real              x;

   for( int i = 0; i < lp->nCols(); ++i )
   {
      rem[i] = 0;

      const SVector& col = lp->colVector(i);

      if (lp->upper(i) != lp->lower(i))
      {
         int upcnt = 0;
         int locnt = 0;

         for( j = 0; (j < col.size()) && (upcnt == 0 || locnt == 0); ++j )
         {            
            x = col.value(j);
            k = col.index(j);

            if (x > 0.0)
            {
               upcnt += (lp->rhs(k) <  infinity) ? 1 : 0;
               locnt += (lp->lhs(k) > -infinity) ? 1 : 0;
            }
            else if (x < 0.0)
            {
               locnt += (lp->rhs(k) <  infinity) ? 1 : 0;
               upcnt += (lp->lhs(k) > -infinity) ? 1 : 0;
            }
         }
         x = lp->maxObj(i);

         if (locnt == 0 && x < 0.0)
         {
            if (lp->lower(i) <= -infinity)
               return 1;           // LP is unbounded
            lp->changeUpper(i, lp->lower(i));
         }
         else if (upcnt == 0 && x > 0.0)
         {
            if (lp->upper(i) >= infinity)
               return 1;           // LP is unbounded
            lp->changeLower(i, lp->upper(i));
         }
         else if (x == 0.0)
         {
            upcnt += (lp->upper(i) <  infinity) ? 1 : 0;
            locnt += (lp->lower(i) > -infinity) ? 1 : 0;

            if (locnt == 0)
            {
               lp->changeUpper(i, infinity);

               for( j = 0; j < col.size(); ++j )
               {
                  if (col.value(j) < 0.0)
                     lp->changeRhs(col.index(j), infinity);
                  else
                     lp->changeLhs(col.index(j), -infinity);
               }
            }
            if (upcnt == 0)
            {
               lp->changeLower(i, -infinity);

               for( j = 0; j < col.size(); ++j )
               {
                  if (col.value(j) > 0.0)
                     lp->changeRhs(col.index(j), infinity);
                  else
                     lp->changeLhs(col.index(j), -infinity);
               }
            }
         }
      }
      // remove fixed variables
      if (lp->upper(i) == lp->lower(i))
      {
         x      = lp->upper(i);
         rem[i] = -1;
         num++;

         if (x != 0.0)
         {
            for( j = 0; j < col.size(); ++j )
            {
               k = col.index(j);

               if (lp->rhs(k) <  infinity)
                  lp->changeRhs(k, lp->rhs(k) - x * col.value(j));
               if (lp->lhs(k) > -infinity)
                  lp->changeLhs(k, lp->lhs(k) - x * col.value(j));
            }
            delta += x * lp->obj(i);
         }
      }
   }
   if (num > 0)
   {
      lp->removeCols(rem.get_ptr());
      assert(lp->isConsistent());

      VERBOSE1({ std::cout << "SPxRedundantSM: removed " << num
                           << " column(s)" << std::endl; });
   }
   return 0;
}


int SPxRedundantSM::treat_rows()
{
   DataArray < int > rem(lp->nRows());
   int               num = 0;
   int               j;
   int               k;
   Real              x;
   Real              y;

   for( int i = 0; i < lp->nRows(); ++i )
   {
      if (lp->rhs(i) < infinity || lp->lhs(i) > -infinity)
      {
         rem[i] = 0;

         const SVector& row = lp->rowVector(i);
         
         Real up    = 0.0;
         Real lo    = 0.0;
         int  upcnt = 0;
         int  locnt = 0;

         for( j = 0; j < row.size(); ++j )
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

         if (  ((GE(lp->rhs(i), up) && upcnt == 0) || lp->rhs(i) >=  infinity)
            && ((LE(lp->lhs(i), lo) && locnt == 0) || lp->lhs(i) <= -infinity))
         {
            rem[i] = -1;
            num++;
         }
         else if ((LT(lp->rhs(i), lo) && locnt == 0)
            ||    (GT(lp->lhs(i), up) && upcnt == 0))
            return -1; // infeasible
#if 0
         else
         {
            /*
                if (LE(lp->lhs(i), lo) && locnt <= 0)
                    lp->changeLhs(i, -infinity);
                else if (GE(lp->rhs(i), up) && upcnt <= 0)
                    lp->changeRhs(i, infinity);
                else
             */
            if (upcnt < 2 || locnt < 2)
            {
               for( j = 0; j < row.size(); ++j )
               {
                  x = row.value(j);
                  k = row.index(j);

                  if (x > 0.0)
                  {
                     if (lp->lhs(i) > -infinity && lp->lower(k) > -infinity
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
                     if (lp->rhs(i) < infinity && lp->upper(k) < infinity
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
                  else if (x < 0.0)
                  {
                     if (lp->lhs(i) >= -infinity && lp->upper(k) < infinity
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
                     if (lp->rhs(i) <= infinity && lp->lower(k) > -infinity
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
#endif
      }
      else
         rem[i] = -1;
   }
   if (num > 0)
   {
      lp->removeRows(rem.get_ptr());
      assert(lp->isConsistent());

      VERBOSE1({ std::cout << "SPxRedundantSM:\tremoved " << num
                           << " row(s)" << std::endl; });
   }
   return 0;
}

int SPxRedundantSM::simplify()
{
   int ret;

      if( (ret = treat_cols()) == 0)
         ret = treat_rows();

   return ret;
}

void SPxRedundantSM::unsimplify()
{
   std::cout << "SPxRedundantSM::unsimplify() not implemented\n";
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
