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
#pragma ident "@(#) $Id: spxrem1sm.cpp,v 1.14 2002/12/17 10:20:06 bzfkocht Exp $"

//#define DEBUGGING 1

#include <iostream>

#include "spxdefines.h"
#include "spxrem1sm.h"
#include "dataarray.h"

namespace soplex
{
// -1 = infeasible
//  0 = ok
//  1 = unbounded
int SPxRem1SM::simplify()
{
   bool  cont;
   int   num;
   Real  up;
   Real  lo;
   DataArray<int> rem;

   do
   {
      /* Handle rows -------------------------------------------------------
       */
      cont = false;
      num  = 0;

      rem.reSize(lp->nRows());

      for(int i = 0; i < lp->nRows(); ++i)
      {
         const SVector& row = lp->rowVector(i);
         rem[i]             = 0;

         // empty row ?
         if (row.size() == 0)
         {
            if (lp->rhs(i) < 0.0 || lp->lhs(i) > 0.0)
               return -1;

            rem[i] = -1;
            num++;
         }
         // unconstraint constraint ?
         else if (lp->rhs(i) >= infinity && lp->lhs(i) <= -infinity)
         {
            rem[i] = -1;
            num++;
         }
         // row singleton
         else if (row.size() == 1)
         {
            Real x = row.value(0);
            int  j = row.index(0);

            if (GT(x, 0.0))
            {
               up = lp->rhs(i) / x;
               lo = lp->lhs(i) / x;
            }
            else if (LT(x, 0.0))
            {
               lo = lp->rhs(i) / x;
               up = lp->lhs(i) / x;
            }
            else if (LT(lp->rhs(i), 0.0) || GT(lp->lhs(i), 0.0))
               return -1;
            else
            {
               lo = lp->lower(j);
               up = lp->upper(j);
            }
            rem[i] = -1;
            num++;

            if (isZero(lo))
               lo = 0.0;
            if (isZero(up))
               up = 0.0;

            if (up < lp->upper(j))
               lp->changeUpper(j, up);
            if (lo > lp->lower(j))
               lp->changeLower(j, lo);
         }
      }
      if (num > 0)
      {
         cont = true;
         lp->removeRows(rem.get_ptr());
         VERBOSE1({ std::cout << "SPxRem1SM:\tremoved " << num
                              << " row(s)" << std::endl; });
         assert(lp->isConsistent());
      }

      /* Handle columns -----------------------------------------------------
       */
      num = 0;
      rem.reSize(lp->nCols());

      for(int i = lp->nCols() - 1; i >= 0; --i)
      {
         const SVector& col = lp->colVector(i);
         rem[i]             = 0;

         // Empty column ? 
         if (col.size() == 0)
         {
            if (GT(lp->maxObj(i), 0.0))
            {
               if (lp->upper(i) >= infinity)
                  return 1;

               delta += lp->upper(i) * lp->obj(i);
            }
            else if (LT(lp->maxObj(i), 0.0))
            {
               if (lp->lower(i) <= -infinity)
                  return 1;

               delta += lp->lower(i) * lp->obj(i);
            }
            rem[i] = -1;
            num++;
         }
         // Fixed column ?
         else if (EQ(lp->upper(i), lp->lower(i)))
         {
            rem[i] = -1;
            num++;

            Real x = lp->upper(i);

            if (isNotZero(x))
            {
               for(int j = 0; j < col.size(); ++j)
               {
                  int k = col.index(j);

                  if (lp->rhs(k) < infinity)
                     lp->changeRhs(k, lp->rhs(k) - x * col.value(j));

                  if (lp->lhs(k) > -infinity)
                     lp->changeLhs(k, lp->lhs(k) - x * col.value(j));
               }
               delta += x * lp->obj(i);
            }
         }
         // Column singleton without objective
         else if (col.size() == 1 && isZero(lp->maxObj(i)))
         {
            Real x = col.value(0);
            int  j = col.index(0);

            if (GT(x, 0.0))
            {
               if (lp->lower(i) > -infinity)
                  up = lp->rhs(j) - lp->lower(i) * x;
               else
                  up = infinity;

               if (lp->upper(i) < infinity)
                  lo = lp->lhs(j) - lp->upper(i) * x;
               else
                  lo = -infinity;
            }
            else if (LT(x, 0.0))
            {
               if (lp->lower(i) > -infinity)
                  lo = lp->lhs(j) - lp->lower(i) * x;
               else
                  lo = -infinity;

               if (lp->upper(i) < infinity)
                  up = lp->rhs(j) - lp->upper(i) * x;
               else
                  up = infinity;
            }
            else
            {
               up = lp->rhs(j);
               lo = lp->lhs(j);
            }

            if (isZero(lo))
               lo = 0.0;
            if (isZero(up))
               up = 0.0;

            lp->changeRange(j, lo, up);
            rem[i] = -1;
            num++;
         }
      }
      if (num > 0)
      {
         cont = true;
         lp->removeCols(rem.get_ptr());
         VERBOSE1({ std::cout << "SPxRem1SM:\tremoved " << num
                              << " column(s)" << std::endl; });
         assert(lp->isConsistent());
      }
   }
   while(cont);

   return 0;
}

void SPxRem1SM::unsimplify()
{
   std::cout << "SPxRem1SM::unsimplify() not implemented\n";
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
