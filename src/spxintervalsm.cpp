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
#pragma ident "@(#) $Id: spxintervalsm.cpp,v 1.2 2003/01/10 12:46:14 bzfkocht Exp $"

//#define DEBUGGING 1

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxintervalsm.h"

namespace soplex
{
SPxSimplifier::Result SPxIntervalSM::simplify(SPxLP& lp, Real eps, Real delta)
{
   METHOD("SPxIntervalSM::simplify()");

   //std::cout << lp << std::endl;

   Real maxval = infinity / 5.0;
   int  nzcnt  = 0;
   int  bdcnt  = 0;
   int  lrcnt  = 0;
   int  ojcnt  = 0;
  
   for(int i = 0; i < lp.nRows(); ++i)
   {
      // LHS 
      Real x = lp.lhs(i);

      if (x != 0.0 && isZero(x, delta))
      {
         lp.changeLhs(i, 0.0);
         lrcnt++;
      }
      else if (x > -infinity && x < -maxval)
      {
         lp.changeLhs(i, -infinity);
         lrcnt++;
      }
      else if (x <  infinity && x >  maxval)
      {
         lp.changeLhs(i,  infinity);
         lrcnt++;
      }

      // RHS
      x = lp.rhs(i);

      if (x != 0.0 && isZero(x, delta))
      {
         lp.changeRhs(i, 0.0);
         lrcnt++;
      }
      else if (x > -infinity && x < -maxval)
      {
         lp.changeRhs(i, -infinity);
         lrcnt++;
      }
      else if (x <  infinity && x >  maxval)
      {
         lp.changeRhs(i,  infinity);
         lrcnt++;
      }
   }
   for(int i = 0; i < lp.nCols(); ++i)
   {
      // lower bound
      Real lo = lp.lower(i);

      if (lo != 0.0 && isZero(lo, delta))
      {
         lp.changeLower(i, 0.0);
         bdcnt++;
      }
      else if (lo > -infinity && lo < -maxval)
      {
         lp.changeLower(i, -infinity);
         bdcnt++;
      }
      else if (lo <  infinity && lo >  maxval)
      {
         lp.changeLower(i,  infinity);
         bdcnt++;
      }

      // upper bound
      Real up = lp.lower(i);

      if (up != 0.0 && isZero(up, delta))
      {
         lp.changeUpper(i, 0.0);
         bdcnt++;
      }
      else if (up > -infinity && up < -maxval)
      {
         lp.changeUpper(i, -infinity);
         bdcnt++;
      }
      else if (up <  infinity && up >  maxval)
      {
         lp.changeUpper(i,  infinity);
         bdcnt++;
      }

      // fixed columns will be eliminated later
      if (NE(lo, up))
      {
         lo = fabs(lo);
         up = fabs(up);

         Real absbnd = lo > up ? lo : up;

         // Anything bigger than delta is ok.
         if (absbnd < 1.0)
            absbnd = 1.0;

         // Non-zeros
         SVector& col = lp.colVector_w(i);
         Real           x;
         int            j = 0;
      
         while(j < col.size())
         {
            x = fabs(col.value(j));

            if (isZero(x, eps) || isZero(x * absbnd, delta))
            {
               SVector& row = lp.rowVector_w(col.index(j));

               // This changes col.size();
               row.remove(row.number(i));
               col.remove(j);           

               VERBOSE3({ std::cout << "removed element x=" << x 
                                    << " absbnd= " << absbnd 
                                    << std::endl; });
               nzcnt++;
            }         
            else
            {
               if (x > maxval)
                  std::cerr << "Warning! Big value " << x << std::endl;

               j++;
            }
         }
      }

      // objective
      Real x = lp.obj(i);

      if (x != 0.0 && isZero(x, delta))
      {
         lp.changeObj(i, 0.0);
         ojcnt++;
      }
      else if (x > -infinity && x < -maxval)
      {
         lp.changeObj(i, -infinity);
         ojcnt++;
      }
      else if (x <  infinity && x >  maxval)
      {
         lp.changeObj(i,  infinity);
         ojcnt++;
      }
   }
   if (nzcnt > 0)
   {
      VERBOSE1({ std::cout << "SPxIntervalSM:\tremoved " << nzcnt
                           << " non-zeros" << std::endl; });
   }
   if (lrcnt > 0)
   {
      VERBOSE1({ std::cout << "SPxIntervalSM:\tcorrected " << lrcnt
                           << " LHS/RHS" << std::endl; });
   }
   if (bdcnt > 0)
   {
      VERBOSE1({ std::cout << "SPxIntervalSM:\tcorrected " << bdcnt
                           << " bounds" << std::endl; });
   }
   if (ojcnt > 0)
   {
      VERBOSE1({ std::cout << "SPxIntervalSM:\tcorrected " << ojcnt
                           << " objective function coefficents" << std::endl; });
   }
   assert(lp.isConsistent());

   //std::cout << lp << std::endl;
   
   return OKAY;
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
