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
#pragma ident "@(#) $Id: spxintervalsm.cpp,v 1.1 2003/01/05 19:03:17 bzfkocht Exp $"

//#define DEBUGGING 1

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxintervalsm.h"

namespace soplex
{
SPxSimplifier::Result SPxIntervalSM::simplify(SPxLP& lp)
{
   METHOD("SPxIntervalSM::simplify()");

   Real maxval = infinity / 5.0;
   int  nzcnt  = 0;
   int  bdcnt  = 0;
   int  lrcnt  = 0;
   int  ojcnt  = 0;
  
   for(int i = 0; i < lp.nRows(); ++i)
   {
      // Non-zeros
      SVector& row = lp.rowVector_w(i);
      Real     x;
      int      j = 0;

      while(j < row.size())
      {
         x = row.value(j);

         if (isZero(x))
         {
            row.remove(j);
            nzcnt++;
         }         
         else
         {
            if (fabs(x) > maxval)
               std::cerr << "Warning! Big value " << x << std::endl;

            j++;
         }
      }
      // LHS 
      x = lp.lhs(i);

      if (x != 0 && isZero(x))
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

      if (x != 0 && isZero(x))
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
      Real x = lp.lower(i);

      if (x != 0 && isZero(x))
      {
         lp.changeLower(i, 0.0);
         bdcnt++;
      }
      else if (x > -infinity && x < -maxval)
      {
         lp.changeLower(i, -infinity);
         bdcnt++;
      }
      else if (x <  infinity && x >  maxval)
      {
         lp.changeLower(i,  infinity);
         bdcnt++;
      }

      // upper bound
      x = lp.lower(i);

      if (x != 0 && isZero(x))
      {
         lp.changeUpper(i, 0.0);
         bdcnt++;
      }
      else if (x > -infinity && x < -maxval)
      {
         lp.changeUpper(i, -infinity);
         bdcnt++;
      }
      else if (x <  infinity && x >  maxval)
      {
         lp.changeUpper(i,  infinity);
         bdcnt++;
      }

      // objective
      x = lp.obj(i);

      if (x != 0 && isZero(x))
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
