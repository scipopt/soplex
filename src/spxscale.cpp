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
#pragma ident "@(#) $Id: spxscale.cpp,v 1.12 2002/03/03 13:50:34 bzfkocht Exp $"

#include <assert.h>
#include <iostream>

#include "spxdefines.h"
#include "spxscale.h"

namespace soplex
{
int SPxScale::simplify()
{
   assert(lp != 0);
   assert(lp->isConsistent());

   Real w;
   Real x;
   Real y;
   Real z;
   int    i; 
   int    j;

   rowscale.reSize(lp->nRows());
   colscale.reSize(lp->nCols());

   if (rowScale)
   {
      for (i = lp->nCols(); i--;)
      {
         SVector& vec = lp->colVector_w(i);

         x = vec.maxAbs();

         if (x > 0)
         {
            y = 1.0 / x;
            colscale[i] = y;
            vec *= y;
            lp->maxObj(i) *= y;
            if (lp->upper(i) < infinity)
               lp->upper(i) *= x;
            if (lp->lower(i) > -infinity)
               lp->lower(i) *= x;
         }
         else
            colscale[i] = 1;
      }

      for (i = lp->nRows(); i--;)
      {
         SVector& vec = lp->rowVector_w(i);
         x = 0;
         z = 1e100;

         for (j = vec.size(); j--;)
         {
            y = vec.value(j) *= colscale[vec.index(j)];
            w = fabs(y);            
            x = (x < w) ? w : x;
            z = (z > w) ? w : z;
         }
         if (x > 0)
         {
            y = 1.0 / x;
            rowscale[i] = y;
            vec *= y;
            if (lp->rhs(i) < infinity)
               lp->rhs(i) *= y;
            if (lp->lhs(i) > -infinity)
               lp->lhs(i) *= y;
         }
         else
            rowscale[i] = 1;
      }
      for (i = lp->nCols(); i--;)
      {
         SVector& vec = lp->colVector_w(i);
         for (j = vec.size(); j--;)
            vec.value(j) *= rowscale[vec.index(j)];
      }
   }
   else
   {
      for (i = lp->nRows(); i--;)
      {
         SVector& vec = lp->rowVector_w(i);
         x = 0;
         for (j = vec.size(); j--;)
         {
            y = vec.value(j);
            if (x < y)
               x = y;
            else if (x < -y)
               x = -y;
         }
         if (x > 0)
         {
            y = 1 / x;
            rowscale[i] = y;
            vec *= y;
            if (lp->rhs(i) < infinity)
               lp->rhs(i) *= y;
            if (lp->lhs(i) > -infinity)
               lp->lhs(i) *= y;
         }
         else
            rowscale[i] = 1;
      }

      for (i = lp->nCols(); i--;)
      {
         SVector& vec = lp->colVector_w(i);
         x = 0;
         for (j = vec.size(); j--;)
         {
            vec.value(j) *= rowscale[vec.index(j)];
            y = vec.value(j);
            if (x < y)
               x = y;
            else if (x < -y)
               x = -y;
         }
         if (x > 0)
         {
            y = 1 / x;
            colscale[i] = y;
            vec *= y;
            lp->maxObj(i) *= y;
            if (lp->upper(i) < infinity)
               lp->upper(i) *= x;
            if (lp->lower(i) > -infinity)
               lp->lower(i) *= x;
         }
         else
            colscale[i] = 1;
      }

      for (i = lp->nRows(); i--;)
      {
         SVector& vec = lp->rowVector_w(i);
         for (j = vec.size(); j--;)
            vec.value(j) *= colscale[vec.index(j)];
      }
   }
   assert(lp->isConsistent());

   return 0;
}

/**@todo We are not unscaling the solution variable values.
 *       these will be reported scaled.
 */
void SPxScale::unsimplify()
{
   assert(lp != 0);
   assert(lp->isConsistent());

   int i;
   int j;

   for (i = lp->nRows(); i--;)
   {
      SVector& vec = lp->rowVector_w(i);
      for (j = vec.size(); j--;)
         vec.value(j) /= colscale[vec.index(j)];
      vec *= 1 / rowscale[i];
      if (lp->rhs(i) < infinity)
         lp->rhs(i) *= rowscale[i];
      if (lp->lhs(i) > -infinity)
         lp->lhs(i) *= rowscale[i];
   }
   for (i = lp->nCols(); i--;)
   {
      SVector& vec = lp->colVector_w(i);
      vec *= 1 / colscale[i];
      for (j = vec.size(); j--;)
         vec.value(j) /= rowscale[vec.index(j)];
      lp->maxObj(i) *= colscale[i];
      if (lp->upper(i) < infinity)
         lp->upper(i) *= 1 / colscale[i];
      if (lp->lower(i) > -infinity)
         lp->lower(i) *= 1 / colscale[i];
   }
   assert(lp->isConsistent());
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
