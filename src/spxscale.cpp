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
#pragma ident "@(#) $Id: spxscale.cpp,v 1.7 2001/11/29 22:52:54 bzfkocht Exp $"

#include <assert.h>
#include <iostream>

#include "spxscale.h"

namespace soplex
{
static const double MinScale = 1e-4;

/** This routine looks for a suitabel scaling factor.
 *  First the maximum value is scaled to 1.0.
 *  If this results in the minimum to drop below #min_scale
 *  The new mimium is scaled up to #min_scale.
 *  The scaling factor returned is the combination of both factors.
 */
static double find_scale(double vmin, double vmax)
{
   assert(vmin <= vmax);

   double fmax  = 1.0 / vmax;
   double fmin  = 1.0;

   if (vmin > 1e-16)
   { 
      vmin *= fmax;

      if (vmin < MinScale)
         fmin = MinScale / vmin;
   }
   return fmax * fmin;  
}

int SPxScale::simplify()
{
   assert(lp != 0);
   assert(lp->isConsistent());

   double w;
   double x;
   double y;
   double z;
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
            y = find_scale(vec.minAbs(), x);
            //y = 1 / x;
            colscale[i] = y;
            vec *= y;
            lp->maxObj(i) *= y;
            if (lp->upper(i) < SPxLP::infinity)
               lp->upper(i) *= x;
            if (lp->lower(i) > -SPxLP::infinity)
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
            y = find_scale(z, x);
            //y = 1 / x;
            rowscale[i] = y;
            vec *= y;
            if (lp->rhs(i) < SPxLP::infinity)
               lp->rhs(i) *= y;
            if (lp->lhs(i) > -SPxLP::infinity)
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
            if (lp->rhs(i) < SPxLP::infinity)
               lp->rhs(i) *= y;
            if (lp->lhs(i) > -SPxLP::infinity)
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
            if (lp->upper(i) < SPxLP::infinity)
               lp->upper(i) *= x;
            if (lp->lower(i) > -SPxLP::infinity)
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
      if (lp->rhs(i) < SPxLP::infinity)
         lp->rhs(i) *= rowscale[i];
      if (lp->lhs(i) > -SPxLP::infinity)
         lp->lhs(i) *= rowscale[i];
   }
   for (i = lp->nCols(); i--;)
   {
      SVector& vec = lp->colVector_w(i);
      vec *= 1 / colscale[i];
      for (j = vec.size(); j--;)
         vec.value(j) /= rowscale[vec.index(j)];
      lp->maxObj(i) *= colscale[i];
      if (lp->upper(i) < SPxLP::infinity)
         lp->upper(i) *= 1 / colscale[i];
      if (lp->lower(i) > -SPxLP::infinity)
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
