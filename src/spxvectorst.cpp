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
#pragma ident "@(#) $Id: spxvectorst.cpp,v 1.2 2001/11/06 23:31:06 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>


/*  and class header files
 */
#include "spxvectorst.h"

namespace soplex
{

#define ABS(x)          ((x)<0   ? -(x) : (x))
#define MIN(x,y)        ((x)<(y) ? (x) : (y))

//@ ----------------------------------------------------------------------------
void SPxVectorST::setupWeights(SoPlex& base)
{
   if (state == PVEC)
   {
      if (vec.dim() != base.nCols())
      {
         SPxWeightST::setupWeights(base);
         return;
      }

      const Vector& obj = base.maxObj();
      double eps = base.epsilon();
      double bias = 10000 * eps;
      double x, y;
      int i;

      for (i = base.nCols(); i--;)
      {
         x = vec[i] - base.SPxLP::lower(i);
         y = base.SPxLP::upper(i) - vec[i];
         if (x < y)
         {
            colWeight[i] = -x - bias * obj[i];
            colUp[i] = 0;
         }
         else
         {
            colWeight[i] = -y + bias * obj[i];
            colUp[i] = 1;
         }
         //@ std::cout << colWeight[i] << " ";
      }
      //@ std::cout << std::endl;
      //@ std::cout << std::endl;


      for (i = base.nRows(); i--;)
      {
         const SVector& row = ((const SoPlex&)base).rowVector(i);
         y = vec * row;
         x = (y - base.lhs(i));
         y = (base.rhs(i) - y);
         if (x < y)
         {
            rowWeight[i] = -x - eps * row.size() - bias * (obj * row);
            rowRight[i] = 0;
         }
         else
         {
            rowWeight[i] = -y - eps * row.size() + bias * (obj * row);
            rowRight[i] = 1;
         }
         //@ std::cout << rowWeight[i] << " ";
      }
      //@ std::cout << std::endl;
   }

   else if (state == DVEC)
   {
      if (vec.dim() != base.nRows())
      {
         SPxWeightST::setupWeights(base);
         return;
      }

      double x, y, len;
      int i, j;
      for (i = base.nRows(); i--;)
         rowWeight[i] += ABS(vec[i]);

      for (i = base.nCols(); i--;)
      {
         const SVector& col = ((const SoPlex&)base).colVector(i);
         for (y = len = 0, j = col.size(); j--;)
         {
            x = col.value(j);
            y += vec[col.index(j)] * x;
            len += x * x;
         }
         if (len > 0)
            colWeight[i] += ABS(y / len - base.maxObj(i));
      }
   }
   else
      SPxWeightST::setupWeights(base);
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
