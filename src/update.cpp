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
#pragma ident "@(#) $Id: update.cpp,v 1.3 2001/12/01 18:21:16 bzfbleya Exp $"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "clutypes.h"
#include "clumembers.h"
#include "cluprotos.h"
#include "cring.h"

namespace soplex
{


/*****************************************************************************/

int updateCLUFactor
(
   CLUFactor* fac,
   int col,
   double* work,
   const int* idx,
   int num
)
{
   int ll, i, j;
   int* lidx;
   double* lval;
   double x, div;
   double maxabs;

   assert(work[col] != 0);
   div = 1 / work[col];
   work[col] = 0;

   ll = fac->makeLvec(num, col);
   //   ll = fac->makeLvec(num, col);
   lval = fac->l.val;
   lidx = fac->l.idx;
   maxabs = fac->maxabs;

   for (i = num - 1; (j = idx[i]) != col; --i)
   {
      lidx[ll] = j;
      lval[ll] = div * work[j];
      work[j] = 0;
      ++ll;
   }

   lidx[ll] = col;
   lval[ll] = 1 - div;
   ++ll;

   for (--i; i >= 0; --i)
   {
      j = idx[i];
      lidx[ll] = j;
      lval[ll] = x = div * work[j];
      work[j] = 0;
      ++ll;
      if (x > maxabs)
         maxabs = x;
      else if (-x > maxabs)
         maxabs = -x;
   }

   fac->maxabs = maxabs;
   return CLU_OK;
}

int updateCLUFactorNoClear
(
   CLUFactor* fac,
   int col,
   const double* work,
   const int* idx,
   int num
)
{
   int ll, i, j;
   int* lidx;
   double* lval;
   double x, div;
   double maxabs;

   assert(work[col] != 0);
   div = 1 / work[col];
   ll = fac->makeLvec(num, col);
   //ll = fac->makeLvec(num, col);
   lval = fac->l.val;
   lidx = fac->l.idx;
   maxabs = fac->maxabs;

   for (i = num - 1; (j = idx[i]) != col; --i)
   {
      lidx[ll] = j;
      lval[ll] = div * work[j];
      ++ll;
   }

   lidx[ll] = col;
   lval[ll] = 1 - div;
   ++ll;

   for (--i; i >= 0; --i)
   {
      j = idx[i];
      lidx[ll] = j;
      lval[ll] = x = div * work[j];
      ++ll;
      if (x > maxabs)
         maxabs = x;
      else if (-x > maxabs)
         maxabs = -x;
   }

   fac->maxabs = maxabs;
   return CLU_OK;
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
