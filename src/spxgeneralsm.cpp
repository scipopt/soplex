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
#pragma ident "@(#) $Id: spxgeneralsm.cpp,v 1.10 2002/01/31 08:19:28 bzfkocht Exp $"

#include <iostream>

#include "real.h"
#include "spxgeneralsm.h"

namespace soplex
{
void SPxGeneralSM::load(SPxLP* p_lp)
{
   lp       = p_lp;
   rem1.load (p_lp);
   redu.load (p_lp);
   aggr.load (p_lp);
   scale.load(p_lp);
}

void SPxGeneralSM::unload()
{
   rem1.unload ();
   redu.unload ();
   aggr.unload ();
   scale.unload();
}

int SPxGeneralSM::simplify()
{
   int i, cnt;
   int rows = lp->nRows();
   int cols = lp->nCols();

   do
   {
      cnt = lp->nRows() + lp->nCols();

      if ((i = rem1.simplify()) != 0) 
         return i;
      if ((i = aggr.simplify()) != 0) 
         return i;
      if ((i = redu.simplify()) != 0) 
         return i;
   }
   while (0.99 * cnt > lp->nRows() + lp->nCols());

   if ((i = scale.simplify()) != 0) 
      return i;

   rows -= lp->nRows();
   cols -= lp->nCols();

   std::cout << "removed " << rows << " rows\n";
   std::cout << "removed " << cols << " columns\n";

   assert(lp->isConsistent());

   return 0;
}

/**@todo This is not correctly implented, since the simplifiers may be
 *      called several times one after the others, this sequence has
 *      to be tracked to make the calls of unsimplify in exactly the
 *      reverse order.
 *      At the moment this is irrelevant because all the unsimplifiers
 *      are not implemented anyway.
 */
void SPxGeneralSM::unsimplify()
{
   //rem1.unsimplify ();
   //aggr.unsimplify ();
   //redu.unsimplify ();
   scale.unsimplify();
}

Real SPxGeneralSM::value(Real x)
{
   return rem1.value(aggr.value(redu.value(scale.value(x))));
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
