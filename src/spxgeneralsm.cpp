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
#pragma ident "@(#) $Id: spxgeneralsm.cpp,v 1.2 2001/11/06 23:31:04 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <stdlib.h>
#include <iostream>


/*  and class header files
 */
#include "spxgeneralsm.h"

namespace soplex
{



//@ ----------------------------------------------------------------------------
void SPxGeneralSM::load(SPxLP* lp)
{
   this->lp = lp;
   rem1.load (lp);
   redu.load (lp);
   aggr.load (lp);
   scale.load(lp);
}

void SPxGeneralSM::unload()
{
   rem1.unload ();
   redu.unload ();
   aggr.unload ();
   scale.unload();
}

SPxLP* SPxGeneralSM::loadedLP() const
{
   return rem1.loadedLP();
}

int SPxGeneralSM::simplify()
{
   int i, cnt;
   int rows = lp->nRows();
   int cols = lp->nCols();

   do
   {
      cnt = lp->nRows() + lp->nCols();
      if ((i = rem1.simplify ()) != 0) return i;
      if ((i = aggr.simplify ()) != 0) return i;
      if ((i = redu.simplify ()) != 0) return i;
   }
   while (0.99*cnt > lp->nRows() + lp->nCols());

   if ((i = scale.simplify()) != 0) return i;

   rows -= lp->nRows();
   cols -= lp->nCols();
   std::cerr << "removed " << rows << " rows\n";
   std::cerr << "removed " << cols << " columns\n";
   assert(lp->isConsistent());

   return 0;
}

void SPxGeneralSM::unsimplify()
{
   rem1.unsimplify ();
   aggr.unsimplify ();
   redu.unsimplify ();
   scale.unsimplify();
}

double SPxGeneralSM::value(double x)
{
   return rem1.value(aggr.value(redu.value(x)));
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
