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
#pragma ident "@(#) $Id: spxdefaultpr.cpp,v 1.1 2001/11/06 16:18:32 bzfkocht Exp $"

/*      \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <assert.h>
#include <iostream>

#define EQ_PREF 1000

/*  and class header files
 */
#include "spxdefaultpr.h"

namespace soplex
{


int SPxDefaultPR::selectLeave(double& best, int start, int incr)
{
   int i, n;
   double x;
   //    const double* up  = thesolver->ubBound();
   //    const double* low = thesolver->lbBound();

   assert(thesolver);
   best = -theeps;
   n = -1;
   for (i = thesolver->dim() - start - 1; i >= 0; i -= incr)
   {
      x = thesolver->fTest()[i];
      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (up[i] == low[i]));
         if (x < best)
         {
            n = i;
            best = x;
         }
      }
   }

   return n;
}

int SPxDefaultPR::selectLeave()
{
   double best;
   return selectLeave(best, 0, 1);
}


SoPlex::Id SPxDefaultPR::selectEnter(double& best, int start1, int incr1,
                                     int start2, int incr2)
{
   SoPlex::Id id;
   int i;
   double x;
   // const SPxBasis::Desc&    ds   = thesolver->basis().desc();

   assert(thesolver);
   best = -theeps;

   for (i = thesolver->dim() - start1 - 1; i >= 0; i -= incr2)
   {
      x = thesolver->coTest()[i];
      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (ds.coStatus(i) == SPxBasis::Desc::P_FREE
         //                || ds.coStatus(i) == SPxBasis::Desc::D_FREE));
         if (x < best)
         {
            id = thesolver->coId(i);
            best = x;
         }
      }
   }

   for (i = thesolver->coDim() - start2 - 1; i >= 0; i -= incr2)
   {
      x = thesolver->test()[i];
      if (x < -theeps)
      {
         // x *= EQ_PREF * (1 + (ds.status(i) == SPxBasis::Desc::P_FREE
         //                || ds.status(i) == SPxBasis::Desc::D_FREE));
         if (x < best)
         {
            id = thesolver->id(i);
            best = x;
         }
      }
   }

   return id;
}

SoPlex::Id SPxDefaultPR::selectEnter()
{
   double best;
   return selectEnter(best, 0, 1, 0, 1);
}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
