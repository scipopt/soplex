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
#pragma ident "@(#) $Id: spxhybridpr.cpp,v 1.9 2001/12/26 12:58:58 bzfkocht Exp $"

#include <iostream>

#include "spxhybridpr.h"
#include "message.h"

namespace soplex
{

int SPxHybridPR::isConsistent() const
{
   if (thesolver != 0 &&
      (thesolver != steep.solver() ||
         thesolver != devex.solver() ||
         thesolver != parmult.solver()))
      return MSGinconsistent("SPxHybridPR");

   return steep.isConsistent()
          && devex.isConsistent()
          && parmult.isConsistent();
}

void SPxHybridPR::load(SoPlex* p_solver)
{
   steep.load(p_solver);
   devex.load(p_solver);
   parmult.load(p_solver);
   thesolver = p_solver;
   setType(p_solver->type());
}

void SPxHybridPR::clear()
{
   steep.clear();
   devex.clear();
   parmult.clear();
   thesolver = 0;
}

void SPxHybridPR::setEpsilon(double eps)
{
   steep.setEpsilon(eps);
   devex.setEpsilon(eps);
   parmult.setEpsilon(eps);
}

void SPxHybridPR::setType(SoPlex::Type tp)
{
   if (tp == SoPlex::LEAVE)
   {
      thepricer = &steep;
      thesolver->setPricing(SoPlex::FULL);
   }
   else
   {
      if (thesolver->dim() > hybridFactor*thesolver->coDim())
      {
         thepricer = &devex;
         thesolver->setPricing(SoPlex::FULL);
#ifndef NDEBUG
         std::cerr << "switching to devex\n";
#endif
      }
      else
      {
         thepricer = &parmult;
         thesolver->setPricing(SoPlex::PARTIAL);
#ifndef NDEBUG
         std::cerr << "switching to partial multiple pricing\n";
#endif
      }
   }
   thepricer->setType(tp);
}

void SPxHybridPR::setRep(SoPlex::Representation rep)
{
   steep.setRep(rep);
   devex.setRep(rep);
   parmult.setRep(rep);
}

int SPxHybridPR::selectLeave()
{
   return thepricer->selectLeave();
}

void SPxHybridPR::left4(int n, SoPlex::Id id)
{
   thepricer->left4(n, id);
}

SoPlex::Id SPxHybridPR::selectEnter()
{
   return thepricer->selectEnter();
}

void SPxHybridPR::entered4(SoPlex::Id id, int n)
{
   thepricer->entered4(id, n);
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
