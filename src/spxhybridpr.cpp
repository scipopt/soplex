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
#pragma ident "@(#) $Id: spxhybridpr.cpp,v 1.3 2001/11/12 16:42:05 bzfpfend Exp $"

/*  \Section{Complex Methods}
 */

/*  Import system include files
 */
#include <stdlib.h>
#include <iostream>


/*  and class header files
 */
#include "spxhybridpr.h"

namespace soplex
{

//@ ----------------------------------------------------------------------------
/*  \Section{Miscellanous Methods}
 */
int SPxHybridPR::isConsistent() const
{
   if (thesolver != 0 &&
        (thesolver != steep.solver() ||
          thesolver != devex.solver() ||
          thesolver != parmult.solver()))
   {
      std::cerr << "ERROR: inconsistency detected in class SPxHybridPR\n";
      return 0;
   }
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

void SPxHybridPR::addedVecs (int n)
{}


void SPxHybridPR::addedCoVecs(int n)
{}


void SPxHybridPR::removedVec(int i)
{}


void SPxHybridPR::removedVecs(const int perm[])
{}


void SPxHybridPR::removedCoVec(int i)
{}


void SPxHybridPR::removedCoVecs(const int perm[])
{}


void SPxHybridPR::changeObj(const Vector& newObj)
{}


void SPxHybridPR::changeObj(int i, double newVal)
{}


void SPxHybridPR::changeLower(const Vector& newLower)
{}


void SPxHybridPR::changeLower(int i, double newLower)
{}


void SPxHybridPR::changeUpper(const Vector& newUpper)
{}


void SPxHybridPR::changeUpper(int i, double newUpper)
{}


void SPxHybridPR::changeLhs(const Vector& newLhs)
{}


void SPxHybridPR::changeLhs(int i, double newLhs)
{}


void SPxHybridPR::changeRhs(const Vector& newRhs)
{}


void SPxHybridPR::changeRhs(int i, double newRhs)
{}


void SPxHybridPR::changeRow(int i, const LPRow& newRow)
{}


void SPxHybridPR::changeCol(int i, const LPCol& newCol)
{}


void SPxHybridPR::changeElement(int i, int j, double val)
{}


void SPxHybridPR::changeSense(SoPlex::Sense sns)
{}
} // namespace soplex

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
