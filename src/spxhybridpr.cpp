/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include "spxdefines.h"
#include "spxhybridpr.h"
#include "spxout.h"

namespace soplex
{

bool SPxHybridPR::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS
   if (thesolver != 0 &&
      (thesolver != steep.solver() ||
         thesolver != devex.solver() ||
         thesolver != parmult.solver()))
      return MSGinconsistent("SPxHybridPR");

   return steep.isConsistent()
          && devex.isConsistent()
          && parmult.isConsistent();
#else
   return true;
#endif
}

void SPxHybridPR::load(SPxSolver* p_solver)
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

void SPxHybridPR::setEpsilon(Real eps)
{
   steep.setEpsilon(eps);
   devex.setEpsilon(eps);
   parmult.setEpsilon(eps);
}

void SPxHybridPR::setType(SPxSolver<R>::Type tp)
{
   if (tp == SPxSolver<R>::LEAVE)
   {
      thepricer = &steep;
      this->thesolver->setPricing(SPxSolver<R>::FULL);
   }
   else
   {
      if (this->thesolver->dim() > hybridFactor * this->thesolver->coDim())
      {
         /**@todo I changed from devex to steepest edge pricing here 
          *       because of numerical difficulties, this should be 
          *       investigated.
          */
         // thepricer = &devex;
         thepricer = &steep;
         this->thesolver->setPricing(SPxSolver<R>::FULL);
      }
      else
      {
         thepricer = &parmult;
         this->thesolver->setPricing(SPxSolver<R>::PARTIAL);
      }
   }
   
   MSG_INFO1( (*this->thesolver->spxout), (*this->thesolver->spxout) << "IPRHYB01 switching to "
                        << thepricer->getName() << std::endl; )

   thepricer->setType(tp);
}

void SPxHybridPR::setRep(SPxSolver<R>::Representation rep)
{
   steep.setRep(rep);
   devex.setRep(rep);
   parmult.setRep(rep);
}

int SPxHybridPR::selectLeave()
{
   return thepricer->selectLeave();
}

void SPxHybridPR::left4(int n, SPxId id)
{
   thepricer->left4(n, id);
}

SPxId SPxHybridPR::selectEnter()
{
   return thepricer->selectEnter();
}

void SPxHybridPR::entered4(SPxId id, int n)
{
   thepricer->entered4(id, n);
}

void SPxHybridPR::addedVecs (int n)
{
   steep.addedVecs(n);
   devex.addedVecs(n);
   parmult.addedVecs(n);
}

void SPxHybridPR::addedCoVecs(int n)
{
   steep.addedCoVecs(n);
   devex.addedCoVecs(n);
   parmult.addedCoVecs(n);
}

} // namespace soplex
