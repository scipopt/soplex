/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

#include "soplex/spxdefines.h"
#include "soplex/spxhybridpr.h"
#include "soplex/spxout.h"

namespace soplex
{
template <>
void SPxHybridPR<Real>::setType(typename SPxSolverBase<Real>::Type tp);

template <>
bool SPxHybridPR<Real>::isConsistent() const
{
#ifdef ENABLE_CONSISTENCY_CHECKS

   if(thesolver != 0 &&
         (thesolver != steep.this->solver() ||
          thesolver != devex.this->solver() ||
          thesolver != parmult.this->solver()))
      return MSGinconsistent("SPxHybridPR");

   return steep.isConsistent()
          && devex.isConsistent()
          && parmult.isConsistent();
#else
   return true;
#endif
}

template <>
void SPxHybridPR<Real>::load(SPxSolverBase<Real>* p_solver)
{
   steep.load(p_solver);
   devex.load(p_solver);
   parmult.load(p_solver);
   this->thesolver = p_solver;
   setType(p_solver->type());
}

template <>
void SPxHybridPR<Real>::clear()
{
   steep.clear();
   devex.clear();
   parmult.clear();
   this->thesolver = 0;
}

template <>
void SPxHybridPR<Real>::setEpsilon(Real eps)
{
   steep.setEpsilon(eps);
   devex.setEpsilon(eps);
   parmult.setEpsilon(eps);
}

template <>
void SPxHybridPR<Real>::setType(typename SPxSolverBase<Real>::Type tp)
{
   if(tp == SPxSolverBase<Real>::LEAVE)
   {
      thepricer = &steep;
      this->thesolver->setPricing(SPxSolverBase<Real>::FULL);
   }
   else
   {
      if(this->thesolver->dim() > hybridFactor * this->thesolver->coDim())
      {
         /**@todo I changed from devex to steepest edge pricing here
          *       because of numerical difficulties, this should be
          *       investigated.
          */
         // thepricer = &devex;
         thepricer = &steep;
         this->thesolver->setPricing(SPxSolverBase<Real>::FULL);
      }
      else
      {
         thepricer = &parmult;
         this->thesolver->setPricing(SPxSolverBase<Real>::PARTIAL);
      }
   }

   MSG_INFO1((*this->thesolver->spxout), (*this->thesolver->spxout) << "IPRHYB01 switching to "
             << thepricer->getName() << std::endl;)

   thepricer->setType(tp);
}

template <>
void SPxHybridPR<Real>::setRep(typename SPxSolverBase<Real>::Representation rep)
{
   steep.setRep(rep);
   devex.setRep(rep);
   parmult.setRep(rep);
}

template <>
int SPxHybridPR<Real>::selectLeave()
{
   return thepricer->selectLeave();
}

template <>
void SPxHybridPR<Real>::left4(int n, SPxId id)
{
   thepricer->left4(n, id);
}

template <>
SPxId SPxHybridPR<Real>::selectEnter()
{
   return thepricer->selectEnter();
}

template <>
void SPxHybridPR<Real>::entered4(SPxId id, int n)
{
   thepricer->entered4(id, n);
}

template <>
void SPxHybridPR<Real>::addedVecs(int n)
{
   steep.addedVecs(n);
   devex.addedVecs(n);
   parmult.addedVecs(n);
}

template <>
void SPxHybridPR<Real>::addedCoVecs(int n)
{
   steep.addedCoVecs(n);
   devex.addedCoVecs(n);
   parmult.addedCoVecs(n);
}

} // namespace soplex
