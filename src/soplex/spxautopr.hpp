/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2020 Konrad-Zuse-Zentrum                            */
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
#include "soplex/spxautopr.h"
#include "soplex/spxout.h"

namespace soplex
{

template <class R>
void SPxAutoPR<R>::load(SPxSolverBase<R>* p_solver)
{
   steep.load(p_solver);
   devex.load(p_solver);
   this->thesolver = p_solver;
   setType(p_solver->type());
}

template <class R>
void SPxAutoPR<R>::clear()
{
   steep.clear();
   devex.clear();
   this->thesolver = nullptr;
}

template <class R>
void SPxAutoPR<R>::setEpsilon(R eps)
{
   steep.setEpsilon(eps);
   devex.setEpsilon(eps);
   this->theeps = eps;
}

template <class R>
void SPxAutoPR<R>::setType(typename SPxSolverBase<R>::Type tp)
{
   activepricer->setType(tp);
}

template <class R>
void SPxAutoPR<R>::setRep(typename SPxSolverBase<R>::Representation rep)
{
   steep.setRep(rep);
   devex.setRep(rep);
}

template <class R>
bool SPxAutoPR<R>::setActivePricer(typename SPxSolverBase<R>::Type type)
{
   // switch to steep as soon as switchIters is reached
   if(activepricer == &devex && this->thesolver->iterations() >= switchIters)
   {
      activepricer = &steep;
      activepricer->setType(type);
      return true;
   }


   // use devex for the iterations < switchIters
   else if(activepricer == &steep && this->thesolver->iterations() < switchIters)
   {
      activepricer = &devex;
      activepricer->setType(type);
      return true;
   }

   return false;
}

template <class R>
int SPxAutoPR<R>::selectLeave()
{
   if(setActivePricer(SPxSolverBase<R>::LEAVE))
      MSG_INFO1((*this->thesolver->spxout),
                (*this->thesolver->spxout) << " --- active pricer: " << activepricer->getName() << std::endl;)

      return activepricer->selectLeave();
}

template <class R>
void SPxAutoPR<R>::left4(int n, SPxId id)
{
   activepricer->left4(n, id);
}

template <class R>
SPxId SPxAutoPR<R>::selectEnter()
{
   if(setActivePricer(SPxSolverBase<R>::ENTER))
      MSG_INFO1((*this->thesolver->spxout),
                (*this->thesolver->spxout) << " --- active pricer: " << activepricer->getName() << std::endl;)

      return activepricer->selectEnter();
}

template <class R>
void SPxAutoPR<R>::entered4(SPxId id, int n)
{
   activepricer->entered4(id, n);
}

} // namespace soplex
