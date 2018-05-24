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

#include "soplex/spxdefines.h"
#include "soplex/spxautopr.h"
#include "soplex/spxout.h"

namespace soplex
{

  template <>
  void SPxAutoPR<Real>::setType(typename SPxSolver<Real>::Type tp);


  template <>
  void SPxAutoPR<Real>::load(SPxSolver<Real>* p_solver)
  {
    steep.load(p_solver);
    devex.load(p_solver);
    this->thesolver = p_solver;
    setType(p_solver->type());
  }

  template <>
  void SPxAutoPR<Real>::clear()
  {
    steep.clear();
    devex.clear();
    this->thesolver = 0;
  }

  template <>
  void SPxAutoPR<Real>::setEpsilon(Real eps)
  {
    steep.setEpsilon(eps);
    devex.setEpsilon(eps);
    this->theeps = eps;
  }

  template <>
  void SPxAutoPR<Real>::setType(typename SPxSolver<Real>::Type tp)
  {
    activepricer->setType(tp);
  }

  template <>
  void SPxAutoPR<Real>::setRep(typename SPxSolver<Real>::Representation rep)
  {
    steep.setRep(rep);
    devex.setRep(rep);
  }

  template <>
  bool SPxAutoPR<Real>::setActivePricer(typename SPxSolver<Real>::Type type)
  {
    // switch to steep as soon as switchIters is reached
    if( activepricer == &devex && this->thesolver->iterations() >= switchIters )
      {
        activepricer = &steep;
        activepricer->setType(type);
        return true;
      }


    // use devex for the iterations < switchIters
    else if( activepricer == &steep && this->thesolver->iterations() < switchIters  )
      {
        activepricer = &devex;
        activepricer->setType(type);
        return true;
      }
  
    return false;
  }

  template <>
  int SPxAutoPR<Real>::selectLeave()
  {
    if( setActivePricer(SPxSolver<Real>::LEAVE) )
      MSG_INFO1( (*this->thesolver->spxout), (*this->thesolver->spxout) << " --- active pricer: " << activepricer->getName() << std::endl; )

        return activepricer->selectLeave();
  }

  template <>
  void SPxAutoPR<Real>::left4(int n, SPxId id)
  {
    activepricer->left4(n, id);
  }

  template <>
  SPxId SPxAutoPR<Real>::selectEnter()
  {
    if( setActivePricer(SPxSolver<Real>::ENTER) )
      MSG_INFO1( (*this->thesolver->spxout), (*this->thesolver->spxout) << " --- active pricer: " << activepricer->getName() << std::endl; )

        return activepricer->selectEnter();
  }

  template <>
  void SPxAutoPR<Real>::entered4(SPxId id, int n)
  {
    activepricer->entered4(id, n);
  }

} // namespace soplex
