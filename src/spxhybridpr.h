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
#pragma ident "@(#) $Id: spxhybridpr.h,v 1.15 2005/09/16 12:42:35 bzfhille Exp $"

/**@file  spxhybridpr.h
 * @brief Hybrid pricer.
 */
#ifndef _SPXHYBRIDPR_H_
#define _SPXHYBRIDPR_H_

#include <assert.h>

#include "spxdefines.h"
#include "spxpricer.h"
#include "spxdevexpr.h"
#include "spxparmultpr.h"
#include "spxsteeppr.h"

namespace soplex
{

/**@brief   Hybrid pricer.
   @ingroup Algo

   The hybrid pricer for SoPlex tries to guess the best pricing strategy to
   use for pricing the loaded LP with the loaded algorithm type and basis
   representation. Currently it does so by switching between #SPxSteepPR,
   #SPxDevexPR and #SPxParMultPR.

   See #SPxPricer for a class documentation.
*/
class SPxHybridPR : public SPxPricer
{
   SPxSteepPR   steep;
   SPxParMultPR parmult;
   SPxDevexPR   devex;

   SPxPricer*   thepricer;

   Real hybridFactor; 

public:
   ///
   virtual void setEpsilon(Real eps);
   ///
   virtual void load(SPxSolver* solver);
   ///
   virtual void clear();
   ///
   virtual void setType(SPxSolver::Type tp);
   ///
   virtual void setRep(SPxSolver::Representation rep);
   ///
   virtual int selectLeave();
   ///
   virtual void left4(int n, SPxId id);
   ///
   virtual SPxId selectEnter();
   ///
   virtual void entered4(SPxId id, int n);
#ifndef NO_CONSISTENCY_CHECKS
   ///
   virtual bool isConsistent() const;
#endif
   ///
   virtual void addedVecs (int n);
   ///
   virtual void addedCoVecs (int n);
   ///
   SPxHybridPR() 
      : SPxPricer("Hybrid")
      , thepricer(0)
      , hybridFactor(3.0) // we want the ParMult pricer
   {}
};

} // namespace soplex
#endif // _SPXHYBRIDPR_H_

//-----------------------------------------------------------------------------
//Emacs Local Variables:
//Emacs mode:c++
//Emacs c-basic-offset:3
//Emacs tab-width:8
//Emacs indent-tabs-mode:nil
//Emacs End:
//-----------------------------------------------------------------------------
