/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <assert.h>

#include "soplex2.h"
#include "statistics.h"

namespace soplex
{
   /// solves real LP with recovery mechanism
   void SoPlex2::_solveRealStable(bool acceptOnlyOptimal)
   {
      ///@todo implement starter
      ///@todo implement auto pricing
      ///@todo implement recovery mechanism
      Real feastol = realParam(SoPlex2::FPFEASTOL);
      Real opttol = realParam(SoPlex2::FPOPTTOL);
      bool solvedFromScratch = false;
      bool inFirstRound = true;
      bool solved = false;

      while( !solved && !_isSolveStopped() && feastol < 0.5 && opttol < 0.5 )
      {
         try
         {
            _solver.setFeastol(feastol);
            _solver.setOpttol(opttol);
            _solver.solve();

            // has LP been solved?
            solved = (_solver.status() == SPxSolver::OPTIMAL || !acceptOnlyOptimal)
               && (_solver.status() == SPxSolver::OPTIMAL || _solver.status() == SPxSolver::INFEASIBLE || _solver.status() == SPxSolver::UNBOUNDED);
         }
         catch( SPxException E )
         {
            solved = false;
         }

         // record statistics
         _statistics->iterations += _solver.iterations();

         // if not, relax tolerances or try from scratch
         if( !solved )
         {
            _slufactor.setMarkowitz(0.99);

            if( _solver.status() == SPxSolver::INFEASIBLE )
            {
               _solver.setType(SPxSolver::ENTER);
               feastol *= 10.0;
            }
            else if( _solver.status() == SPxSolver::UNBOUNDED )
            {
               _solver.setType(SPxSolver::LEAVE);
               opttol *= 10.0;
            }

            if( !inFirstRound && !solvedFromScratch )
            {
               _solver.reLoad();
               solvedFromScratch = true;
            }
         }

         inFirstRound = false;
      }
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
