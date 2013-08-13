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
   SPxSolver::Status SoPlex2::_solveRealStable(bool acceptUnbounded, bool acceptInfeasible, VectorReal& primal, VectorReal& dual, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols)
   {
      SPxSolver::Status result = SPxSolver::UNKNOWN;

      bool solved = false;
      bool solvedFromScratch = false;
      bool initialSolve = true;
      bool increasedMarkowitz = false;
      bool relaxedTolerances = false;
      bool tightenedTolerances = false;
      bool switchedRatiotester = false;
      bool switchedPricer = false;

      int ratiotester = intParam(SoPlex2::RATIOTESTER);
      int pricer = intParam(SoPlex2::PRICER);

      while( !_isSolveStopped() )
      {
         assert(!increasedMarkowitz || GE(_slufactor.markowitz(), 0.9));

         try
         {
            MSG_INFO1( spxout << std::endl );

            _solver.solve();
            result = _solver.status();

            MSG_INFO1( spxout << std::endl );

            solved = (result == SPxSolver::OPTIMAL)
               || (result == SPxSolver::INFEASIBLE && acceptInfeasible)
               || (result == SPxSolver::UNBOUNDED && acceptUnbounded);

            if( result == SPxSolver::OPTIMAL )
            {
               assert(basisStatusRows.size() == _solver.nRows());
               assert(basisStatusCols.size() == _solver.nCols());

               _solver.getPrimal(primal);
               _solver.getDual(dual);
               _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());
            }
         }
         catch( ... )
         {
            result = SPxSolver::ERROR;
            solved = false;
         }

         // record statistics
         _statistics->iterations += _solver.iterations();

         if( solved )
            break;

         if( initialSolve )
         {
            MSG_INFO1( spxout << "Numerical troubles during floating-point solve." << std::endl );
            initialSolve = false;
         }

         if( !increasedMarkowitz )
         {
            MSG_INFO1( spxout << "Increasing Markowitz threshold." << std::endl );

            _slufactor.setMarkowitz(0.9);
            increasedMarkowitz = true;
            try
            {
               _solver.factorize();
               continue;
            }
            catch( ... )
            {
               MSG_DEBUG( spxout << std::endl << "Factorization failed." << std::endl );
            }
         }

         if( !solvedFromScratch )
         {
            MSG_INFO1( spxout << "Solving from scratch." << std::endl );

            _solver.reLoad();
            solvedFromScratch = true;
            continue;
         }

         setIntParam(SoPlex2::RATIOTESTER, ratiotester);
         setIntParam(SoPlex2::PRICER, pricer);

         if( !relaxedTolerances && result == SPxSolver::INFEASIBLE )
         {
            MSG_INFO1( spxout << "Relaxing tolerances." << std::endl );

            _solver.setType(_solver.rep() == SPxSolver::COLUMN ? SPxSolver::ENTER : SPxSolver::LEAVE);
            _solver.setDelta(_solver.feastol() * 1e3 > 1e-3 ? 1e-3 : _solver.feastol() * 1e3);
            relaxedTolerances = _solver.feastol() >= 1e-3;
            solvedFromScratch = false;
            continue;
         }

         if( !tightenedTolerances && result != SPxSolver::INFEASIBLE )
         {
            MSG_INFO1( spxout << "Tightening tolerances." << std::endl );

            _solver.setType(_solver.rep() == SPxSolver::COLUMN ? SPxSolver::LEAVE : SPxSolver::ENTER);
            _solver.setDelta(_solver.feastol() * 1e-3 < 1e-9 ? 1e-9 : _solver.feastol() * 1e-3);
            tightenedTolerances = _solver.feastol() <= 1e-9;
            solvedFromScratch = false;
            continue;
         }

         if( !switchedRatiotester )
         {
            MSG_INFO1( spxout << "Switching ratio test." << std::endl );

            _solver.setType(_solver.type() == SPxSolver::LEAVE ? SPxSolver::ENTER : SPxSolver::LEAVE);
            _solver.setTester(_solver.ratiotester() != (SPxRatioTester*)&_ratiotesterTextbook ? (SPxRatioTester*)&_ratiotesterTextbook : (SPxRatioTester*)&_ratiotesterFast);
            switchedRatiotester = true;
            solvedFromScratch = false;
            continue;
         }

         if( !switchedPricer )
         {
            MSG_INFO1( spxout << "Switching pricer." << std::endl );

            _solver.setType(_solver.type() == SPxSolver::LEAVE ? SPxSolver::ENTER : SPxSolver::LEAVE);
            _solver.setPricer(_solver.pricer() != (SPxPricer*)&_pricerDevex ? (SPxPricer*)&_pricerDevex : (SPxPricer*)&_pricerSteep);
            switchedPricer = true;
            solvedFromScratch = false;
            continue;
         }

         MSG_INFO1( spxout << "Giving up." << std::endl );

         break;
      }

      _solver.setFeastol(realParam(SoPlex2::FPFEASTOL));
      _solver.setOpttol(realParam(SoPlex2::FPOPTTOL));

      return result;
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
