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
   /// solves real LP
   SPxSolver::Status SoPlex2::_solveReal(bool fromscratch, VectorReal& primal, VectorReal& dual, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols)
   {
      assert(_isConsistent());

      assert(_solver.nRows() == numRowsRational());
      assert(_solver.nCols() == numColsRational());
      assert(primal.dim() == numColsRational());
      assert(dual.dim() == numRowsRational());

      SPxSolver::Status result = SPxSolver::UNKNOWN;

      if( fromscratch || !_hasBasisRational )
      {
         _enableSimplifierAndScalers();
         _solver.setTerminationValue(realParam(SoPlex2::INFTY));
      }
      else
      {
         _disableSimplifierAndScalers();
         ///@todo implement for both objective senses
         _solver.setTerminationValue(intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE
            ? realParam(SoPlex2::OBJLIMIT_UPPER) : realParam(SoPlex2::INFTY));
      }

      // start timing
      _statistics->syncTime.start();

      SPxLPRational rationalLP(_solver);

      // stop timing
      _statistics->syncTime.stop();

      try
      {
         // apply scaling before the simplification
         if( _firstScaler != 0 )
         {
            _firstScaler->scale(_solver);
         }

         // apply problem simplification
         SPxSimplifier::Result simplificationStatus = SPxSimplifier::OKAY;
         if( _simplifier != 0 )
         {
            simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlex2::EPSILON_ZERO), realParam(SoPlex2::FPFEASTOL), realParam(SoPlex2::FPOPTTOL));
         }

         // apply scaling after the simplification
         if( _secondScaler != 0 && simplificationStatus == SPxSimplifier::OKAY )
         {
            _secondScaler->scale(_solver);
         }

         // run the simplex method if problem has not been solved by the simplifier
         if( simplificationStatus == SPxSimplifier::OKAY )
         {
            MSG_INFO1( spxout << std::endl );

            _solver.solve();

            MSG_INFO1( spxout << std::endl );
         }

         ///@todo move to private helper methods
         // evaluate status flag
         if( simplificationStatus == SPxSimplifier::INFEASIBLE )
            result = SPxSolver::INFEASIBLE;
         else if( simplificationStatus == SPxSimplifier::DUAL_INFEASIBLE )
            result = SPxSolver::INForUNBD;
         else if( simplificationStatus == SPxSimplifier::UNBOUNDED )
            result = SPxSolver::UNBOUNDED;
         else if( simplificationStatus == SPxSimplifier::VANISHED || simplificationStatus == SPxSimplifier::OKAY )
         {
            result = simplificationStatus == SPxSimplifier::VANISHED ? SPxSolver::OPTIMAL : _solver.status();

            ///@todo move to private helper methods
            // process result
            switch( result )
            {
            case SPxSolver::OPTIMAL:
               // unsimplify if simplifier is active and LP is solved to optimality; this must be done here and not at solution
               // query, because we want to have the basis for the original problem
               if( _simplifier != 0 )
               {
                  assert(!_simplifier->isUnsimplified());
                  assert(simplificationStatus == SPxSimplifier::VANISHED || simplificationStatus == SPxSimplifier::OKAY);

                  bool vanished = simplificationStatus == SPxSimplifier::VANISHED;

                  // get solution vectors for transformed problem
                  DVectorReal tmpPrimal(vanished ? 0 : _solver.nCols());
                  DVectorReal tmpSlacks(vanished ? 0 : _solver.nRows());
                  DVectorReal tmpDual(vanished ? 0 : _solver.nRows());
                  DVectorReal tmpRedcost(vanished ? 0 : _solver.nCols());

                  if( !vanished )
                  {
                     assert(_solver.status() == SPxSolver::OPTIMAL);

                     _solver.getPrimal(tmpPrimal);
                     _solver.getSlacks(tmpSlacks);
                     _solver.getDual(tmpDual);
                     _solver.getRedCost(tmpRedcost);

                     // unscale vectors w.r.t. second scaler
                     if( _secondScaler != 0 )
                     {
                        _secondScaler->unscalePrimal(tmpPrimal);
                        _secondScaler->unscaleSlacks(tmpSlacks);
                        _secondScaler->unscaleDual(tmpDual);
                        _secondScaler->unscaleRedCost(tmpRedcost);
                     }

                     // get basis of transformed problem
                     _basisStatusRowsReal.reSize(_solver.nRows());
                     _basisStatusColsReal.reSize(_solver.nCols());
                     _solver.getBasis(_basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());
                  }

                  ///@todo catch exception
                  _simplifier->unsimplify(tmpPrimal, tmpDual, tmpSlacks, tmpRedcost, _basisStatusRowsReal.get_ptr(), _basisStatusColsReal.get_ptr());

                  // store basis for original problem
                  _simplifier->getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());

                  primal = _simplifier->unsimplifiedPrimal();
                  dual = _simplifier->unsimplifiedDual();

                  // unscale vectors w.r.t. first scaler
                  if( _firstScaler != 0 )
                  {
                     _firstScaler->unscalePrimal(primal);
                     _firstScaler->unscaleDual(dual);
                  }
               }
               // if the original problem is not in the solver because of scaling, we also need to store the basis
               else
               {
                  _solver.getPrimal(primal);
                  _solver.getDual(dual);

                  // unscale vectors w.r.t. second scaler
                  if( _secondScaler != 0 )
                  {
                     _secondScaler->unscalePrimal(primal);
                     _secondScaler->unscaleDual(dual);
                  }

                  // get basis of transformed problem
                  _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());

                  // unscale vectors w.r.t. first scaler
                  if( _firstScaler != 0 )
                  {
                     _firstScaler->unscalePrimal(primal);
                     _firstScaler->unscaleDual(dual);
                  }
               }
               break;

            case SPxSolver::ABORT_CYCLING:
            case SPxSolver::ABORT_TIME:
            case SPxSolver::ABORT_ITER:
            case SPxSolver::ABORT_VALUE:
            case SPxSolver::REGULAR:
            case SPxSolver::RUNNING:
            case SPxSolver::UNBOUNDED:
            case SPxSolver::INFEASIBLE:
               // if simplifier is active we cannot return a Farkas ray currently
               if( _simplifier != 0 )
                  break;

               // return Farkas ray as dual solution
               _solver.getDualfarkas(dual);

               // unscale vectors w.r.t. second scaler
               if( _secondScaler != 0 )
                  _secondScaler->unscaleDual(dual);

               // unscale vectors w.r.t. first scaler
               if( _firstScaler != 0 )
                  _firstScaler->unscaleDual(dual);

               // if the original problem is not in the solver because of scaling, we also need to store the basis
               _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());

               break;

            case SPxSolver::INForUNBD:
            case SPxSolver::SINGULAR:
            default:
               _hasBasisReal = false;
               break;
            }
         }
      }
      catch( ... )
      {
         MSG_INFO1( spxout << "Exception thrown during floating-point solve.\n" );
         result = SPxSolver::ERROR;
      }

      // record statistics
      _statistics->iterations += _solver.iterations();

      if( result == SPxSolver::OPTIMAL )
      {
         // check violation in rational LP
         DVectorRational primalRational(primal);
         DVectorRational activity = rationalLP.computePrimalActivity(primalRational);
         Rational maxBoundViolation = 0;
         Rational maxConstraintViolation = 0;

         for( int i = rationalLP.nCols() - 1; i >= 0; i-- )
         {
            Rational viol = rationalLP.lower(i) - primalRational[i];
            if( viol > maxBoundViolation )
               maxBoundViolation = viol;

            viol = primalRational[i] - rationalLP.upper(i);
            if( viol > maxBoundViolation )
               maxBoundViolation = viol;
         }

         for( int i = rationalLP.nRows() - 1; i >= 0; i-- )
         {
            Rational viol = rationalLP.lhs(i) - activity[i];
            if( viol > maxConstraintViolation )
               maxConstraintViolation = viol;

            viol = activity[i] - rationalLP.rhs(i);
            if( viol > maxConstraintViolation )
               maxConstraintViolation = viol;
         }

         Rational violation = (maxBoundViolation > maxConstraintViolation ? maxBoundViolation : maxConstraintViolation);

         if( violation > double(realParam(SoPlex2::FPFEASTOL)) )
         {
            MSG_INFO1( spxout << "Warning: Floating-point solution violates bounds and rows by up to " << rationalToString(violation) << ".\n" );
         }
      }

      // copy rounded rational LP to real LP
      if( _simplifier != 0 || _firstScaler != 0 || _secondScaler != 0 )
         _solver.loadLP((SPxLPReal)(rationalLP));

      return result;
   }



   /// solves real LP with recovery mechanism
   SPxSolver::Status SoPlex2::_solveRealStable(bool acceptUnbounded, bool acceptInfeasible, VectorReal& primal, VectorReal& dual, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols)
   {
      SPxSolver::Status result = SPxSolver::UNKNOWN;

      bool fromScratch = false;
      bool solved = false;
      bool solvedFromScratch = false;
      bool initialSolve = true;
      bool increasedMarkowitz = false;
      bool relaxedTolerances = false;
      bool tightenedTolerances = false;
      bool switchedScalers = false;
      bool switchedSimplifier = false;
      bool switchedRatiotester = false;
      bool switchedPricer = false;

      int ratiotester = intParam(SoPlex2::RATIOTESTER);
      int pricer = intParam(SoPlex2::PRICER);
      int simplifier = intParam(SoPlex2::SIMPLIFIER);
      int scaler_before_simplifier = intParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER);
      int scaler_after_simplifier = intParam(SoPlex2::SCALER_AFTER_SIMPLIFIER);

      setIntParam(SoPlex2::SIMPLIFIER, SoPlex2::SIMPLIFIER_OFF);

      while( !_isSolveStopped() )
      {
         assert(!increasedMarkowitz || GE(_slufactor.markowitz(), 0.9));

         result = _solveReal(fromScratch, primal, dual, basisStatusRows, basisStatusCols);

         solved = (result == SPxSolver::OPTIMAL)
            || (result == SPxSolver::INFEASIBLE && acceptInfeasible)
            || (result == SPxSolver::UNBOUNDED && acceptUnbounded);

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

            fromScratch = true;
            _solver.reLoad();

            solvedFromScratch = true;
            continue;
         }

         setIntParam(SoPlex2::RATIOTESTER, ratiotester);
         setIntParam(SoPlex2::PRICER, pricer);

         if( !switchedScalers )
         {
            MSG_INFO1( spxout << "Switching scaling." << std::endl );

            if( scaler_before_simplifier == int(SoPlex2::SCALER_OFF) && scaler_after_simplifier == int(SoPlex2::SCALER_OFF) )
            {
               setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_BIEQUI);
               setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_BIEQUI);
            }
            else
            {
               setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, SoPlex2::SCALER_OFF);
               setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, SoPlex2::SCALER_OFF);
            }

            fromScratch = true;
            _solver.reLoad();

            solvedFromScratch = true;
            switchedScalers = true;
            continue;
         }

         if( !switchedSimplifier )
         {
            MSG_INFO1( spxout << "Switching simplification." << std::endl );

            if( simplifier == int(SoPlex2::SIMPLIFIER_OFF) )
               setIntParam(SoPlex2::SIMPLIFIER, SoPlex2::SIMPLIFIER_AUTO);
            else
               setIntParam(SoPlex2::SIMPLIFIER, SoPlex2::SIMPLIFIER_OFF);

            fromScratch = true;
            _solver.reLoad();

            solvedFromScratch = true;
            switchedSimplifier = true;
            continue;
         }

         setIntParam(SoPlex2::SIMPLIFIER, SoPlex2::SIMPLIFIER_OFF);

         if( !relaxedTolerances )
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

      setIntParam(SoPlex2::RATIOTESTER, ratiotester);
      setIntParam(SoPlex2::PRICER, pricer);
      setIntParam(SoPlex2::SIMPLIFIER, simplifier);
      setIntParam(SoPlex2::SCALER_BEFORE_SIMPLIFIER, scaler_before_simplifier);
      setIntParam(SoPlex2::SCALER_AFTER_SIMPLIFIER, scaler_after_simplifier);

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
