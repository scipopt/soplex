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
   /// solves rational LP
   void SoPlex2::_solveRational()
   {
      bool hasUnboundedRay = false;
      bool infeasibilityNotCertified = false;
      bool unboundednessNotCertified = false;

      _statusRational = SPxSolver::UNKNOWN;

      // ensure that the solver has the original problem
      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;
      }

      // deactivate objective limit in floating-point solver
      if( realParam(SoPlex2::OBJLIMIT_LOWER) > -realParam(SoPlex2::INFTY) || realParam(SoPlex2::OBJLIMIT_UPPER) < realParam(SoPlex2::INFTY) )
      {
         MSG_INFO2( spxout << "Deactivating objective limit." << std::endl );
      }

      _solver.setTerminationValue(realParam(SoPlex2::INFTY));

      // introduce slack variables to transform inequality constraints into equations
      _transformEquality();

      do
      {
         bool primalFeasible = false;
         bool dualFeasible = false;
         bool infeasible = false;
         bool unbounded = false;
         bool stopped = false;
         bool error = false;

         // solve problem with iterative refinement and recovery mechanism
         _performOptIRStable(_solRational, !unboundednessNotCertified, !infeasibilityNotCertified,
            primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

         // case: an unrecoverable error occured
         if( error )
         {
            _statusRational = SPxSolver::ERROR;
            break;
         }
         // case: stopped due to some limit
         else if( stopped )
         {
            _statusRational = SPxSolver::ABORT_TIME;
            break;
         }
         // case: unboundedness detected for the first time
         else if( unbounded && !unboundednessNotCertified )
         {
            SolRational solUnbounded;

            _performUnboundedIRStable(solUnbounded, hasUnboundedRay, stopped, error);

            assert(!hasUnboundedRay || solUnbounded.hasPrimalray());
            assert(!solUnbounded.hasPrimalray() || hasUnboundedRay);

            if( error )
            {
               MSG_INFO1( spxout << "Error while testing for unboundedness." << std::endl );
               _statusRational = SPxSolver::ERROR;
               break;
            }

            if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "Dual infeasible.  Primal unbounded ray available." << std::endl );
            }
            else
            {
               MSG_INFO1( spxout << "Dual feasible.  Rejecting primal unboundedness." << std::endl );
            }

            unboundednessNotCertified = !hasUnboundedRay;

            if( stopped )
            {
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }

            _performFeasIRStable(_solRational, infeasible, stopped, error);

            if( hasUnboundedRay )
            {
               _solRational._primalray = solUnbounded._primalray;
               _solRational._hasPrimalray = true;
            }

            if( error )
            {
               MSG_INFO1( spxout << "Error while testing for feasibility." << std::endl );
               _statusRational = SPxSolver::ERROR;
               break;
            }
            else if( stopped )
            {
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "Primal infeasible.  Dual Farkas ray available." << std::endl );
               _statusRational = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "Primal feasible and unbounded." << std::endl );
               _statusRational = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "Primal feasible and bounded." << std::endl );
               continue;
            }
         }
         // case: infeasibility detected
         else if( infeasible && !infeasibilityNotCertified )
         {
            _performFeasIRStable(_solRational, infeasible, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "Error while testing for infeasibility." << std::endl );
               _statusRational = SPxSolver::ERROR;
               break;
            }

            infeasibilityNotCertified = !infeasible;

            if( stopped )
            {
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "Primal infeasible.  Dual Farkas ray available." << std::endl );
               _statusRational = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "Primal feasible and unbounded." << std::endl );
               _statusRational = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "Primal feasible.  Optimizing again." << std::endl );
               continue;
            }
         }
         else if( primalFeasible && dualFeasible )
         {
            MSG_INFO1( spxout << "Solved to optimality." << std::endl );
            _statusRational = SPxSolver::OPTIMAL;
            break;
         }
         else
         {
            MSG_INFO1( spxout << "Terminating without success." << std::endl );
            break;
         }
      }
      while( !_isSolveStopped() );

      if( _isSolveStopped() )
         _statusRational = SPxSolver::ABORT_TIME;

      ///@todo set status to ABORT_VALUE if optimal solution exceeds objective limit

      // restore original problem
      _untransformEquality(_solRational);
   }



   /// solves current problem with iterative refinement and recovery mechanism
   void SoPlex2::_performOptIRStable(SolRational& sol, bool acceptUnbounded, bool acceptInfeasible, bool& primalFeasible, bool& dualFeasible, bool& infeasible, bool& unbounded, bool& stopped, bool& error)
   {
      primalFeasible = false;
      dualFeasible = false;
      infeasible = false;
      unbounded = false;
      stopped = false;
      error = false;

      // set working tolerances in floating-point solver
      _solver.setFeastol(realParam(SoPlex2::FPFEASTOL));
      _solver.setOpttol(realParam(SoPlex2::FPOPTTOL));

      // solve original LP
      MSG_INFO1( spxout << "Initial floating-point solve . . ." << std::endl );

      if( intParam(SoPlex2::ITERLIMIT) >= 0 )
         _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT) - _statistics->iterations);

      if( realParam(SoPlex2::TIMELIMIT) < realParam(SoPlex2::INFTY) )
         _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT) - _statistics->solvingTime.userTime());

      if( _hasBasisRational )
         _solver.setBasis(_basisStatusRowsRational.get_const_ptr(), _basisStatusColsRational.get_const_ptr());

      _solveRealStable(acceptUnbounded, acceptInfeasible);

      _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT));
      _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT));

      // evaluate result
      switch( _solver.status() )
      {
      case SPxSolver::OPTIMAL:
         MSG_INFO1( spxout << "Floating-point optimal." << std::endl );
         break;
      case SPxSolver::INFEASIBLE:
         MSG_INFO1( spxout << "Floating-point infeasible." << std::endl );
         infeasible = true;
         return;
      case SPxSolver::UNBOUNDED:
         MSG_INFO1( spxout << "Floating-point unbounded." << std::endl );
         unbounded = true;
         return;
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
         stopped = true;
         return;
      default:
         error = true;
         return;
      }

      // declare vectors and variables
      DVectorRational modLower(numColsRational());
      DVectorRational modUpper(numColsRational());
      DVectorRational modSide(numRowsRational());
      DVectorRational modObj(numColsRational());
      DVectorReal primalReal(numColsRational());
      DVectorReal dualReal(numRowsRational());

      Rational boundsViolation;
      Rational sideViolation;
      Rational redcostViolation;
      Rational primalScale;
      Rational dualScale;
      Rational maxScale;

      // get floating-point solution of original LP
      _solver.getPrimal(primalReal);
      _solver.getDual(dualReal);

      // store floating-point solution of original LP as current rational solution
      sol._primal = primalReal;
      sol._slacks = _rationalLP->computePrimalActivity(sol._primal);
      sol._hasPrimal = true;

      sol._dual = dualReal;
      sol._redcost = _rationalLP->computeDualActivity(sol._dual) + _rationalLP->maxObj();
      sol._redcost *= -1;
      sol._hasDual = true;

      _basisStatusRowsRational.reSize(numRowsRational());
      _basisStatusColsRational.reSize(numColsRational());
      (void)_solver.getBasis(_basisStatusRowsRational.get_ptr(), _basisStatusColsRational.get_ptr());
      _hasBasisRational = true;

      // initial scaling factors are one
      primalScale = 1;
      dualScale = 1;

      // refinement loop
      do
      {
         assert(_solver.status() == SPxSolver::OPTIMAL);

         MSG_DEBUG( spxout << "Computing violations." << std::endl );

         // compute violation of bounds
         boundsViolation = 0;

         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            // lower bound
            modLower[c] = lowerRational(c);

            if( modLower[c] > -realParam(SoPlex2::INFTY) )
               modLower[c] -= sol._primal[c];

            if( modLower[c] > boundsViolation )
               boundsViolation = modLower[c];

            // upper bound
            modUpper[c] = upperRational(c);

            if( modUpper[c] < realParam(SoPlex2::INFTY) )
               modUpper[c] -= sol._primal[c];

            if( modUpper[c] < -boundsViolation )
               boundsViolation = -modUpper[c];
         }

         // compute violation of sides
         modSide = sol._slacks;
         sideViolation = 0;

         for( int r = numRowsRational() - 1; r >= 0; r-- )
         {
            assert(lhsRational(r) == rhsRational(r));

            modSide[r] *= -1;
            modSide[r] += rhsRational(r);

            if( modSide[r] > sideViolation )
               sideViolation = modSide[r];
            else if( modSide[r] < -sideViolation )
               sideViolation = -modSide[r];
         }

         // compute reduced costs and reduced cost violation
         modObj = _rationalLP->computeDualActivity(sol._dual);
         redcostViolation = 0;

         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            modObj[c] *= -1;
            modObj[c] += objRational(c);

            SPxSolver::VarStatus basisStatus = _basisStatusColsRational[c];

            if( basisStatus != SPxSolver::ON_UPPER && basisStatus != SPxSolver::FIXED && modObj[c] < -redcostViolation )
               redcostViolation = -modObj[c];

            if( basisStatus != SPxSolver::ON_LOWER && basisStatus != SPxSolver::FIXED && modObj[c] > redcostViolation )
               redcostViolation = modObj[c];
         }

         // output violations; the reduced cost violations for artificially introduced slack columns are actually violations of the dual multipliers
         MSG_INFO2( spxout
            << "Max. bound violation = " << rationalToString(boundsViolation) << std::endl
            << "Max. row violation = " << rationalToString(sideViolation) << std::endl
            << "Max. dual violation = " << rationalToString(redcostViolation) << std::endl );

         // terminate if tolerances are satisfied
         primalFeasible = (boundsViolation <= rationalParam(SoPlex2::FEASTOL) && sideViolation <= rationalParam(SoPlex2::FEASTOL));
         dualFeasible = (redcostViolation <= rationalParam(SoPlex2::OPTTOL));
         if( primalFeasible && dualFeasible )
         {
            MSG_INFO1( spxout << "Tolerances reached." << std::endl );
            return;
         }

         // terminate if some limit is reached
         if( _isSolveStopped() )
         {
            stopped = true;
            return;
         }

         ///@todo try rational reconstruction at geometric frequency

         // otherwise start refinement
         _statistics->refinements++;

         // compute primal scaling factor; limit increase in scaling by tolerance used in floating point solve
         maxScale = primalScale * Rational(realParam(SoPlex2::MAXSCALEINCR));

         primalScale = boundsViolation > sideViolation ? boundsViolation : sideViolation;
         assert(primalScale >= 0);

         if( primalScale > 0 )
         {
            primalScale = Rational(1) / primalScale;
            if( primalScale > maxScale )
               primalScale = maxScale;
         }
         else
            primalScale = maxScale;

         if( primalScale < 1 )
            primalScale = 1;

         MSG_INFO2( spxout << "Scaling primal by " << rationalToString(primalScale) << "." << std::endl );

         // compute dual scaling factor; limit increase in scaling by tolerance used in floating point solve
         maxScale = dualScale * Rational(realParam(SoPlex2::MAXSCALEINCR));

         dualScale = redcostViolation;
         assert(dualScale >= 0);

         if( dualScale > 0 )
         {
            dualScale = Rational(1) / dualScale;
            if( dualScale > maxScale )
               dualScale = maxScale;
         }
         else
            dualScale = maxScale;

         if( dualScale < 1 )
            dualScale = 1;

         MSG_INFO2( spxout << "Scaling dual by " << rationalToString(dualScale) << "." << std::endl );

         // perform primal and dual scaling
         modLower *= primalScale;
         modUpper *= primalScale;
         modSide *= primalScale;
         modObj *= dualScale;

         // apply scaled bounds, side, and objective function
         _solver.changeBounds(DVectorReal(modLower), DVectorReal(modUpper));
         _solver.changeRange(DVectorReal(modSide), DVectorReal(modSide));
         _solver.changeObj(DVectorReal(modObj));

         MSG_INFO1( spxout << "Refined floating-point solve . . ." << std::endl );

         // load basis
         _solver.setBasis(_basisStatusRowsRational.get_const_ptr(), _basisStatusColsRational.get_const_ptr());

         // solve modified problem
         if( intParam(SoPlex2::ITERLIMIT) >= 0 )
            _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT) - _statistics->iterations);

         if( realParam(SoPlex2::TIMELIMIT) < realParam(SoPlex2::INFTY) )
            _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT) - _statistics->solvingTime.userTime());

         _solveRealStable(acceptUnbounded, acceptInfeasible);

         _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT));
         _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT));

         // remember whether we moved to a new basis
         if( _solver.iterations() == 0 )
            _statistics->stallRefinements++;

         // evaluate result; if modified problem was not solved to optimality, stop refinement
         switch( _solver.status() )
         {
         case SPxSolver::OPTIMAL:
            MSG_INFO1( spxout << "Floating-point optimal." << std::endl );
            break;
         case SPxSolver::INFEASIBLE:
            MSG_INFO1( spxout << "Floating-point infeasible." << std::endl );
            infeasible = true;
            return;
         case SPxSolver::UNBOUNDED:
            MSG_INFO1( spxout << "Floating-point unbounded." << std::endl );
            unbounded = true;
            return;
         case SPxSolver::ABORT_TIME:
         case SPxSolver::ABORT_ITER:
            stopped = true;
            return;
         default:
            error = true;
            return;
         }

         // get floating point solution of modified problem
         _solver.getPrimal(primalReal);
         _solver.getDual(dualReal);

         // get basis
         (void)_solver.getBasis(_basisStatusRowsRational.get_ptr(), _basisStatusColsRational.get_ptr());
         assert(_hasBasisRational);

         // correct primal solution
         MSG_DEBUG( spxout << "Correcting primal solution." << std::endl );

         int numAdjustedBounds = 0;
         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            sol._primal[c] += Rational(primalReal[c]) / primalScale;

            // force values of nonbasic variables to bounds
            SPxSolver::VarStatus basisStatus = _basisStatusColsRational[c];

            if( basisStatus == SPxSolver::ON_LOWER && sol._primal[c] != lowerRational(c) )
            {
               sol._primal[c] = lowerRational(c);
               numAdjustedBounds++;
            }
            else if( basisStatus == SPxSolver::ON_UPPER && sol._primal[c] != upperRational(c) )
            {
               sol._primal[c] = upperRational(c);
               numAdjustedBounds++;
            }
            else if( basisStatus == SPxSolver::FIXED )
            {
               assert(lowerRational(c) == upperRational(c));

               if( sol._primal[c] != lowerRational(c) )
               {
                  sol._primal[c] = lowerRational(c);
                  numAdjustedBounds++;
               }
            }
            else if( basisStatus == SPxSolver::ZERO && sol._primal[c] != 0 )
            {
               sol._primal[c] = 0;
               numAdjustedBounds++;
            }
         }

         // correct dual solution
         MSG_DEBUG( spxout << "Correcting dual solution." << std::endl );

         for( int r = numRowsRational() - 1; r >= 0; r-- )
         {
            sol._dual[r] += Rational(dualReal[r]) / dualScale;
         }

         // recompute slack and reduced cost values
         sol._slacks = _rationalLP->computePrimalActivity(sol._primal);
         sol._redcost = _rationalLP->computeDualActivity(sol._dual) + _rationalLP->maxObj();
         sol._redcost *= -1;

         assert(sol._hasPrimal);
         assert(sol._hasDual);

         if( numAdjustedBounds > 0 )
         {
            MSG_INFO2( spxout << "Adjusted " << numAdjustedBounds << " nonbasic variables to bounds." << std::endl );
         }
      }
      while( true );

      // reset tolerances in floating-point solver
      _solver.setFeastol(rationalParam(SoPlex2::FEASTOL));
      _solver.setOpttol(rationalParam(SoPlex2::OPTTOL));
   }



   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void SoPlex2::_performUnboundedIRStable(SolRational& sol, bool& hasUnboundedRay, bool& stopped, bool& error)
   {
      bool primalFeasible;
      bool dualFeasible;
      bool infeasible;
      bool unbounded;

      // move objective function to constraints and adjust sides and bounds
      _transformUnbounded();

      // invalidate solution
      sol._invalidate();

      // perform iterative refinement
      _performOptIRStable(sol, false, false, primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

      // the unbounded problem should always be solved to optimality
      if( error || unbounded || infeasible || !primalFeasible || !dualFeasible )
      {
         sol._invalidate();
         hasUnboundedRay = false;
         stopped = false;
         error = true;
      }
      // or stopped due to some limit
      else if( stopped )
      {
         sol._invalidate();
         hasUnboundedRay = false;
         error = false;
      }
      else
      {
         Rational tau = sol._primal[numColsRational() - 1];

         MSG_DEBUG( spxout << "tau = " << tau << " (roughly " << rationalToString(tau) << ")" << std::endl );

         assert(tau <= 1 + 2 * rationalParam(SoPlex2::FEASTOL));
         assert(tau >= -rationalParam(SoPlex2::FEASTOL));

         // because the right-hand side and all bounds (but tau's upper bound) are zero, tau should be approximately
         // zero if basic; otherwise 0 or 1
         error = !(tau >= 1 || tau < rationalParam(SoPlex2::FEASTOL));
         assert(!error);

         hasUnboundedRay = (tau >= 1);

         sol._hasDual = false;
         if( !hasUnboundedRay )
            sol._hasPrimal = false;
      }

      // restore problem
      _untransformUnbounded(sol, hasUnboundedRay);
   }



   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void SoPlex2::_performFeasIRStable(SolRational& sol, bool& hasDualfarkas, bool& stopped, bool& error)
   {
      bool primalFeasible;
      bool dualFeasible;
      bool infeasible;
      bool unbounded;

      // remove objective function, shift, homogenize
      _transformFeasibility();

      // invalidate solution
      sol._invalidate();

      // perform iterative refinement
      _performOptIRStable(sol, false, false, primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

      // the feasibility problem should always be solved to optimality
      if( error || unbounded || infeasible || !primalFeasible || !dualFeasible )
      {
         sol._invalidate();
         hasDualfarkas = false;
         stopped = false;
         error = true;
      }
      // or stopped due to some limit
      else if( stopped )
      {
         sol._invalidate();
         hasDualfarkas = false;
         error = false;
      }
      else
      {
         Rational tau = sol._primal[numColsRational() - 1];

         MSG_DEBUG( spxout << "tau = " << tau << " (roughly " << rationalToString(tau) << ")" << std::endl );

         assert(tau >= -rationalParam(SoPlex2::FEASTOL));
         assert(tau <= 1 + rationalParam(SoPlex2::FEASTOL));

         error = (tau < -rationalParam(SoPlex2::FEASTOL) || tau > 1 + rationalParam(SoPlex2::FEASTOL));
         hasDualfarkas = (tau < 1);

         if( hasDualfarkas )
            sol._hasPrimal = false;
         else
            sol._hasDual = false;
      }

      // restore problem
      _untransformFeasibility(sol, hasDualfarkas);
   }



   /// introduces slack variables to transform inequality constraints into equations for both rational and real LP,
   /// which should be in sync
   void SoPlex2::_transformEquality()
   {
      MSG_DEBUG( spxout << "Transforming rows to equation form." << std::endl );

      // start timing
      _statistics->transformTime.start();

      MSG_DEBUG( _realLP->writeFile("beforeTransEqu.lp", 0, 0, 0) );

      // transform LP to minimization problem
      if( intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MAXIMIZE )
      {
         assert(_rationalLP->spxSense() == SPxLPRational::MAXIMIZE);
         assert(_realLP->spxSense() == SPxLPReal::MAXIMIZE);

         _rationalLP->changeObj(-(_rationalLP->maxObj()));
         _rationalLP->changeSense(SPxLPRational::MINIMIZE);

         _realLP->changeObj(-(_realLP->maxObj()));
         _realLP->changeSense(SPxLPReal::MINIMIZE);
      }

      // clear array of slack columns
      _slackCols.clear();

      // add artificial slack variables to convert inequality to equality constraints
      for( int i = 0; i < numRowsRational(); i++ )
      {
         if( lhsRational(i) != rhsRational(i) )
         {
            _slackCols.add(0.0, -rhsRational(i), DSVectorRational(UnitVector(i)), -lhsRational(i));
            _rationalLP->changeRange(i, 0.0, 0.0);
            _realLP->changeRange(i, 0.0, 0.0);
         }
      }

      _rationalLP->addCols(_slackCols);
      _realLP->addCols(_slackCols);

      // adjust basis
      if( _hasBasisRational )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numColsRational() - _slackCols.num() + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());

            switch( _basisStatusRowsRational[row] )
            {
            case SPxSolver::ON_LOWER:
               _basisStatusColsRational.append(SPxSolver::ON_UPPER);
               break;
            case SPxSolver::ON_UPPER:
               _basisStatusColsRational.append(SPxSolver::ON_LOWER);
               break;
            case SPxSolver::BASIC:
            case SPxSolver::FIXED:
            default:
               _basisStatusColsRational.append(_basisStatusRowsRational[row]);
               break;
            }

            _basisStatusRowsRational[row] = SPxSolver::FIXED;
         }
      }

      MSG_DEBUG( _realLP->writeFile("afterTransEqu.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime.stop();

      if( _slackCols.num() > 0 )
      {
         MSG_INFO1( spxout << "Added " << _slackCols.num() << " slack columns to transform rows to equality form." << std::endl );
      }
   }



   /// restores original problem
   void SoPlex2::_untransformEquality(SolRational& sol)
   {
      // start timing
      _statistics->transformTime.start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeUntransEqu.lp", 0, 0, 0) );

      int numCols = numColsRational();
      int numOrigCols = numColsRational() - _slackCols.num();

      // adjust solution
      if( sol.hasPrimal() )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());

            sol._slacks[row] -= sol._primal[col];
         }

         sol._primal.reDim(numOrigCols);
      }

      if( sol.hasPrimalray() )
      {
         sol._primalray.reDim(numOrigCols);
      }

      if( sol.hasDual() )
      {
         sol._redcost.reDim(numOrigCols);
      }

      // adjust basis
      if( _hasBasisRational )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());
            assert(_basisStatusRowsRational[row] == SPxSolver::FIXED || _basisStatusRowsRational[row] == SPxSolver::BASIC);

            if( _basisStatusRowsRational[row] == SPxSolver::FIXED )
            {
               switch( _basisStatusColsRational[col] )
               {
               case SPxSolver::ON_LOWER:
                  _basisStatusRowsRational[row] = SPxSolver::ON_UPPER;
                  break;
               case SPxSolver::ON_UPPER:
                  _basisStatusRowsRational[row] = SPxSolver::ON_LOWER;
                  break;
               case SPxSolver::BASIC:
               case SPxSolver::FIXED:
               default:
                  _basisStatusRowsRational[row] = _basisStatusColsRational[col];
                  break;
               }
            }
         }

         _basisStatusColsRational.reSize(numOrigCols);
      }

      // restore sides and remove slack columns
      for( int i = 0; i < _slackCols.num(); i++ )
      {
         int col = numOrigCols + i;
         int row = _slackCols.colVector(i).index(0);

         _rationalLP->changeRange(row, -upperRational(col), -lowerRational(col));
      }

      _rationalLP->removeColRange(numOrigCols, numCols - 1);
      _realLP->removeColRange(numOrigCols, numCols - 1);

      // restore original objective sense
      if( intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MAXIMIZE )
      {
         assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
         assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

         _rationalLP->changeObj(_rationalLP->maxObj());
         _rationalLP->changeSense(SPxLPRational::MAXIMIZE);
         _realLP->changeSense(SPxLPReal::MAXIMIZE);
      }

      // restore bounds and objective coefficients in real LP
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         _realLP->changeBounds(c, (Real)lowerRational(c), (Real)upperRational(c));
         _realLP->changeObj(c, (Real)objRational(c));
      }

      // restore sides in real LP
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         _realLP->changeRange(r, (Real)lhsRational(r), (Real)rhsRational(r));
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterUntransEqu.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime.stop();
   }



   /// transforms LP to unboundedness problem by moving the objective function to the constraints, changing right-hand
   /// side and bounds to zero, and adding an auxiliary variable for the decrease in the objective function
   void SoPlex2::_transformUnbounded()
   {
      assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
      assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

      MSG_INFO1( spxout << "Setting up LP to compute primal unbounded ray." << std::endl );

      // start timing
      _statistics->transformTime.start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeTransUnbounded.lp", 0, 0, 0) );

      // store right-hand side and bounds
      _unboundedSide = _rationalLP->rhs();
      _unboundedLower = _rationalLP->lower();
      _unboundedUpper = _rationalLP->upper();

      // make right-hand side zero
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         _rationalLP->changeRange(r, 0, 0);
         _realLP->changeRange(r, 0.0, 0.0);
      }

      // transform objective function to constraint and add auxiliary variable
      int numOrigCols = numColsRational();
      DSVectorRational obj(numOrigCols + 1);
      obj = _rationalLP->maxObj();
      obj *= -1;
      obj.add(numOrigCols, 1);
      _rationalLP->addRow(LPRowRational(0, obj, 0));
      _realLP->addRow(LPRowReal(0, DSVectorReal(obj), 0));

      assert(numColsRational() == numOrigCols + 1);

      // set objective coefficient and bounds for auxiliary variable
      _rationalLP->changeObj(numOrigCols, -1);
      _realLP->changeObj(numOrigCols, -1.0);

      _rationalLP->changeBounds(numOrigCols, 0, 1);
      _realLP->changeBounds(numOrigCols, 0.0, 1.0);

      // set objective coefficients to zero and adjust bounds for problem variables
      for( int c = numColsRational() - 2; c >= 0; c-- )
      {
         _rationalLP->changeObj(c, 0);
         _realLP->changeObj(c, 0.0);

         if( lowerRational(c) > -realParam(SoPlex2::INFTY) )
         {
            _rationalLP->changeLower(c, 0);
            _realLP->changeLower(c, 0.0);
         }

         if( upperRational(c) < realParam(SoPlex2::INFTY) )
         {
            _rationalLP->changeUpper(c, 0);
            _realLP->changeUpper(c, 0.0);
         }
      }

      // adjust basis
      if( _hasBasisRational )
      {
         _basisStatusColsRational.append(SPxSolver::ON_UPPER);
         _basisStatusRowsRational.append(SPxSolver::BASIC);
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterTransUnbounded.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime.stop();
   }



   /// undoes transformation to unboundedness problem
   void SoPlex2::_untransformUnbounded(SolRational& sol, bool unbounded)
   {
      // start timing
      _statistics->transformTime.start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeUntransUnbounded.lp", 0, 0, 0) );

      int numOrigCols = numColsRational() - 1;
      int numOrigRows = numRowsRational() - 1;

      // adjust solution and basis
      if( unbounded )
      {
         assert(sol._primal[numOrigCols] >= 1);

         sol._hasPrimal = false;
         sol._hasPrimalray = true;
         sol._hasDual = false;
         sol._hasDualfarkas = false;

         if( sol._primal[numOrigCols] != 1 )
            sol._primal /= sol._primal[numOrigCols];

         sol._primalray = sol._primal;
         sol._primalray.reDim(numOrigCols);

         _hasBasisRational = (_basisStatusColsRational[numOrigCols] != SPxSolver::BASIC && _basisStatusRowsRational[numOrigRows] == SPxSolver::BASIC);
         _basisStatusColsRational.reSize(numOrigCols);
         _basisStatusColsRational.reSize(numOrigRows);
      }
      else
      {
         sol._invalidate();
         _hasBasisRational = false;
      }

      // restore objective function
      const SVectorRational& rowVector = _rationalLP->rowVector(numOrigRows);
      DVectorRational obj(numOrigCols + 1);

      obj.clear();
      for( int i = rowVector.size() - 1; i >= 0; i-- )
         obj[rowVector.index(i)] = rowVector.value(i);

      DVectorReal objReal(obj);
      _rationalLP->changeObj(obj);
      _realLP->changeObj(objReal);

      // remove objective function constraint and auxiliary variable
      _rationalLP->removeRow(numOrigRows);
      _realLP->removeRow(numOrigRows);

      _rationalLP->removeCol(numOrigCols);
      _realLP->removeCol(numOrigCols);

      // restore right-hand side and bounds
      DVectorReal vectorReal(_unboundedSide);

      vectorReal = _unboundedSide;
      _rationalLP->changeRange(_unboundedSide, _unboundedSide);
      _realLP->changeRange(vectorReal, vectorReal);

      vectorReal = _unboundedLower;
      _rationalLP->changeLower(_unboundedLower);
      _realLP->changeLower(vectorReal);

      vectorReal = _unboundedUpper;
      _rationalLP->changeUpper(_unboundedUpper);
      _realLP->changeUpper(vectorReal);

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterUntransUnbounded.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime.stop();
   }



   /// transforms LP to feasibility problem by removing the objective function, shifting variables, and homogenizing the
   /// right-hand side
   void SoPlex2::_transformFeasibility()
   {
      assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
      assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

      MSG_INFO1( spxout << "Setting up LP to test for feasibility." << std::endl );

      // start timing
      _statistics->transformTime.start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeTransFeas.lp", 0, 0, 0) );

      // store objective function
      _feasObj.reDim(numColsRational());
      _rationalLP->getObj(_feasObj);

      // set objective coefficients to zero; shift primal space such as to guarantee that the zero solution is within
      // the bounds
      DVectorRational shiftedSide(rhsRational());
      _feasShiftValues.reDim(numColsRational());

      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         _rationalLP->changeObj(c, 0);
         _realLP->changeObj(c, 0.0);

         if( lowerRational(c) > 0 )
         {
            shiftedSide -= (colVectorRational(c) * lowerRational(c));
            _feasShiftValues[c] = lowerRational(c);
            _rationalLP->changeBounds(c, 0, upperRational(c) < realParam(SoPlex2::INFTY) ? upperRational(c) - lowerRational(c) : upperRational(c));
            _realLP->changeBounds(c, 0.0, (Real)upperRational(c));
         }
         else if( upperRational(c) < 0 )
         {
            shiftedSide -= (colVectorRational(c) * upperRational(c));
            _feasShiftValues[c] = upperRational(c);
            _rationalLP->changeBounds(c, lowerRational(c) > -realParam(SoPlex2::INFTY) ? lowerRational(c) - upperRational(c) : lowerRational(c), 0);
            _realLP->changeBounds(c, (Real)lowerRational(c), 0.0);
         }
         else
         {
            _feasShiftValues[c] = 0;
            _realLP->changeBounds(c, (Real)lowerRational(c), (Real)upperRational(c));
         }

         assert(lowerReal(c) <= upperReal(c));
      }

      // homogenize right-hand side
      SPxColId id;
      shiftedSide *= -1;
      _rationalLP->addCol(id, LPColRational(-1, DSVectorRational(shiftedSide), 1, 0));
      _realLP->addCol(id, LPColReal(-1.0, DSVectorReal(shiftedSide), 1.0, 0.0));

      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         _rationalLP->changeRange(r, 0, 0);
         _realLP->changeRange(r, 0.0, 0.0);
      }

      // adjust basis
      if( _hasBasisRational )
      {
         _basisStatusColsRational.append(SPxSolver::ON_UPPER);
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterTransFeas.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime.stop();
   }



   /// undoes transformation to feasibility problem
   void SoPlex2::_untransformFeasibility(SolRational& sol, bool infeasible)
   {
      // start timing
      _statistics->transformTime.start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeUntransFeas.lp", 0, 0, 0) );

      int numOrigCols = numColsRational() - 1;

      // adjust solution and basis
      if( infeasible )
      {
         assert(sol._hasDual);
         assert(sol._primal[numOrigCols] < 1);

         sol._hasPrimal = false;
         sol._hasPrimalray = false;
         sol._hasDual = false;
         sol._hasDualfarkas = true;

         sol._dualfarkas = sol._dual;

         _hasBasisRational = false;
         _basisStatusColsRational.reSize(numOrigCols);
      }
      else if( sol._hasPrimal )
      {
         assert(sol._primal[numOrigCols] >= 1);

         sol._hasPrimalray = false;
         sol._hasDual = false;
         sol._hasDualfarkas = false;

         if( sol._primal[numOrigCols] != 1 )
            sol._primal /= sol._primal[numOrigCols];

         sol._primal.reDim(numOrigCols);
         sol._slacks -= _rationalLP->colVector(numOrigCols);

         _hasBasisRational = (_basisStatusColsRational[numOrigCols] != SPxSolver::BASIC);
         _basisStatusColsRational.reSize(numOrigCols);
      }

      // unshift primal space and restore objective coefficients
      const SVectorRational& colVector = _rationalLP->colVector(numOrigCols);
      DVectorRational shiftedSide(numRowsRational());

      shiftedSide.clear();
      for( int i = colVector.size() - 1; i >= 0; i-- )
         shiftedSide[colVector.index(i)] = -colVector.value(i);

      for( int c = numOrigCols - 1; c >= 0; c-- )
      {
         Rational value = _feasShiftValues[c];

         assert(value == 0 || lowerRational(c) == 0 || upperRational(c) == 0);

         shiftedSide += (colVectorRational(c) * value);
         if( value != 0 )
         {
            _rationalLP->changeBounds(c, lowerRational(c) > -realParam(SoPlex2::INFTY) ? lowerRational(c) + value : lowerRational(c),
               upperRational(c) < realParam(SoPlex2::INFTY) ? upperRational(c) + value : upperRational(c));
         }

         _realLP->changeBounds(c, (Real)lowerRational(c), (Real)upperRational(c));

         _rationalLP->changeObj(c, _feasObj[c]);
         _realLP->changeObj(c, _feasObj[c]);

         assert(lowerReal(c) <= upperReal(c));
      }

      // restore right-hand side
      DVectorReal shiftedSideReal(shiftedSide);
      _rationalLP->changeRange(shiftedSide, shiftedSide);
      _realLP->changeRange(shiftedSideReal, shiftedSideReal);

      // remove last column
      _rationalLP->removeCol(numOrigCols);
      _realLP->removeCol(numOrigCols);

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterUntransFeas.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime.stop();

      assert(!sol._hasPrimal || sol._slacks == _rationalLP->computePrimalActivity(sol._primal));
   }

   /// compute radius of infeasibility box implied by an approximate Farkas' proof
   void SoPlex2::_computeInfeasBox(SolRational& sol)
   {
      // Given constraints of the form  b_l <= Ax <= b_r a farkas proof y should
      // satisfy y^TA = 0 and (y(+)^T b_l - y(-)^T b_r) > 0 where (+),(-) denote
      // the positive and negative parts of y.  If y is approximate, it may not
      // satisfy y^TA = 0 exactly, but the proof is still valid as long as the
      // following holds for all potentially feasible x:
      //
      //        y^TAx < (y(+)^T b_l - y(-)^T b_r)              (*)
      //
      // we may therefore calculate y^TA and (y(+)^T b_l - y(-)^T b_r) exactly
      // and check if the upper and lower bounds on x imply that all feasible x
      // satisfy (*), and if not then compute bounds on x to guarantee (*).
      // The simplest way to do this is to compute
      //        B = (y(+)^T b_l - y(-)^T b_r) / Sum_i(|(y^TA)_i|)
      // noting that if every component of x has |x_i| < B, then (*) holds.
      //
      // Note: A bound tighter than B can be computed by using the available
      // variable bound information.
      // The speed of this method can be further improved by using interval
      // arithmetic for all computations.
      // For related information see Sec. 4 of Neumaier and Shcherbina (2004 MPA)

      assert(sol.hasDualfarkas());

      Rational ytransb; // stores (y(+)^T b_l - y(-)^T b_r)
      DVectorRational ytransA(numColsRational());
      DVectorRational y(numRowsRational());
      Rational B;
      Rational temp;

      y = sol._dualfarkas;

      for( int r = 0; r < numRowsRational(); r++ )
      {
         ytransA.multAdd(y[r], _rationalLP->rowVector(r));
         (y[r] > 0) ? ytransb += y[r] * lhsRational(r) : ytransb += y[r] * rhsRational(r);
      }
      spxout << "ytransb = " << rationalToString(ytransb) << std::endl;

      temp = 0;
      for( int c = 0; c < numColsRational(); c++ )
      {
         (ytransA[c] > 0) ? temp += ytransA[c] * upperRational(c): temp += ytransA[c] * lowerRational(c);
      }
      spxout << "max ytransA*[x_l,x_u] = " << rationalToString(temp) << std::endl;

      if( temp < ytransb )
      {
         spxout << "Primal infeasible.  Farkas' proof verified exactly." << std::endl;
         return;
      }

      if( ytransb <= 0 )
      {
         spxout << "No infeasibility box available, Farkas' proof is no good." << std::endl;
         return;
      }

      temp = 0;
      for( int c = 0; c < numColsRational(); c++ )
      {
         temp += (ytransA[c] > 0) ? ytransA[c]: -ytransA[c];
      }
      assert(temp !=0);

      B = ytransb/temp;
      spxout << "Infeasibility box available, no feas. solutions with any compoents less than "
             << rationalToString(B) << std::endl;
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
