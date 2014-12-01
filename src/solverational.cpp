/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef SOPLEX_LEGACY
#include <iostream>
#include <assert.h>

#include "soplex.h"
#include "statistics.h"

namespace soplex
{
   /// solves rational LP
   void SoPlex::_solveRational()
   {
      bool hasUnboundedRay = false;
      bool infeasibilityNotCertified = false;
      bool unboundednessNotCertified = false;

      // start timing
      _statistics->solvingTime->start();
      _statistics->preprocessingTime->start();

      // remember that last solve was rational
      _lastSolveMode = SOLVEMODE_RATIONAL;

      // ensure that the solver has the original problem
      if( !_isRealLPLoaded )
      {
         assert(_realLP != &_solver);

         _solver.loadLP(*_realLP);
         spx_free(_realLP);
         _realLP = &_solver;
         _isRealLPLoaded = true;
      }
      // during the rational solve, we always store basis information in the basis arrays
      else if( _hasBasis )
      {
         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());
         _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
      }

      // deactivate objective limit in floating-point solver
      if( realParam(SoPlex::OBJLIMIT_LOWER) > -realParam(SoPlex::INFTY) || realParam(SoPlex::OBJLIMIT_UPPER) < realParam(SoPlex::INFTY) )
      {
         MSG_INFO2( spxout << "Deactivating objective limit.\n" );
      }

      _solver.setTerminationValue(realParam(SoPlex::INFTY));

      // introduce slack variables to transform inequality constraints into equations
      _transformEquality();
      _statistics->preprocessingTime->stop();

      // apply lifting to reduce range of nonzero matrix coefficients
      if( boolParam(SoPlex::LIFTING) )
         _lift();

      _statistics->preprocessingTime->stop();

      do
      {
         bool primalFeasible = false;
         bool dualFeasible = false;
         bool infeasible = false;
         bool unbounded = false;
         bool stopped = false;
         bool error = false;

         // solve problem with iterative refinement and recovery mechanism
         _performOptIRStable(_solRational, !unboundednessNotCertified, !infeasibilityNotCertified, 0,
            primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

         // case: an unrecoverable error occured
         if( error )
         {
            _status = SPxSolver::ERROR;
            break;
         }
         // case: stopped due to some limit
         else if( stopped )
         {
            _status = SPxSolver::ABORT_TIME;
            break;
         }
         // case: unboundedness detected for the first time
         else if( unbounded && !unboundednessNotCertified )
         {
            SolRational solUnbounded;

            _performUnboundedIRStable(solUnbounded, hasUnboundedRay, stopped, error);

            assert(!hasUnboundedRay || solUnbounded.hasPrimalRay());
            assert(!solUnbounded.hasPrimalRay() || hasUnboundedRay);

            if( error )
            {
               MSG_INFO1( spxout << "Error while testing for unboundedness.\n" );
               _status = SPxSolver::ERROR;
               break;
            }

            if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "Dual infeasible.  Primal unbounded ray available.\n" );
            }
            else
            {
               MSG_INFO1( spxout << "Dual feasible.  Rejecting primal unboundedness.\n" );
            }

            unboundednessNotCertified = !hasUnboundedRay;

            if( stopped )
            {
               _status = SPxSolver::ABORT_TIME;
               break;
            }

            _performFeasIRStable(_solRational, infeasible, stopped, error);

            if( hasUnboundedRay )
            {
               _solRational._primalRay = solUnbounded._primalRay;
               _solRational._hasPrimalRay = true;
            }

            if( error )
            {
               MSG_INFO1( spxout << "Error while testing for feasibility.\n" );
               _status = SPxSolver::ERROR;
               break;
            }
            else if( stopped )
            {
               _status = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "Primal infeasible.  Dual Farkas ray available.\n" );
               _status = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "Primal feasible and unbounded.\n" );
               _status = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "Primal feasible and bounded.\n" );
               continue;
            }
         }
         // case: infeasibility detected
         else if( infeasible && !infeasibilityNotCertified )
         {
            _performFeasIRStable(_solRational, infeasible, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "Error while testing for infeasibility.\n" );
               _status = SPxSolver::ERROR;
               break;
            }

            infeasibilityNotCertified = !infeasible;

            if( stopped )
            {
               _status = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "Primal infeasible.  Dual Farkas ray available.\n" );
               _status = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "Primal feasible and unbounded.\n" );
               _status = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "Primal feasible.  Optimizing again.\n" );
               continue;
            }
         }
         else if( primalFeasible && dualFeasible )
         {
            MSG_INFO1( spxout << "Solved to optimality.\n" );
            _status = SPxSolver::OPTIMAL;
            break;
         }
         else
         {
            MSG_INFO1( spxout << "Terminating without success.\n" );
            break;
         }
      }
      while( !_isSolveStopped() );

      ///@todo set status to ABORT_VALUE if optimal solution exceeds objective limit
      if( _isSolveStopped() )
         _status = SPxSolver::ABORT_TIME;

      if( _status == SPxSolver::OPTIMAL || _status == SPxSolver::INFEASIBLE || _status == SPxSolver::UNBOUNDED )
         _hasSolRational = true;

      // undo lifting
      if( boolParam(SoPlex::LIFTING) )
         _project(_solRational);

      // if the problem has been found to be infeasible and an approximate Farkas proof is available, we compute a
      // scaled unit box around the origin that provably contains no feasible solution; this currently only works for
      // equality form
#if 0
      if( _solRational.hasDualFarkas() )
         _computeInfeasBox(_solRational, false);
#endif

      // restore original problem
      _untransformEquality(_solRational);

      // since the real LP is loaded in the solver, we need to also pass the basis information to the solver if
      // available
      if( _hasBasis )
      {
         assert(_isRealLPLoaded);
         _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }

      // stop timing
      _statistics->solvingTime->stop();
   }



   /// solves current problem with iterative refinement and recovery mechanism
   void SoPlex::_performOptIRStable(SolRational& sol, bool acceptUnbounded, bool acceptInfeasible, int minRounds, bool& primalFeasible, bool& dualFeasible, bool& infeasible, bool& unbounded, bool& stopped, bool& error)
   {
      _statistics->rationalTime->start();
      primalFeasible = false;
      dualFeasible = false;
      infeasible = false;
      unbounded = false;
      stopped = false;
      error = false;

      // set working tolerances in floating-point solver
      _solver.setFeastol(realParam(SoPlex::FPFEASTOL));
      _solver.setOpttol(realParam(SoPlex::FPOPTTOL));

      // declare vectors and variables
      SPxSolver::Status result = SPxSolver::UNKNOWN;

      DVectorRational modLower(numColsRational());
      DVectorRational modUpper(numColsRational());
      DVectorRational modSide(numRowsRational());
      DVectorRational modObj(numColsRational());
      DVectorReal primalReal(numColsRational());
      DVectorReal dualReal(numRowsRational());

      Rational boundsViolation;
      Rational sideViolation;
      Rational redCostViolation;
      Rational primalScale;
      Rational dualScale;
      Rational maxScale;

      // solve original LP
      MSG_INFO1( spxout << "Initial floating-point solve . . .\n" );

      if( _hasBasis )
      {
         assert(_basisStatusRows.size() == numRowsRational());
         assert(_basisStatusCols.size() == numColsRational());
         _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
      }

      result = _solveRealStable(acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols);
      _statistics->rationalTime->stop();

      // evaluate result
      switch( result )
      {
      case SPxSolver::OPTIMAL:
         MSG_INFO1( spxout << "Floating-point optimal.\n" );
         break;
      case SPxSolver::INFEASIBLE:
         MSG_INFO1( spxout << "Floating-point infeasible.\n" );
         // the floating-point solve returns a Farkas ray if and only if the simplifier was not used, which is exactly
         // the case when a basis could be returned
         if( _hasBasis )
         {
            sol._dualFarkas = dualReal;
            sol._hasDualFarkas = true;
         }
         else
            sol._hasDualFarkas = false;
         infeasible = true;
         return;
      case SPxSolver::UNBOUNDED:
         MSG_INFO1( spxout << "Floating-point unbounded.\n" );
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

      _statistics->rationalTime->start();
      // corresponds to a minimization problem
      sol._primal = primalReal;
      _rationalLP->computePrimalActivity(sol._primal, sol._slacks);
      sol._hasPrimal = true;
      sol._primalObjVal = (sol._primal * _rationalLP->maxObj()) * -1;

      sol._dual = dualReal;
      _rationalLP->computeDualActivity(sol._dual, sol._redCost);
      sol._redCost += _rationalLP->maxObj();
      sol._redCost *= -1;
      sol._hasDual = true;
      sol._dualObjVal = sol._primalObjVal;

      _hasBasis = true;

      // initial scaling factors are one
      primalScale = 1;
      dualScale = 1;

      // control progress
      Real bestViolation = realParam(SoPlex::INFTY);
      int numFailedRefinements = 0;

      // refinement loop
      do
      {
         // decrement minRounds counter
         minRounds--;

         MSG_DEBUG( spxout << "Computing violations.\n" );

         // compute violation of bounds
         boundsViolation = 0;
         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            // lower bound
            modLower[c] = lowerRational(c);

            if( modLower[c] > -realParam(SoPlex::INFTY) )
               modLower[c] -= sol._primal[c];

            if( modLower[c] > boundsViolation )
               boundsViolation = modLower[c];

            // upper bound
            modUpper[c] = upperRational(c);

            if( modUpper[c] < realParam(SoPlex::INFTY) )
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
         _rationalLP->computeDualActivity(sol._dual, modObj);
         redCostViolation = 0;

         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            modObj[c] *= -1;
            modObj[c] += objRational(c);

            SPxSolver::VarStatus basisStatusCol = _basisStatusCols[c];

            if( basisStatusCol != SPxSolver::ON_UPPER && basisStatusCol != SPxSolver::FIXED && modObj[c] < -redCostViolation )
               redCostViolation = -modObj[c];

            if( basisStatusCol != SPxSolver::ON_LOWER && basisStatusCol != SPxSolver::FIXED && modObj[c] > redCostViolation )
               redCostViolation = modObj[c];
         }

         // output violations; the reduced cost violations for artificially introduced slack columns are actually violations of the dual multipliers
         MSG_INFO1( spxout
            << "Max. bound violation = " << rationalToString(boundsViolation) << "\n"
            << "Max. row violation = " << rationalToString(sideViolation) << "\n"
            << "Max. dual violation = " << rationalToString(redCostViolation) << "\n" );

         // terminate if tolerances are satisfied
         primalFeasible = (boundsViolation <= realParam(SoPlex::FEASTOL) && sideViolation <= realParam(SoPlex::FEASTOL));
         dualFeasible = (redCostViolation <= realParam(SoPlex::OPTTOL));
         if( primalFeasible && dualFeasible )
         {
            if( minRounds < 0 )
            {
               MSG_INFO1( spxout << "Tolerances reached.\n" );
               return;
            }
            else
            {
               MSG_INFO1( spxout << "Tolerances reached but minRounds forcing additional refinement rounds.\n" );
            }
         }

         // terminate if some limit is reached
         if( _isSolveStopped() )
         {
            stopped = true;
            return;
         }

         // check progress
         Rational sumMaxViolation = boundsViolation + sideViolation + redCostViolation;
         if( double(sumMaxViolation) > double(0.9 * bestViolation) )
         {
            MSG_INFO2( spxout << "Refinement failed to reduce violation significantly.\n" );
            numFailedRefinements++;
         }

         if( double(sumMaxViolation) < double(bestViolation) )
            bestViolation = Real(sumMaxViolation);

         if( numFailedRefinements >= 15 )
         {
            MSG_INFO1( spxout << "Giving up after 15 refinements without significantly increased precision.\n" );
            error = true;
            return;
         }

         ///@todo try rational reconstruction at geometric frequency

         // otherwise start refinement
         _statistics->refinements++;

         // compute primal scaling factor; limit increase in scaling by tolerance used in floating point solve
         maxScale = primalScale * Rational(realParam(SoPlex::MAXSCALEINCR));

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

         MSG_INFO2( spxout << "Scaling primal by " << rationalToString(primalScale) << ".\n" );

         // compute dual scaling factor; limit increase in scaling by tolerance used in floating point solve
         maxScale = dualScale * Rational(realParam(SoPlex::MAXSCALEINCR));

         dualScale = redCostViolation;
         assert(dualScale >= Rational(0));

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

         MSG_INFO2( spxout << "Scaling dual by " << rationalToString(dualScale) << ".\n" );

         // perform primal and dual scaling
         modLower *= primalScale;
         modUpper *= primalScale;
         modSide *= primalScale;
         modObj *= dualScale;

         // apply scaled bounds, side, and objective function
         _solver.changeBounds(DVectorReal(modLower), DVectorReal(modUpper));
         _solver.changeRange(DVectorReal(modSide), DVectorReal(modSide));
         _solver.changeObj(DVectorReal(modObj));

         MSG_INFO1( spxout << "Refined floating-point solve . . .\n" );

         // load basis
         _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());

         // solve modified problem
         result = _solveRealStable(acceptUnbounded, acceptInfeasible, primalReal, dualReal, _basisStatusRows, _basisStatusCols);
         _statistics->rationalTime->stop();

         // remember whether we moved to a new basis
         if( _solver.iterations() == 0 )
            _statistics->stallRefinements++;
         else
            _statistics->pivotRefinements = _statistics->refinements;

         // evaluate result; if modified problem was not solved to optimality, stop refinement
         switch( result )
         {
         case SPxSolver::OPTIMAL:
            MSG_INFO1( spxout << "Floating-point optimal.\n" );
            break;
         case SPxSolver::INFEASIBLE:
            MSG_INFO1( spxout << "Floating-point infeasible.\n" );
            sol._dualFarkas = dualReal;
            sol._hasDualFarkas = true;
            infeasible = true;
            return;
         case SPxSolver::UNBOUNDED:
            MSG_INFO1( spxout << "Floating-point unbounded.\n" );
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

         _statistics->rationalTime->start();
         MSG_DEBUG( spxout << "Correcting primal solution.\n" );

         int numAdjustedBounds = 0;
         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            sol._primal[c] += Rational(primalReal[c]) / primalScale;

            // force values of nonbasic variables to bounds
            SPxSolver::VarStatus basisStatusCol = _basisStatusCols[c];

            if( basisStatusCol == SPxSolver::ON_LOWER && sol._primal[c] != lowerRational(c) )
            {
               sol._primal[c] = lowerRational(c);
               numAdjustedBounds++;
            }
            else if( basisStatusCol == SPxSolver::ON_UPPER && sol._primal[c] != upperRational(c) )
            {
               sol._primal[c] = upperRational(c);
               numAdjustedBounds++;
            }
            else if( basisStatusCol == SPxSolver::FIXED )
            {
               assert(lowerRational(c) == upperRational(c));

               if( sol._primal[c] != lowerRational(c) )
               {
                  sol._primal[c] = lowerRational(c);
                  numAdjustedBounds++;
               }
            }
            else if( basisStatusCol == SPxSolver::ZERO && sol._primal[c] != 0 )
            {
               sol._primal[c] = 0;
               numAdjustedBounds++;
            }
         }

         // correct dual solution
         MSG_DEBUG( spxout << "Correcting dual solution.\n" );

         for( int r = numRowsRational() - 1; r >= 0; r-- )
         {
            sol._dual[r] += Rational(dualReal[r]) / dualScale;
         }

         // recompute slack and reduced cost values
         _rationalLP->computePrimalActivity(sol._primal, sol._slacks);
         _rationalLP->computeDualActivity(sol._dual, sol._redCost);
         sol._redCost += _rationalLP->maxObj();
         sol._redCost *= -1;

         assert(sol._hasPrimal);
         assert(sol._hasDual);

         if( numAdjustedBounds > 0 )
         {
            MSG_INFO2( spxout << "Adjusted " << numAdjustedBounds << " nonbasic variables to bounds.\n" );
         }
      }
      while( true );

      // compute objective values
      assert(sol._hasPrimal == sol._hasDual);
      if( sol._hasPrimal )
      {
         sol._primalObjVal = (sol._primal * _rationalLP->maxObj()) * -1;
         sol._dualObjVal = sol._primalObjVal;
      }

      // reset tolerances in floating-point solver
      _solver.setFeastol(realParam(SoPlex::FEASTOL));
      _solver.setOpttol(realParam(SoPlex::OPTTOL));
      _statistics->rationalTime->stop();
   }



   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void SoPlex::_performUnboundedIRStable(SolRational& sol, bool& hasUnboundedRay, bool& stopped, bool& error)
   {
      bool primalFeasible;
      bool dualFeasible;
      bool infeasible;
      bool unbounded;

      // move objective function to constraints and adjust sides and bounds
      _transformUnbounded();

      // invalidate solution
      sol.invalidate();

      // remember current number of refinements
      int oldRefinements = _statistics->refinements;

      // perform iterative refinement
      _performOptIRStable(sol, false, false, 0, primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

      // update unbounded refinement counter
      _statistics->unbdRefinements += _statistics->refinements - oldRefinements;

      // the unbounded problem should always be solved to optimality
      if( error || unbounded || infeasible || !primalFeasible || !dualFeasible )
      {
         sol.invalidate();
         hasUnboundedRay = false;
         stopped = false;
         error = true;
      }
      // or stopped due to some limit
      else if( stopped )
      {
         sol.invalidate();
         hasUnboundedRay = false;
         error = false;
      }
      else
      {
         Rational tau = sol._primal[numColsRational() - 1];

         MSG_DEBUG( spxout << "tau = " << tau << " (roughly " << rationalToString(tau) << ")\n" );

         assert(tau <= 1.0 + 2.0 * realParam(SoPlex::FEASTOL));
         assert(tau >= -realParam(SoPlex::FEASTOL));

         // because the right-hand side and all bounds (but tau's upper bound) are zero, tau should be approximately
         // zero if basic; otherwise 0 or 1
         error = !(tau >= Rational(1) || tau < realParam(SoPlex::FEASTOL));
         assert(!error);

         hasUnboundedRay = (tau >= Rational(1));

         sol._hasDual = false;
         if( !hasUnboundedRay )
            sol._hasPrimal = false;
      }

      // restore problem
      _untransformUnbounded(sol, hasUnboundedRay);
   }



   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void SoPlex::_performFeasIRStable(SolRational& sol, bool& withDualFarkas, bool& stopped, bool& error)
   {
      bool primalFeasible;
      bool dualFeasible;
      bool infeasible;
      bool unbounded;
      bool success = false;
      error = false;

#if 0
      ///@todo check whether approximate Farkas proof can be used
      _computeInfeasBox(_solRational, false);
      ///@todo if approx Farkas proof is good enough then exit without doing any transformation
#endif

      // remove objective function, shift, homogenize
      _transformFeasibility();

      // invalidate solution
      sol.invalidate();

      do
      {
         // remember current number of refinements
         int oldRefinements = _statistics->refinements;

         // perform iterative refinement
         _performOptIRStable(sol, false, false, 0, primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

         // update feasible refinement counter
         _statistics->feasRefinements += _statistics->refinements - oldRefinements;

         // the feasibility problem should always be solved to optimality
         if( error || unbounded || infeasible || !primalFeasible || !dualFeasible )
         {
            sol.invalidate();
            withDualFarkas = false;
            stopped = false;
            error = true;
         }
         // or stopped due to some limit
         else if( stopped )
         {
            sol.invalidate();
            withDualFarkas = false;
            error = false;
         }
         // else we should have either a refined Farkas proof or an approximate feasible solution to the original
         else
         {
            Rational tau = sol._primal[numColsRational() - 1];

            MSG_DEBUG( spxout << "tau = " << tau << " (roughly " << rationalToString(tau) << ")\n" );

            assert(tau >= -realParam(SoPlex::FEASTOL));
            assert(tau <= 1.0 + realParam(SoPlex::FEASTOL));

            error = (tau < -realParam(SoPlex::FEASTOL) || tau > Rational(1) + realParam(SoPlex::FEASTOL));
            withDualFarkas = (tau < Rational(1)); ///@todo shouldn't this use a tolerance? like (tau < 1-eps)? or even (tau <1/2)?

            if( withDualFarkas )
            {
               _solRational._hasDualFarkas = true;
               _solRational._dualFarkas = _solRational._dual;

#if 0
               // check if we can compute sufficiently large Farkas box
               _computeInfeasBox(_solRational, true);
#endif
               if( true )//@todo check if computeInfeasBox found a sufficient box
               {

                  success = true;
                  sol._hasPrimal = false;
               }
            }
            else
            {
               sol._hasDual = false;
               success = true; //successfully found approximate feasible solution
            }
         }
      }
      while(!error && !success);

      // restore problem
      _untransformFeasibility(sol, withDualFarkas);
   }



   /// reduces matrix coefficient in absolute value by the lifting procedure of Thiele et al. 2013
   void SoPlex::_lift()
   {
      MSG_DEBUG( spxout << "Reducing matrix coefficients by lifting.\n" );

      // start timing
      _statistics->transformTime->start();

      MSG_DEBUG( _realLP->writeFile("beforeLift.lp", 0, 0, 0) );

      // remember unlifted state
      _beforeLiftCols = numColsRational();
      _beforeLiftRows = numRowsRational();

      // allocate vector memory
      DSVectorRational colVector;
      SVectorRational::Element liftingRowMem[2];
      SVectorRational liftingRowVector(2, liftingRowMem);

      // search each column for large nonzeros entries
      const Rational maxValue = realParam(SoPlex::LIFTMAXVAL);

      for( int i = 0; i < numColsRational(); i++ )
      {
         MSG_DEBUG( spxout << "in lifting: examining column " << i << "\n" );

         // get column vector
         colVector = colVectorRational(i);

         bool addedLiftingRow = false;
         int liftingColumnIndex = -1;

         // go through nonzero entries of the column
         for( int k = colVector.size() - 1; k >= 0; k-- )
         {
            Rational value = colVector.value(k);

            if( abs(value) > maxValue )
            {
               MSG_DEBUG( spxout << "   --> nonzero " << k << " has value " << rationalToString(value) << " in row " << colVector.index(k) << "\n" );

               // add new column equal to maxValue times original column
               if( !addedLiftingRow )
               {
                  MSG_DEBUG( spxout << "            --> adding lifting row\n" );

                  assert(liftingRowVector.size() == 0);

                  liftingColumnIndex = numColsRational();
                  liftingRowVector.add(i, maxValue);
                  liftingRowVector.add(liftingColumnIndex, -1);

                  _rationalLP->addRow(LPRowRational(0, liftingRowVector, 0));
                  _realLP->addRow(LPRowReal(0.0, DSVectorReal(liftingRowVector), 0.0));

                  assert(liftingColumnIndex == numColsRational() - 1);
                  assert(liftingColumnIndex == numColsReal() - 1);

                  _rationalLP->changeBounds(liftingColumnIndex, -realParam(SoPlex::INFTY), realParam(SoPlex::INFTY));
                  _realLP->changeBounds(liftingColumnIndex, -realParam(SoPlex::INFTY), realParam(SoPlex::INFTY));

                  liftingRowVector.clear();
                  addedLiftingRow = true;
               }

               // get row index
               int rowIndex = colVector.index(k);
               assert(rowIndex >= 0);
               assert(rowIndex < _beforeLiftRows);
               assert(liftingColumnIndex == numColsRational() - 1);

               MSG_DEBUG( spxout << "            --> changing matrix\n" );

               // remove nonzero from original column
               _rationalLP->changeElement(rowIndex, i, 0);
               _realLP->changeElement(rowIndex, i, 0.0);

               // add nonzero divided by maxValue to new column
               Rational newValue = value / maxValue;
               _rationalLP->changeElement(rowIndex, liftingColumnIndex, newValue);
               _realLP->changeElement(rowIndex, liftingColumnIndex, Real(newValue));
            }
         }
      }

      // search each column for small nonzeros entries
      const Rational minValue = realParam(SoPlex::LIFTMINVAL);

      for( int i = 0; i < numColsRational(); i++ )
      {
         MSG_DEBUG( spxout << "in lifting: examining column " << i << "\n" );

         // get column vector
         colVector = colVectorRational(i);

         bool addedLiftingRow = false;
         int liftingColumnIndex = -1;

         // go through nonzero entries of the column
         for( int k = colVector.size() - 1; k >= 0; k-- )
         {
            Rational value = colVector.value(k);

            if( abs(value) < minValue )
            {
               MSG_DEBUG( spxout << "   --> nonzero " << k << " has value " << rationalToString(value) << " in row " << colVector.index(k) << "\n" );

               // add new column equal to maxValue times original column
               if( !addedLiftingRow )
               {
                  MSG_DEBUG( spxout << "            --> adding lifting row\n" );

                  assert(liftingRowVector.size() == 0);

                  liftingColumnIndex = numColsRational();
                  liftingRowVector.add(i, minValue);
                  liftingRowVector.add(liftingColumnIndex, -1);

                  _rationalLP->addRow(LPRowRational(0, liftingRowVector, 0));
                  _realLP->addRow(LPRowReal(0.0, DSVectorReal(liftingRowVector), 0.0));

                  assert(liftingColumnIndex == numColsRational() - 1);
                  assert(liftingColumnIndex == numColsReal() - 1);

                  _rationalLP->changeBounds(liftingColumnIndex, -realParam(SoPlex::INFTY), realParam(SoPlex::INFTY));
                  _realLP->changeBounds(liftingColumnIndex, -realParam(SoPlex::INFTY), realParam(SoPlex::INFTY));

                  liftingRowVector.clear();
                  addedLiftingRow = true;
               }

               // get row index
               int rowIndex = colVector.index(k);
               assert(rowIndex >= 0);
               assert(rowIndex < _beforeLiftRows);
               assert(liftingColumnIndex == numColsRational() - 1);

               MSG_DEBUG( spxout << "            --> changing matrix\n" );

               // remove nonzero from original column
               _rationalLP->changeElement(rowIndex, i, 0);
               _realLP->changeElement(rowIndex, i, 0.0);

               // add nonzero divided by maxValue to new column
               Rational newValue = value / minValue;
               _rationalLP->changeElement(rowIndex, liftingColumnIndex, newValue);
               _realLP->changeElement(rowIndex, liftingColumnIndex, Real(newValue));
            }
         }
      }

      // adjust basis
      if( _hasBasis )
      {
         assert(numColsRational() >= _beforeLiftCols);
         assert(numRowsRational() >= _beforeLiftRows);

         _basisStatusCols.append(numColsRational() - _beforeLiftCols, SPxSolver::BASIC);
         _basisStatusRows.append(numRowsRational() - _beforeLiftRows, SPxSolver::FIXED);
      }

      MSG_DEBUG( _realLP->writeFile("afterLift.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();

      if( numColsRational() > _beforeLiftCols || numRowsRational() > _beforeLiftRows )
      {
         MSG_INFO1( spxout << "Added " << numColsRational() - _beforeLiftCols << " columns and "
            << numRowsRational() - _beforeLiftRows << " rows to reduce large matrix coefficients\n." );
      }
   }



   /// undoes lifting
   void SoPlex::_project(SolRational& sol)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeProject.lp", 0, 0, 0) );

      assert(numColsRational() >= _beforeLiftCols);
      assert(numRowsRational() >= _beforeLiftRows);

      // shrink rational LP to original size
      _rationalLP->removeColRange(_beforeLiftCols, numColsRational() - 1);
      _rationalLP->removeRowRange(_beforeLiftRows, numRowsRational() - 1);

      // shrink real LP to original size
      _realLP->removeColRange(_beforeLiftCols, numColsReal() - 1);
      _realLP->removeRowRange(_beforeLiftRows, numRowsReal() - 1);

      // adjust solution
      if( sol.hasPrimal() )
      {
         sol._primal.reDim(_beforeLiftCols);
         sol._slacks.reDim(_beforeLiftRows);
      }

      if( sol.hasPrimalRay() )
      {
         sol._primalRay.reDim(_beforeLiftCols);
      }

      ///@todo if we know the mapping between original and lifting columns, we simply need to add the reduced cost of
      ///      the lifting column to the reduced cost of the original column; this is not implemented now, because for
      ///      optimal solutions the reduced costs of the lifting columns are zero
      const Rational maxValue = realParam(SoPlex::LIFTMAXVAL);

      for( int i = _beforeLiftCols; i < numColsRational() && sol._hasDual; i++ )
      {
         if( abs(maxValue * sol._redCost[i]) > realParam(SoPlex::OPTTOL) )
         {
            MSG_INFO1( spxout << "Warning: lost dual solution during project phase.\n" );
            sol._hasDual = false;
         }
      }

      if( sol.hasDual() )
      {
         sol._redCost.reDim(_beforeLiftCols);
         sol._dual.reDim(_beforeLiftRows);
      }

      if( sol.hasDualFarkas() )
      {
         sol._dualFarkas.reDim(_beforeLiftRows);
      }

      // adjust basis
      for( int i = _beforeLiftCols; i < numColsRational() && _hasBasis; i++ )
      {
         if( _basisStatusCols[i] != SPxSolver::BASIC )
         {
            MSG_INFO1( spxout << "Warning: lost basis during project phase because of nonbasic lifting column.\n" );
            _hasBasis = false;
         }
      }

      for( int i = _beforeLiftRows; i < numRowsRational() && _hasBasis; i++ )
      {
         if( _basisStatusRows[i] == SPxSolver::BASIC )
         {
            MSG_INFO1( spxout << "Warning: lost basis during project phase because of basic lifting row.\n" );
            _hasBasis = false;
         }
      }

      if( _hasBasis )
      {
         _basisStatusCols.reSize(_beforeLiftCols);
         _basisStatusRows.reSize(_beforeLiftRows);
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterProject.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// introduces slack variables to transform inequality constraints into equations for both rational and real LP,
   /// which should be in sync
   void SoPlex::_transformEquality()
   {
      MSG_DEBUG( spxout << "Transforming rows to equation form.\n" );

      // start timing
      _statistics->transformTime->start();

      MSG_DEBUG( _realLP->writeFile("beforeTransEqu.lp", 0, 0, 0) );

      // transform LP to minimization problem
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MAXIMIZE )
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
      if( _hasBasis )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());

            switch( _basisStatusRows[row] )
            {
            case SPxSolver::ON_LOWER:
               _basisStatusCols.append(SPxSolver::ON_UPPER);
               break;
            case SPxSolver::ON_UPPER:
               _basisStatusCols.append(SPxSolver::ON_LOWER);
               break;
            case SPxSolver::BASIC:
            case SPxSolver::FIXED:
            default:
               _basisStatusCols.append(_basisStatusRows[row]);
               break;
            }

            _basisStatusRows[row] = SPxSolver::FIXED;
         }
      }

      MSG_DEBUG( _realLP->writeFile("afterTransEqu.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();

      if( _slackCols.num() > 0 )
      {
         MSG_INFO1( spxout << "Added " << _slackCols.num() << " slack columns to transform rows to equality form.\n" );
      }
   }



   /// restores original problem
   void SoPlex::_untransformEquality(SolRational& sol)
   {
      // start timing
      _statistics->transformTime->start();

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

      if( sol.hasPrimalRay() )
      {
         sol._primalRay.reDim(numOrigCols);
      }

      if( sol.hasDual() )
      {
         sol._redCost.reDim(numOrigCols);
      }

      // adjust basis
      if( _hasBasis )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());
            assert(_basisStatusRows[row] == SPxSolver::FIXED || _basisStatusRows[row] == SPxSolver::BASIC);

            if( _basisStatusRows[row] == SPxSolver::FIXED )
            {
               switch( _basisStatusCols[col] )
               {
               case SPxSolver::ON_LOWER:
                  _basisStatusRows[row] = SPxSolver::ON_UPPER;
                  break;
               case SPxSolver::ON_UPPER:
                  _basisStatusRows[row] = SPxSolver::ON_LOWER;
                  break;
               case SPxSolver::BASIC:
               case SPxSolver::FIXED:
               default:
                  _basisStatusRows[row] = _basisStatusCols[col];
                  break;
               }
            }
         }

         _basisStatusCols.reSize(numOrigCols);
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
      if( intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MAXIMIZE )
      {
         assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
         assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

         _rationalLP->changeObj(_rationalLP->maxObj());
         _rationalLP->changeSense(SPxLPRational::MAXIMIZE);
         _realLP->changeSense(SPxLPReal::MAXIMIZE);

         sol._primalObjVal *= -1;
         sol._dualObjVal *= -1;
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
      _statistics->transformTime->stop();
   }



   /// transforms LP to unboundedness problem by moving the objective function to the constraints, changing right-hand
   /// side and bounds to zero, and adding an auxiliary variable for the decrease in the objective function
   void SoPlex::_transformUnbounded()
   {
      assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
      assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

      MSG_INFO1( spxout << "Setting up LP to compute primal unbounded ray.\n" );

      // start timing
      _statistics->transformTime->start();

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

         if( double(lowerRational(c)) > double(-realParam(SoPlex::INFTY)) )
         {
            _rationalLP->changeLower(c, 0);
            _realLP->changeLower(c, 0.0);
         }

         if( double(upperRational(c)) < double(realParam(SoPlex::INFTY)) )
         {
            _rationalLP->changeUpper(c, 0);
            _realLP->changeUpper(c, 0.0);
         }
      }

      // adjust basis
      if( _hasBasis )
      {
         _basisStatusCols.append(SPxSolver::ON_UPPER);
         _basisStatusRows.append(SPxSolver::BASIC);
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterTransUnbounded.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// undoes transformation to unboundedness problem
   void SoPlex::_untransformUnbounded(SolRational& sol, bool unbounded)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeUntransUnbounded.lp", 0, 0, 0) );

      int numOrigCols = numColsRational() - 1;
      int numOrigRows = numRowsRational() - 1;

      // adjust solution and basis
      if( unbounded )
      {
         assert(sol._primal[numOrigCols] >= Rational(1));

         sol._hasPrimal = false;
         sol._hasPrimalRay = true;
         sol._hasDual = false;
         sol._hasDualFarkas = false;

         if( sol._primal[numOrigCols] != Rational(1) )
            sol._primal /= sol._primal[numOrigCols];

         sol._primalRay = sol._primal;
         sol._primalRay.reDim(numOrigCols);

         _hasBasis = (_basisStatusCols[numOrigCols] != SPxSolver::BASIC && _basisStatusRows[numOrigRows] == SPxSolver::BASIC);
         _basisStatusCols.reSize(numOrigCols);
         _basisStatusCols.reSize(numOrigRows);
      }
      else
      {
         sol.invalidate();
         _hasBasis = false;
      }

      // restore objective function
      const SVectorRational& rowVector = _rationalLP->rowVector(numOrigRows);
      DVectorRational objCoefs(numOrigCols + 1);

      objCoefs.clear();
      for( int i = rowVector.size() - 1; i >= 0; i-- )
         objCoefs[rowVector.index(i)] = rowVector.value(i);

      DVectorReal objCoefsReal(objCoefs);
      _rationalLP->changeObj(objCoefs);
      _realLP->changeObj(objCoefsReal);

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
      _statistics->transformTime->stop();
   }



   /// transforms LP to feasibility problem by removing the objective function, shifting variables, and homogenizing the
   /// right-hand side
   void SoPlex::_transformFeasibility()
   {
      assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
      assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

      MSG_INFO1( spxout << "Setting up LP to test for feasibility.\n" );

      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeTransFeas.lp", 0, 0, 0) );

      // store objective function
      _feasObj.reDim(numColsRational());
      _rationalLP->getObj(_feasObj);

      // store right-hand side and bounds
      _feasSide = rhsRational();
      _feasLower = lowerRational();
      _feasUpper = upperRational();

      // set objective coefficients to zero; shift primal space such as to guarantee that the zero solution is within
      // the bounds
      DVectorRational shiftedSide(rhsRational());

      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         _rationalLP->changeObj(c, 0);
         _realLP->changeObj(c, 0.0);

         if( lowerRational(c) > Rational(0) )
         {
            shiftedSide -= (colVectorRational(c) * lowerRational(c));
            _rationalLP->changeBounds(c, 0, double(upperRational(c)) < double(realParam(SoPlex::INFTY)) ? upperRational(c) - lowerRational(c) : upperRational(c));
            _realLP->changeBounds(c, 0.0, (Real)upperRational(c));
         }
         else if( upperRational(c) < Rational(0) )
         {
            shiftedSide -= (colVectorRational(c) * upperRational(c));
            _rationalLP->changeBounds(c, double(lowerRational(c)) > double(-realParam(SoPlex::INFTY)) ? lowerRational(c) - upperRational(c) : lowerRational(c), 0);
            _realLP->changeBounds(c, (Real)lowerRational(c), 0.0);
         }
         else
         {
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
      if( _hasBasis )
      {
         _basisStatusCols.append(SPxSolver::ON_UPPER);
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterTransFeas.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// undoes transformation to feasibility problem
   void SoPlex::_untransformFeasibility(SolRational& sol, bool infeasible)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("beforeUntransFeas.lp", 0, 0, 0) );

      int numOrigCols = numColsRational() - 1;

      // adjust solution and basis
      if( infeasible )
      {
         assert(sol._hasDual);
         assert(sol._primal[numOrigCols] < Rational(1));

         sol._hasPrimal = false;
         sol._hasPrimalRay = false;
         sol._hasDual = false;
         sol._hasDualFarkas = true;

         sol._dualFarkas = sol._dual;

         _hasBasis = false;
         _basisStatusCols.reSize(numOrigCols);
      }
      else if( sol._hasPrimal )
      {
         assert(sol._primal[numOrigCols] >= Rational(1));

         sol._hasPrimalRay = false;
         sol._hasDual = false;
         sol._hasDualFarkas = false;

         if( sol._primal[numOrigCols] != Rational(1) )
            sol._primal /= sol._primal[numOrigCols];

         sol._primal.reDim(numOrigCols);
         sol._slacks -= _rationalLP->colVector(numOrigCols);

         _hasBasis = (_basisStatusCols[numOrigCols] != SPxSolver::BASIC);
         _basisStatusCols.reSize(numOrigCols);
      }
      else
         _hasBasis = false;

      // unshift primal space and restore objective coefficients
      for( int c = numOrigCols - 1; c >= 0; c-- )
      {
         assert(double(upperRational(c)) >= double(realParam(SoPlex::INFTY)) || double(lowerRational(c)) <= double(-realParam(SoPlex::INFTY))
            || _feasLower[c] - lowerRational(c) == _feasUpper[c] - upperRational(c));

         _rationalLP->changeBounds(c, _feasLower[c], _feasUpper[c]);
         _realLP->changeBounds(c, (Real)lowerRational(c), (Real)upperRational(c));

         _rationalLP->changeObj(c, _feasObj[c]);
         _realLP->changeObj(c, Real(_feasObj[c]));

         assert(lowerReal(c) <= upperReal(c));
      }

      // restore right-hand side
      DVectorReal feasSideReal(_feasSide);
      _rationalLP->changeRange(_feasSide, _feasSide);
      _realLP->changeRange(feasSideReal, feasSideReal);

      // remove last column
      _rationalLP->removeCol(numOrigCols);
      _realLP->removeCol(numOrigCols);

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFile("afterUntransFeas.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();

#ifndef NDEBUG
      if( sol._hasPrimal )
      {
         DVectorRational activity;
         _rationalLP->computePrimalActivity(sol._primal, activity);
         assert(sol._slacks == activity);
      }
#endif
   }

   /** computes radius of infeasibility box implied by an approximate Farkas' proof

    Given constraints of the form \f$ lhs <= Ax <= rhs \f$, a farkas proof y should satisfy \f$ y^T A = 0 \f$ and
    \f$ y_+^T lhs - y_-^T rhs > 0 \f$, where \f$ y_+, y_- \f$ denote the positive and negative parts of \f$ y \f$.
    If \f$ y \f$ is approximate, it may not satisfy \f$ y^T A = 0 \f$ exactly, but the proof is still valid as long
    as the following holds for all potentially feasible \f$ x \f$:

    \f[
       y^T Ax < (y_+^T lhs - y_-^T rhs)              (*)
    \f]

    we may therefore calculate \f$ y^T A \f$ and \f$ y_+^T lhs - y_-^T rhs \f$ exactly and check if the upper and lower
    bounds on \f$ x \f$ imply that all feasible \f$ x \f$ satisfy (*), and if not then compute bounds on \f$ x \f$ to
    guarantee (*).  The simplest way to do this is to compute

    \f[
       B = (y_+^T lhs - y_-^T rhs) / \sum_i(|(y^T A)_i|)
    \f]

    noting that if every component of \f$ x \f$ has \f$ |x_i| < B \f$, then (*) holds.

    \f$ B \f$ can be increased by iteratively including variable bounds smaller than \f$ B \f$.  The speed of this
    method can be further improved by using interval arithmetic for all computations.  For related information see
    Sec. 4 of Neumaier and Shcherbina, Mathematical Programming A, 2004.

    Set transformed to true if this method is called after _transformFeasibility().
   */
   void SoPlex::_computeInfeasBox(SolRational& sol, bool transformed)
   {
      assert(sol.hasDualFarkas());

      const VectorRational& lower = transformed ? _feasLower : lowerRational();
      const VectorRational& upper = transformed ? _feasUpper : upperRational();
      const VectorRational& lhs = transformed ? _feasSide : lhsRational();
      const VectorRational& rhs = transformed ? _feasSide : rhsRational();
      const VectorRational& y = sol._dualFarkas;

      const int numRows = numRowsRational();
      const int numCols = transformed ? numColsRational() - 1 : numColsRational();

      SSVectorRational ytransA(numColsRational());
      Rational ytransb;
      Rational temp;

      // prepare ytransA and ytransb; since we want exact arithmetic, we set the zero threshold of the semi-sparse
      // vector to zero
      ytransA.setEpsilon(0);
      ytransA.clear();
      ytransb = 0;

      ///@todo this currently works only if all constraints are equations
      // aggregate rows and sides using the multipliers of the Farkas ray
      for( int r = 0; r < numRows; r++ )
      {
         ytransA += y[r] * _rationalLP->rowVector(r);
         ytransb += y[r] * (y[r] > Rational(0) ? lhs[r] : rhs[r]);
      }

      // if we work on the feasibility problem, we ignore the last column
      if( transformed )
         ytransA.reDim(numCols);

      MSG_DEBUG( spxout << "ytransb = " << rationalToString(ytransb) << "\n" );

      // if we choose minus ytransb as vector of multipliers for the bound constraints on the variables, we obtain an
      // exactly feasible dual solution for the LP with zero objective function; we aggregate the bounds of the
      // variables accordingly and store its negation in temp
      temp = 0;
      bool isTempFinite = true;
      for( int c = 0; c < numCols && isTempFinite; c++ )
      {
         const Rational& minusRedCost = ytransA[c];

         if( minusRedCost > Rational(0) )
         {
            if( double(upper[c]) < double(realParam(SoPlex::INFTY)) )
               temp += minusRedCost * upper[c];
            else
               isTempFinite = false;
         }
         else if( minusRedCost < Rational(0) )
         {
            if( double(lower[c]) > double(-realParam(SoPlex::INFTY)) )
               temp += minusRedCost * lower[c];
            else
               isTempFinite = false;
         }
      }

      MSG_DEBUG( spxout << "max ytransA*[x_l,x_u] = " << (isTempFinite ? rationalToString(temp) : "infinite") << "\n" );

      // ytransb - temp is the increase in the dual objective along the Farkas ray; if this is positive, the dual is
      // unbounded and certifies primal infeasibility
      if( isTempFinite && temp < ytransb )
      {
         MSG_INFO1( spxout << "Farkas infeasibility proof verified exactly. (1)\n" );
         return;
      }

      // ensure that array of nonzero elements in ytransA is available
      assert(ytransA.isSetup());
      ytransA.setup();

      // if ytransb is negative, try to make it zero by including a positive lower bound or a negative upper bound
      if( ytransb < Rational(0) )
      {
         for( int c = 0; c < numCols; c++ )
         {
            if( lower[c] > Rational(0) )
            {
               ytransA.setValue(c, ytransA[c] - ytransb / lower[c]);
               ytransb = 0;
               break;
            }
            else if( upper[c] < Rational(0) )
            {
               ytransA.setValue(c, ytransA[c] - ytransb / upper[c]);
               ytransb = 0;
               break;
            }
         }
      }

      // if ytransb is still zero then the zero solution is inside the bounds and cannot be cut off by the Farkas
      // constraint; in this case, we cannot compute a Farkas box
      if( ytransb < Rational(0) )
      {
         MSG_INFO1( spxout << "Approximate Farkas proof to weak.  Could not compute Farkas box. (1)\n" );
         return;
      }

      // compute the one norm of ytransA
      temp = 0;
      const int size = ytransA.size();
      for( int n = 0; n < size; n++ )
         temp += abs(ytransA.value(n));

      // if the one norm is zero then ytransA is zero the Farkas proof should have been verified above
      assert(temp != Rational(0));

      // initialize variables in loop: size of Farkas box B, flag whether B has been increased, and number of current
      // nonzero in ytransA
      Rational B = ytransb / temp;
      bool success = false;
      int n = 0;

      // loop through nonzeros of ytransA
      MSG_DEBUG( spxout << "B = " << rationalToString(B) << "\n" );
      assert(ytransb >= Rational(0));

      while( true )
      {
         // if all nonzeros have been inspected once without increasing B, we abort; otherwise, we start another round
         if( n >= ytransA.size() )
         {
            if( !success )
               break;

            success = false;
            n = 0;
         }

         // get Farkas multiplier of the bound constraint as minus the value in ytransA
         const Rational& minusRedCost = ytransA.value(n);
         int colIdx = ytransA.index(n);

         // if the multiplier is positive we inspect the lower bound: if it is finite and within the Farkas box, we can
         // increase B by including it in the Farkas proof
         if( minusRedCost < Rational(0) && lower[colIdx] > -B && double(lower[colIdx]) > double(-realParam(SoPlex::INFTY)) )
         {
            ytransA.clearNum(n);
            ytransb -= minusRedCost * lower[colIdx];
            temp += minusRedCost;

            assert(ytransb >= Rational(0));
            assert(temp >= Rational(0));
            assert(temp == Rational(0) || ytransb / temp > B);

            // if ytransA and ytransb are zero, we have 0^T x >= 0 and cannot compute a Farkas box
            if( temp == Rational(0) && ytransb == Rational(0) )
            {
               MSG_INFO1( spxout << "Approximate Farkas proof to weak.  Could not compute Farkas box. (2)\n" );
               return;
            }
            // if ytransb is positive and ytransA is zero, we have 0^T x > 0, proving infeasibility
            else if( temp == Rational(0) )
            {
               assert(ytransb > Rational(0));
               MSG_INFO1( spxout << "Farkas infeasibility proof verified exactly. (2)\n" );
               return;
            }
            else
            {
               B = ytransb / temp;
               MSG_DEBUG( spxout << "B = " << rationalToString(B) << "\n" );
            }

            success = true;
         }
         // if the multiplier is negative we inspect the upper bound: if it is finite and within the Farkas box, we can
         // increase B by including it in the Farkas proof
         else if( minusRedCost > Rational(0) && upper[colIdx] < B && double(upper[colIdx]) < double(realParam(SoPlex::INFTY)) )
         {
            ytransA.clearNum(n);
            ytransb -= minusRedCost * upper[colIdx];
            temp -= minusRedCost;

            assert(ytransb >= Rational(0));
            assert(temp >= Rational(0));
            assert(temp == Rational(0) || ytransb / temp > B);

            // if ytransA and ytransb are zero, we have 0^T x >= 0 and cannot compute a Farkas box
            if( temp == Rational(0) && ytransb == Rational(0) )
            {
               MSG_INFO1( spxout << "Approximate Farkas proof to weak.  Could not compute Farkas box. (2)\n" );
               return;
            }
            // if ytransb is positive and ytransA is zero, we have 0^T x > 0, proving infeasibility
            else if( temp == Rational(0) )
            {
               assert(ytransb > Rational(0));
               MSG_INFO1( spxout << "Farkas infeasibility proof verified exactly. (2)\n" );
               return;
            }
            else
            {
               B = ytransb / temp;
               MSG_DEBUG( spxout << "B = " << rationalToString(B) << "\n" );
            }

            success = true;
         }
         // the multiplier is zero, we can ignore the bound constraints on this variable
         else if( minusRedCost == Rational(0) )
            ytransA.clearNum(n);
         // currently this bound cannot be used to increase B; we will check it again in the next round, because B might
         // have increased by then
         else
            n++;
      }

      if( B > Rational(0) )
      {
         MSG_INFO1( spxout << "Computed Farkas box: provably no feasible solutions with components less than "
            << rationalToString(B) << " in absolute value.\n" );
      }
   }



   /// solves real LP during iterative refinement
   SPxSolver::Status SoPlex::_solveRealForRational(bool fromscratch, VectorReal& primal, VectorReal& dual, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols)
   {
      assert(_isConsistent());

      assert(_solver.nRows() == numRowsRational());
      assert(_solver.nCols() == numColsRational());
      assert(primal.dim() == numColsRational());
      assert(dual.dim() == numRowsRational());

      SPxSolver::Status result = SPxSolver::UNKNOWN;

      if( fromscratch || !_hasBasis )
         _enableSimplifierAndScaler();
      else
         _disableSimplifierAndScaler();

      // start timing
      _statistics->syncTime->start();

      // if preprocessing is applied, we need to restore the original LP at the end
      SPxLPRational* rationalLP = 0;
      if( _simplifier != 0 || _scaler != 0 )
      {
         spx_alloc(rationalLP);
         rationalLP = new (rationalLP) SPxLPRational(_solver);
      }

      // stop timing
      _statistics->syncTime->stop();

      try
      {
         // apply problem simplification
         SPxSimplifier::Result simplificationStatus = SPxSimplifier::OKAY;
         if( _simplifier != 0 )
         {
            simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlex::EPSILON_ZERO), realParam(SoPlex::FPFEASTOL), realParam(SoPlex::FPOPTTOL));
         }

         // apply scaling after the simplification
         if( _scaler != 0 && simplificationStatus == SPxSimplifier::OKAY )
            _scaler->scale(_solver);

         // run the simplex method if problem has not been solved by the simplifier
         if( simplificationStatus == SPxSimplifier::OKAY )
         {
            MSG_INFO1( spxout << std::endl );

            _solveRealLPAndRecordStatistics();

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
                  DVectorReal tmpRedCost(vanished ? 0 : _solver.nCols());

                  if( !vanished )
                  {
                     assert(_solver.status() == SPxSolver::OPTIMAL);

                     _solver.getPrimal(tmpPrimal);
                     _solver.getSlacks(tmpSlacks);
                     _solver.getDual(tmpDual);
                     _solver.getRedCost(tmpRedCost);

                     // unscale vectors
                     if( _scaler != 0 )
                     {
                        _scaler->unscalePrimal(tmpPrimal);
                        _scaler->unscaleSlacks(tmpSlacks);
                        _scaler->unscaleDual(tmpDual);
                        _scaler->unscaleRedCost(tmpRedCost);
                     }

                     // get basis of transformed problem
                     _basisStatusRows.reSize(_solver.nRows());
                     _basisStatusCols.reSize(_solver.nCols());
                     _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
                  }

                  ///@todo catch exception
                  _simplifier->unsimplify(tmpPrimal, tmpDual, tmpSlacks, tmpRedCost, _basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

                  // store basis for original problem
                  basisStatusRows.reSize(numRowsReal());
                  basisStatusCols.reSize(numColsReal());
                  _simplifier->getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());

                  primal = _simplifier->unsimplifiedPrimal();
                  dual = _simplifier->unsimplifiedDual();
               }
               // if the original problem is not in the solver because of scaling, we also need to store the basis
               else
               {
                  _solver.getPrimal(primal);
                  _solver.getDual(dual);

                  // unscale vectors
                  if( _scaler != 0 )
                  {
                     _scaler->unscalePrimal(primal);
                     _scaler->unscaleDual(dual);
                  }

                  // get basis of transformed problem
                  basisStatusRows.reSize(_solver.nRows());
                  basisStatusCols.reSize(_solver.nCols());
                  _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());
               }
               break;

            case SPxSolver::ABORT_CYCLING:
            case SPxSolver::ABORT_TIME:
            case SPxSolver::ABORT_ITER:
            case SPxSolver::ABORT_VALUE:
            case SPxSolver::REGULAR:
            case SPxSolver::RUNNING:
            case SPxSolver::UNBOUNDED:
               break;
            case SPxSolver::INFEASIBLE:
               // if simplifier is active we cannot return a Farkas ray currently
               if( _simplifier != 0 )
                  break;

               // return Farkas ray as dual solution
               _solver.getDualfarkas(dual);

               // unscale vectors
               if( _scaler != 0 )
                  _scaler->unscaleDual(dual);

               // if the original problem is not in the solver because of scaling, we also need to store the basis
               basisStatusRows.reSize(_solver.nRows());
               basisStatusCols.reSize(_solver.nCols());
               _solver.getBasis(basisStatusRows.get_ptr(), basisStatusCols.get_ptr());

               break;

            case SPxSolver::INForUNBD:
            case SPxSolver::SINGULAR:
            default:
               _hasBasis = false;
               break;
            }
         }
      }
      catch( ... )
      {
         MSG_INFO1( spxout << "Exception thrown during floating-point solve.\n" );
         result = SPxSolver::ERROR;
      }

      // restore original LP if necessary
      if( _simplifier != 0 || _scaler != 0 )
      {
         assert(rationalLP != 0);
         _solver.loadLP((SPxLPReal)(*rationalLP));
         rationalLP->~SPxLPRational();
         spx_free(rationalLP);
      }

      return result;
   }



   /// solves real LP with recovery mechanism
   SPxSolver::Status SoPlex::_solveRealStable(bool acceptUnbounded, bool acceptInfeasible, VectorReal& primal, VectorReal& dual, DataArray< SPxSolver::VarStatus >& basisStatusRows, DataArray< SPxSolver::VarStatus >& basisStatusCols)
   {
      SPxSolver::Status result = SPxSolver::UNKNOWN;

      bool fromScratch = false;
      bool solved = false;
      bool solvedFromScratch = false;
      bool initialSolve = true;
      bool increasedMarkowitz = false;
      bool relaxedTolerances = false;
      bool tightenedTolerances = false;
      bool switchedScaler = false;
      bool switchedSimplifier = false;
      bool switchedRatiotester = false;
      bool switchedPricer = false;

      int ratiotester = intParam(SoPlex::RATIOTESTER);
      int pricer = intParam(SoPlex::PRICER);
      int simplifier = intParam(SoPlex::SIMPLIFIER);
      int scaler = intParam(SoPlex::SCALER);

      setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);

      while( !_isSolveStopped() )
      {
         assert(!increasedMarkowitz || GE(_slufactor.markowitz(), 0.9));

         result = _solveRealForRational(fromScratch, primal, dual, basisStatusRows, basisStatusCols);

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

         setIntParam(SoPlex::RATIOTESTER, ratiotester);
         setIntParam(SoPlex::PRICER, pricer);

         if( !switchedScaler )
         {
            MSG_INFO1( spxout << "Switching scaling." << std::endl );

            if( scaler == int(SoPlex::SCALER_OFF) )
               setIntParam(SoPlex::SCALER, SoPlex::SCALER_BIEQUI);
            else
               setIntParam(SoPlex::SCALER, SoPlex::SCALER_OFF);

            fromScratch = true;
            _solver.reLoad();

            solvedFromScratch = true;
            switchedScaler = true;
            continue;
         }

         if( !switchedSimplifier )
         {
            MSG_INFO1( spxout << "Switching simplification." << std::endl );

            if( simplifier == int(SoPlex::SIMPLIFIER_OFF) )
               setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_AUTO);
            else
               setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);

            fromScratch = true;
            _solver.reLoad();

            solvedFromScratch = true;
            switchedSimplifier = true;
            continue;
         }

         setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);

         if( !relaxedTolerances )
         {
            MSG_INFO1( spxout << "Relaxing tolerances." << std::endl );

            _solver.setType(_solver.rep() == SPxSolver::COLUMN ? SPxSolver::ENTER : SPxSolver::LEAVE);
            _solver.setDelta((_solver.feastol() * 1e3 > 1e-3) ? 1e-3 : (_solver.feastol() * 1e3));
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

      _solver.setFeastol(realParam(SoPlex::FPFEASTOL));
      _solver.setOpttol(realParam(SoPlex::FPOPTTOL));

      setIntParam(SoPlex::RATIOTESTER, ratiotester);
      setIntParam(SoPlex::PRICER, pricer);
      setIntParam(SoPlex::SIMPLIFIER, simplifier);
      setIntParam(SoPlex::SCALER, scaler);

      return result;
   }
} // namespace soplex
#endif
