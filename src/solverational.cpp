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
         _performOptIRStable(primalFeasible, dualFeasible, infeasible, unbounded, stopped, error);

         // case: an unrecoverable error occured
         if( error )
         {
            MSG_INFO1( spxout << "failed during feasibility problem\n" );
            _statusRational = SPxSolver::ERROR;
            break;
         }
         // case: stopped due to some limit
         else if( stopped )
         {
            MSG_INFO1( spxout << "stopped solving\n" );
            _statusRational = SPxSolver::ABORT_TIME;
            break;
         }
         // case: unboundedness detected for the first time
         else if( unbounded && !unboundednessNotCertified )
         {
            _performUnboundedIRStable(hasUnboundedRay, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "failed while trying to compute primal unbounded ray\n" );
               _statusRational = SPxSolver::ERROR;
               break;
            }

            assert(!hasUnboundedRay || _solRational.hasPrimalray());
            assert(!_solRational.hasPrimalray() || hasUnboundedRay);

            unboundednessNotCertified = !hasUnboundedRay;

            if( stopped )
            {
               MSG_INFO1( spxout << "stopped solving\n" );
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }

            _performFeasIRStable(infeasible, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "failed during feasibility problem\n" );
               _statusRational = SPxSolver::ERROR;
               break;
            }
            else if( stopped )
            {
               MSG_INFO1( spxout << "stopped solving\n" );
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "proved infeasiblity\n" );
               _statusRational = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "proved primal unboundedness\n" );
               _statusRational = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "unboundedness was not confirmed\n" );
               continue;
            }
         }
         // case: infeasibility detected
         else if( infeasible && !infeasibilityNotCertified )
         {
            _performFeasIRStable(infeasible, stopped, error);

            if( error )
            {
               MSG_INFO1( spxout << "failed during feasibility problem\n" );
               _statusRational = SPxSolver::ERROR;
               break;
            }

            infeasibilityNotCertified = !infeasible;

            if( stopped )
            {
               MSG_INFO1( spxout << "stopped solving\n" );
               _statusRational = SPxSolver::ABORT_TIME;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout << "proved infeasiblity\n" );
               _statusRational = SPxSolver::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout << "proved primal unboundedness\n" );
               _statusRational = SPxSolver::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout << "infeasibility was not confirmed\n" );
               continue;
            }
         }
         else if( primalFeasible && dualFeasible )
         {
            MSG_INFO1( spxout << "solved to optimality\n" );
            _statusRational = SPxSolver::OPTIMAL;
            break;
         }   
         else
         {
            MSG_INFO1( spxout << "an unknown error occured\n" );
            continue;
         }
      }
      while( !_isSolveStopped() );

      if( _isSolveStopped() )
         _statusRational = SPxSolver::ABORT_TIME;

      // restore original problem
      _untransformEquality();
   }



   /// introduces slack variables to transform inequality constraints into equations for both rational and real LP,
   /// which should be in sync
   void SoPlex2::_transformEquality()
   {
      // start timing
      _statistics->transformTime.start();

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

      MSG_DEBUG( spxout << "adding slack columns\n" );

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

      // stop timing
      _statistics->transformTime.stop();
   }



   /// restores original problem
   void SoPlex2::_untransformEquality()
   {
      // start timing
      _statistics->transformTime.start();

      int numCols = numColsRational();
      int numOrigCols = numColsRational() - _slackCols.num();

      // adjust solution
      if( _solRational.hasPrimal() )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());

            _solRational._slacks[row] -= _solRational._primal[col];
         }

         _solRational._primal.reDim(numOrigCols);
      }

      if( _solRational.hasPrimalray() )
      {
         _solRational._primalray.reDim(numOrigCols);
      }

      if( _solRational.hasDual() )
      {
         _solRational._redcost.reDim(numOrigCols);
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
         _realLP->changeRange(row, -upperRational(col), -lowerRational(col));
      }

      _rationalLP->removeColRange(numOrigCols, numCols - 1);
      _realLP->removeColRange(numOrigCols, numCols - 1);

      // transform LP to minimization problem
      if( intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MAXIMIZE )
      {
         assert(_rationalLP->spxSense() == SPxLPRational::MINIMIZE);
         assert(_realLP->spxSense() == SPxLPReal::MINIMIZE);

         _rationalLP->changeObj(_rationalLP->maxObj());
         _rationalLP->changeSense(SPxLPRational::MAXIMIZE);

         _realLP->changeObj(_realLP->maxObj());
         _realLP->changeSense(SPxLPReal::MAXIMIZE);
      }

      // stop timing
      _statistics->transformTime.stop();
   }



   /// solves current problem with iterative refinement and recovery mechanism
   void SoPlex2::_performOptIRStable(bool& primalFeasible, bool& dualFeasible, bool& infeasible, bool& unbounded, bool& stopped, bool& error)
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
      if( intParam(SoPlex2::ITERLIMIT) >= 0 )
         _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT) - _statistics->iterations);

      if( realParam(SoPlex2::TIMELIMIT) < realParam(SoPlex2::INFTY) )
         _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT) - _statistics->solvingTime.userTime());

      if( _hasBasisRational )
         _solver.setBasis(_basisStatusRowsRational.get_const_ptr(), _basisStatusColsRational.get_const_ptr());

      _solveRealStable();

      _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT));
      _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT));

      // evaluate result
      switch( _solver.status() )
      {
      case SPxSolver::OPTIMAL:
         break;
      case SPxSolver::INFEASIBLE:
         infeasible = true;
         return;
      case SPxSolver::UNBOUNDED:
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
      _solRational._primal = primalReal;
      _solRational._slacks = _rationalLP->computePrimalActivity(_solRational._primal);
      _solRational._hasPrimal = true;

      _solRational._dual = dualReal;
      _solRational._redcost = _rationalLP->computeDualActivity(_solRational._dual) + _rationalLP->maxObj();
      _solRational._redcost *= -1;
      _solRational._hasDual = true;

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

         MSG_DEBUG( spxout << "computing violations exactly . . .\n" );

         // compute violation of bounds
         boundsViolation = 0;

         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            // lower bound
            modLower[c] = lowerRational(c);

            if( modLower[c] > -realParam(SoPlex2::INFTY) )
               modLower[c] -= _solRational._primal[c];

            if( modLower[c] > boundsViolation )
               boundsViolation = modLower[c];

            // upper bound
            modUpper[c] = upperRational(c);

            if( modUpper[c] < realParam(SoPlex2::INFTY) )
               modUpper[c] -= _solRational._primal[c];

            if( modUpper[c] < -boundsViolation )
               boundsViolation = -modUpper[c];
         }

         // compute violation of sides
         modSide = _solRational._slacks;//.reDim(numRowsRational());
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
         modObj = _rationalLP->computeDualActivity(_solRational._dual);
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
         MSG_INFO1( spxout << std::endl );
         MSG_INFO1( spxout << "maximum violation: bounds=" << rationalToString(boundsViolation) << ", rows=" << rationalToString(sideViolation) << ", duals/redcosts=" << rationalToString(redcostViolation) << "\n" );

         ///@todo at this point the vectors primal_ex, modSide, dual_ex, modObj are computed and we may exit?

         // terminate if tolerances are satisfied
         primalFeasible = (boundsViolation <= rationalParam(SoPlex2::FEASTOL) && sideViolation <= rationalParam(SoPlex2::FEASTOL));
         dualFeasible = (redcostViolation <= rationalParam(SoPlex2::OPTTOL));
         if( primalFeasible && dualFeasible )
         {
            MSG_INFO1( spxout << "refinement finished: tolerances reached\n\n" );
            return;
         }

         // terminate if some limit is reached
         if( _isSolveStopped() )
         {
            MSG_INFO1( spxout << "refinement finished: limit reached\n\n" );
            stopped = true;
            return;
         }

         // otherwise start refinement
         _statistics->refinements++;

         MSG_INFO1( spxout << "starting refinement round " << _statistics->refinements + 1 << ": " );

         // compute primal scaling factor; limit increase in scaling by tolerance used in floating point solve
         maxScale = primalScale / Rational(realParam(SoPlex2::FPFEASTOL));

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

         MSG_INFO1( spxout << "scaling primal by " << rationalToString(primalScale) );

         // compute dual scaling factor; limit increase in scaling by tolerance used in floating point solve
         maxScale = dualScale / Rational(realParam(SoPlex2::FPOPTTOL));

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

         MSG_INFO1( spxout << ", scaling dual by " << rationalToString(dualScale) << " . . .\n\n" );

         // perform primal and dual scaling
         modLower *= primalScale;
         modUpper *= primalScale;
         modSide *= primalScale;
         modObj *= dualScale;

         // apply scaled bounds, side, and objective function
         _solver.changeBounds(DVectorReal(modLower), DVectorReal(modUpper));
         _solver.changeRange(DVectorReal(modSide), DVectorReal(modSide));
         _solver.changeObj(DVectorReal(modObj));

         MSG_DEBUG( spxout << "solving modified problem . . .\n" );

         // load basis
         _solver.setBasis(_basisStatusRowsRational.get_const_ptr(), _basisStatusColsRational.get_const_ptr());

         // solve modified problem
         if( intParam(SoPlex2::ITERLIMIT) >= 0 )
            _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT) - _statistics->iterations);

         if( realParam(SoPlex2::TIMELIMIT) < realParam(SoPlex2::INFTY) )
            _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT) - _statistics->solvingTime.userTime());

         _solveRealStable();

         _solver.setTerminationTime(realParam(SoPlex2::TIMELIMIT));
         _solver.setTerminationIter(intParam(SoPlex2::ITERLIMIT));

         // remember whether we moved to a new basis
         if( _solver.iterations() == 0 )
            _statistics->stallRefinements++;

         // evaluate result; if modified problem was not solved to optimality, stop refinement
         switch( _solver.status() )
         {
         case SPxSolver::OPTIMAL:
            break;
         case SPxSolver::INFEASIBLE:
            infeasible = true;
            return;
         case SPxSolver::UNBOUNDED:
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
         MSG_DEBUG( spxout << "correcting primal solution . . ." );

         int nadjusted = 0;
         for( int c = numColsRational() - 1; c >= 0; c-- )
         {
            _solRational._primal[c] += Rational(primalReal[c]) / primalScale;

            // force values of nonbasic variables to bounds
            SPxSolver::VarStatus basisStatus = _basisStatusColsRational[c];

            if( basisStatus == SPxSolver::ON_LOWER && _solRational._primal[c] != lowerRational(c) )
            {
               _solRational._primal[c] = lowerRational(c);
               nadjusted++;
            }
            else if( basisStatus == SPxSolver::ON_UPPER && _solRational._primal[c] != upperRational(c) )
            {
               _solRational._primal[c] = upperRational(c);
               nadjusted++;
            }
            else if( basisStatus == SPxSolver::FIXED )
            {
               assert(lowerRational(c) == upperRational(c));

               if( _solRational._primal[c] != lowerRational(c) )
               {
                  _solRational._primal[c] = lowerRational(c);
                  nadjusted++;
               }
            }
            else if( basisStatus == SPxSolver::ZERO && _solRational._primal[c] != 0 )
            {
               _solRational._primal[c] = 0;
               nadjusted++;
            }
         }

         // correct dual solution
         MSG_DEBUG( spxout << "correcting dual solution . . .\n" );

         for( int r = numRowsRational() - 1; r >= 0; r-- )
         {
            _solRational._dual[r] += Rational(dualReal[r]) / dualScale;
         }

         // recompute slack and reduced cost values
         _solRational._slacks = _rationalLP->computePrimalActivity(_solRational._primal);
         _solRational._redcost = _rationalLP->computeDualActivity(_solRational._dual) + _rationalLP->maxObj();
         _solRational._redcost *= -1;

         assert(_solRational._hasPrimal);
         assert(_solRational._hasDual);

         MSG_DEBUG( spxout << " adjusted " << nadjusted << " nonbasic variables to their bounds\n" );
      }
      while( true );

      // reset tolerances in floating-point solver
      _solver.setFeastol(rationalParam(SoPlex2::FEASTOL));
      _solver.setOpttol(rationalParam(SoPlex2::OPTTOL));
   }



   /// performs iterative refinement on the auxiliary problem for testing unboundedness
   void SoPlex2::_performUnboundedIRStable(bool& hasUnboundedRay, bool& stopped, bool& error)
   {
      hasUnboundedRay = false;
      stopped = false;
      error = false;
   }



   /// performs iterative refinement on the auxiliary problem for testing feasibility
   void SoPlex2::_performFeasIRStable(bool& infeasible, bool& stopped, bool& error)
   {
      infeasible = false;
      stopped = false;
      error = false;
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
