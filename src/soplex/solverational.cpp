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
#include <assert.h>

#include "soplex.h"
#include "soplex/statistics.h"
#include "soplex/slufactor_rational.h"
#include "soplex/ratrecon.h"

namespace soplex
{

   /// solves rational LP
  template <>
  void SoPlexBase<Real>::_optimizeRational()
   {
      bool hasUnboundedRay = false;
      bool infeasibilityNotCertified = false;
      bool unboundednessNotCertified = false;

      // start timing
      _statistics->solvingTime->start();
      _statistics->preprocessingTime->start();

      // remember that last solve was rational
      _lastSolveMode = SOLVEMODE_RATIONAL;

      // ensure that the solver has the original problemo
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
         _basisStatusRows.reSize(numRows());
         _basisStatusCols.reSize(numCols());
         _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr(), _basisStatusRows.size(), _basisStatusCols.size());
      }

      // store objective, bounds, and sides of real LP in case they will be modified during iterative refinement
      _storeLPReal();

      // deactivate objective limit in floating-point solver
      if( realParam(SoPlexBase<Real>::OBJLIMIT_LOWER) > -realParam(SoPlexBase<Real>::INFTY) || realParam(SoPlexBase<Real>::OBJLIMIT_UPPER) < realParam(SoPlexBase<Real>::INFTY) )
      {
         MSG_INFO2( spxout, spxout << "Deactivating objective limit.\n" );
      }

      _solver.setTerminationValue(realParam(SoPlexBase<Real>::INFTY));

      _statistics->preprocessingTime->stop();

      // apply lifting to reduce range of nonzero matrix coefficients
      if( boolParam(SoPlexBase<Real>::LIFTING) )
         _lift();

      // force column representation
      ///@todo implement row objectives with row representation
      int oldRepresentation = intParam(SoPlexBase<Real>::REPRESENTATION);
      setIntParam(SoPlexBase<Real>::REPRESENTATION, SoPlexBase<Real>::REPRESENTATION_COLUMN);

      // force ratio test (avoid bound flipping)
      int oldRatiotester = intParam(SoPlexBase<Real>::RATIOTESTER);
      setIntParam(SoPlexBase<Real>::RATIOTESTER, SoPlexBase<Real>::RATIOTESTER_FAST);

      ///@todo implement handling of row objectives in Cplex interface
#ifdef SOPLEX_WITH_CPX
      int oldEqtrans = boolParam(SoPlexBase<Real>::EQTRANS);
      setBoolParam(SoPlexBase<Real>::EQTRANS, true);
#endif

      // introduce slack variables to transform inequality constraints into equations
      if( boolParam(SoPlexBase<Real>::EQTRANS) )
         _transformEquality();

      _storedBasis = false;

      bool stoppedTime;
      bool stoppedIter;
      do
      {
         bool primalFeasible = false;
         bool dualFeasible = false;
         bool infeasible = false;
         bool unbounded = false;
         bool error = false;
         stoppedTime = false;
         stoppedIter = false;

         // solve problem with iterative refinement and recovery mechanism
         _performOptIRStable(_solRational, !unboundednessNotCertified, !infeasibilityNotCertified, 0,
            primalFeasible, dualFeasible, infeasible, unbounded, stoppedTime, stoppedIter, error);

         // case: an unrecoverable error occured
         if( error )
         {
            _status = SPxSolverBase<Real>::ERROR;
            break;
         }
         // case: stopped due to some limit
         else if( stoppedTime )
         {
            _status = SPxSolverBase<Real>::ABORT_TIME;
            break;
         }
         else if(  stoppedIter )
         {
            _status = SPxSolverBase<Real>::ABORT_ITER;
            break;
         }
         // case: unboundedness detected for the first time
         else if( unbounded && !unboundednessNotCertified )
         {
            SolRational solUnbounded;

            _performUnboundedIRStable(solUnbounded, hasUnboundedRay, stoppedTime, stoppedIter, error);

            assert(!hasUnboundedRay || solUnbounded.hasPrimalRay());
            assert(!solUnbounded.hasPrimalRay() || hasUnboundedRay);

            if( error )
            {
               MSG_INFO1( spxout, spxout << "Error while testing for unboundedness.\n" );
               _status = SPxSolverBase<Real>::ERROR;
               break;
            }

            if( hasUnboundedRay )
            {
               MSG_INFO1( spxout, spxout << "Dual infeasible.  Primal unbounded ray available.\n" );
            }
            else
            {
               MSG_INFO1( spxout, spxout << "Dual feasible.  Rejecting primal unboundedness.\n" );
            }

            unboundednessNotCertified = !hasUnboundedRay;

            if( stoppedTime )
            {
               _status = SPxSolverBase<Real>::ABORT_TIME;
               break;
            }
            else if( stoppedIter )
            {
               _status = SPxSolverBase<Real>::ABORT_ITER;
               break;
            }

            _performFeasIRStable(_solRational, infeasible, stoppedTime, stoppedIter, error);

            ///@todo this should be stored already earlier, possible switch use solRational above and solFeas here
            if( hasUnboundedRay )
            {
               _solRational._primalRay = solUnbounded._primalRay;
               _solRational._hasPrimalRay = true;
            }

            if( error )
            {
               MSG_INFO1( spxout, spxout << "Error while testing for feasibility.\n" );
               _status = SPxSolverBase<Real>::ERROR;
               break;
            }
            else if( stoppedTime )
            {
               _status = SPxSolverBase<Real>::ABORT_TIME;
               break;
            }
            else if( stoppedIter )
            {
               _status = SPxSolverBase<Real>::ABORT_ITER;
               break;
            }
            else if( infeasible )
            {
               MSG_INFO1( spxout, spxout << "Primal infeasible.  Dual Farkas ray available.\n" );
               _status = SPxSolverBase<Real>::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout, spxout << "Primal feasible and unbounded.\n" );
               _status = SPxSolverBase<Real>::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout, spxout << "Primal feasible and bounded.\n" );
               continue;
            }
         }
         // case: infeasibility detected
         else if( infeasible && !infeasibilityNotCertified )
         {
            _storeBasis();

            _performFeasIRStable(_solRational, infeasible, stoppedTime, stoppedIter, error);

            if( error )
            {
               MSG_INFO1( spxout, spxout << "Error while testing for infeasibility.\n" );
               _status = SPxSolverBase<Real>::ERROR;
               _restoreBasis();
               break;
            }

            infeasibilityNotCertified = !infeasible;

            if( stoppedTime )
            {
               _status = SPxSolverBase<Real>::ABORT_TIME;
               _restoreBasis();
               break;
            }
            else if( stoppedIter )
            {
               _status = SPxSolverBase<Real>::ABORT_ITER;
               _restoreBasis();
               break;
            }

            if( infeasible && boolParam(SoPlexBase<Real>::TESTDUALINF) )
            {
               SolRational solUnbounded;

               _performUnboundedIRStable(solUnbounded, hasUnboundedRay, stoppedTime, stoppedIter, error);

               assert(!hasUnboundedRay || solUnbounded.hasPrimalRay());
               assert(!solUnbounded.hasPrimalRay() || hasUnboundedRay);

               if( error )
               {
                  MSG_INFO1( spxout, spxout << "Error while testing for dual infeasibility.\n" );
                  _status = SPxSolverBase<Real>::ERROR;
                  _restoreBasis();
                  break;
               }

               if( hasUnboundedRay )
               {
                  MSG_INFO1( spxout, spxout << "Dual infeasible.  Primal unbounded ray available.\n" );
                  _solRational._primalRay = solUnbounded._primalRay;
                  _solRational._hasPrimalRay = true;
               }
               else if( solUnbounded._isDualFeasible )
               {
                  MSG_INFO1( spxout, spxout << "Dual feasible.  Storing dual multipliers.\n" );
                  _solRational._dual = solUnbounded._dual;
                  _solRational._redCost = solUnbounded._redCost;
                  _solRational._isDualFeasible = true;
               }
               else
               {
                  assert(false);
                  MSG_INFO1( spxout, spxout << "Not dual infeasible.\n" );
               }
            }

            _restoreBasis();

            if( infeasible )
            {
               MSG_INFO1( spxout, spxout << "Primal infeasible.  Dual Farkas ray available.\n" );
               _status = SPxSolverBase<Real>::INFEASIBLE;
               break;
            }
            else if( hasUnboundedRay )
            {
               MSG_INFO1( spxout, spxout << "Primal feasible and unbounded.\n" );
               _status = SPxSolverBase<Real>::UNBOUNDED;
               break;
            }
            else
            {
               MSG_INFO1( spxout, spxout << "Primal feasible.  Optimizing again.\n" );
               continue;
            }
         }
         else if( primalFeasible && dualFeasible )
         {
            MSG_INFO1( spxout, spxout << "Solved to optimality.\n" );
            _status = SPxSolverBase<Real>::OPTIMAL;
            break;
         }
         else
         {
            MSG_INFO1( spxout, spxout << "Terminating without success.\n" );
            break;
         }
      }
      while( !_isSolveStopped(stoppedTime, stoppedIter) );

      ///@todo set status to ABORT_VALUE if optimal solution exceeds objective limit

      if( _status == SPxSolverBase<Real>::OPTIMAL || _status == SPxSolverBase<Real>::INFEASIBLE || _status == SPxSolverBase<Real>::UNBOUNDED )
         _hasSolRational = true;

      // restore original problem
      if( boolParam(SoPlexBase<Real>::EQTRANS) )
         _untransformEquality(_solRational);

#ifdef SOPLEX_WITH_CPX
      setBoolParam(SoPlexBase<Real>::EQTRANS, oldEqtrans);
#endif

      // reset representation and ratio test
      setIntParam(SoPlexBase<Real>::REPRESENTATION, oldRepresentation);
      setIntParam(SoPlexBase<Real>::RATIOTESTER, oldRatiotester);

      // undo lifting
      if( boolParam(SoPlexBase<Real>::LIFTING) )
         _project(_solRational);

      // restore objective, bounds, and sides of real LP in case they have been modified during iterative refinement
      _restoreLPReal();

      // since the real LP is loaded in the solver, we need to also pass the basis information to the solver if
      // available
      if( _hasBasis )
      {
         assert(_isRealLPLoaded);
         _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         _hasBasis = (_solver.basis().status() > SPxBasisBase<Real>::NO_PROBLEM);
      }

      // stop timing
      _statistics->solvingTime->stop();
   }



  


   /// performs iterative refinement on the auxiliary problem for testing unboundedness
  template <class R>
  void SoPlexBase<R>::_performUnboundedIRStable(
      SolRational& sol,
      bool& hasUnboundedRay,
      bool& stopped,
      bool& stoppedIter,
      bool& error)
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
      _performOptIRStable(sol, false, false, 0, primalFeasible, dualFeasible, infeasible, unbounded, stopped, stoppedIter, error);

      // update unbounded refinement counter
      _statistics->unbdRefinements += _statistics->refinements - oldRefinements;

      // stopped due to some limit
      if( stopped )
      {
         sol.invalidate();
         hasUnboundedRay = false;
         error = false;
      }
      // the unbounded problem should always be solved to optimality
      else if( error || unbounded || infeasible || !primalFeasible || !dualFeasible )
      {
         sol.invalidate();
         hasUnboundedRay = false;
         stopped = false;
         error = true;
      }
      else
      {
         const Rational& tau = sol._primal[numColsRational() - 1];

         MSG_DEBUG( std::cout << "tau = " << tau << " (roughly " << rationalToString(tau) << ")\n" );

         assert(tau <= 1.0 + 2.0 * realParam(SoPlexBase<R>::FEASTOL));
         assert(tau >= -realParam(SoPlexBase<R>::FEASTOL));

         // because the right-hand side and all bounds (but tau's upper bound) are zero, tau should be approximately
         // zero if basic; otherwise at its upper bound 1
         error = !(tau >= _rationalPosone || tau <= _rationalFeastol);
         assert(!error);

         hasUnboundedRay = (tau >= 1);
      }

      // restore problem
      _untransformUnbounded(sol, hasUnboundedRay);
   }



   /// performs iterative refinement on the auxiliary problem for testing feasibility
  template <class R>
  void SoPlexBase<R>::_performFeasIRStable(
      SolRational& sol,
      bool& withDualFarkas,
      bool& stopped,
      bool& stoppedIter,
      bool& error)
   {
      bool primalFeasible;
      bool dualFeasible;
      bool infeasible;
      bool unbounded;
      bool success = false;
      error = false;

#if 0
      // if the problem has been found to be infeasible and an approximate Farkas proof is available, we compute a
      // scaled unit box around the origin that provably contains no feasible solution; this currently only works for
      // equality form
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
         _performOptIRStable(sol, false, false, 0, primalFeasible, dualFeasible, infeasible, unbounded, stopped, stoppedIter, error);

         // update feasible refinement counter
         _statistics->feasRefinements += _statistics->refinements - oldRefinements;

         // stopped due to some limit
         if( stopped )
         {
            sol.invalidate();
            withDualFarkas = false;
            error = false;
         }
         // the feasibility problem should always be solved to optimality
         else if( error || unbounded || infeasible || !primalFeasible || !dualFeasible )
         {
            sol.invalidate();
            withDualFarkas = false;
            stopped = false;
            error = true;
         }
         // else we should have either a refined Farkas proof or an approximate feasible solution to the original
         else
         {
            const Rational& tau = sol._primal[numColsRational() - 1];

            MSG_DEBUG( std::cout << "tau = " << tau << " (roughly " << rationalToString(tau) << ")\n" );

            assert(tau >= -realParam(SoPlexBase<R>::FEASTOL));
            assert(tau <= 1.0 + realParam(SoPlexBase<R>::FEASTOL));

            error = (tau < -_rationalFeastol || tau > _rationalPosone + _rationalFeastol);
            withDualFarkas = (tau < _rationalPosone);

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
                  sol._isPrimalFeasible= false;
               }
            }
            else
            {
               sol._isDualFeasible = false;
               success = true; //successfully found approximate feasible solution
            }
         }
      }
      while(!error && !success && !stopped);

      // restore problem
      _untransformFeasibility(sol, withDualFarkas);
   }



   /// reduces matrix coefficient in absolute value by the lifting procedure of Thiele et al. 2013
  template <class R>
   void SoPlexBase<R>::_lift()
   {
      MSG_DEBUG( std::cout << "Reducing matrix coefficients by lifting.\n" );

      // start timing
      _statistics->transformTime->start();

      MSG_DEBUG( _realLP->writeFileLPBase("beforeLift.lp", 0, 0, 0) );

      // remember unlifted state
      _beforeLiftCols = numColsRational();
      _beforeLiftRows = numRowsRational();

      // allocate vector memory
      DSVectorRational colVector;
      SVectorRational::Element liftingRowMem[2];
      SVectorRational liftingRowVector(2, liftingRowMem);

      // search each column for large nonzeros entries
      const Rational maxValue = realParam(SoPlexBase<R>::LIFTMAXVAL);

      for( int i = 0; i < numColsRational(); i++ )
      {
         MSG_DEBUG( std::cout << "in lifting: examining column " << i << "\n" );

         // get column vector
         colVector = colVectorRational(i);

         bool addedLiftingRow = false;
         int liftingColumnIndex = -1;

         // go through nonzero entries of the column
         for( int k = colVector.size() - 1; k >= 0; k-- )
         {
            const Rational& value = colVector.value(k);

            if( spxAbs(value) > maxValue )
            {
               MSG_DEBUG( std::cout << "   --> nonzero " << k << " has value " << rationalToString(value) << " in row " << colVector.index(k) << "\n" );

               // add new column equal to maxValue times original column
               if( !addedLiftingRow )
               {
                  MSG_DEBUG( std::cout << "            --> adding lifting row\n" );

                  assert(liftingRowVector.size() == 0);

                  liftingColumnIndex = numColsRational();
                  liftingRowVector.add(i, maxValue);
                  liftingRowVector.add(liftingColumnIndex, -1);

                  _rationalLP->addRow(LPRowRational(0, liftingRowVector, 0));
                  _realLP->addRow(LPRowReal(0.0, DSVectorReal(liftingRowVector), 0.0));

                  assert(liftingColumnIndex == numColsRational() - 1);
                  assert(liftingColumnIndex == numCols() - 1);

                  _rationalLP->changeBounds(liftingColumnIndex, _rationalNegInfty, _rationalPosInfty);
                  _realLP->changeBounds(liftingColumnIndex, -realParam(SoPlexBase<R>::INFTY), realParam(SoPlexBase<R>::INFTY));

                  liftingRowVector.clear();
                  addedLiftingRow = true;
               }

               // get row index
               int rowIndex = colVector.index(k);
               assert(rowIndex >= 0);
               assert(rowIndex < _beforeLiftRows);
               assert(liftingColumnIndex == numColsRational() - 1);

               MSG_DEBUG( std::cout << "            --> changing matrix\n" );

               // remove nonzero from original column
               _rationalLP->changeElement(rowIndex, i, 0);
               _realLP->changeElement(rowIndex, i, 0.0);

               // add nonzero divided by maxValue to new column
               Rational newValue(value);
               newValue /= maxValue;
               _rationalLP->changeElement(rowIndex, liftingColumnIndex, newValue);
               _realLP->changeElement(rowIndex, liftingColumnIndex, Real(newValue));
            }
         }
      }

      // search each column for small nonzeros entries
      const Rational minValue = realParam(SoPlexBase<R>::LIFTMINVAL);

      for( int i = 0; i < numColsRational(); i++ )
      {
         MSG_DEBUG( std::cout << "in lifting: examining column " << i << "\n" );

         // get column vector
         colVector = colVectorRational(i);

         bool addedLiftingRow = false;
         int liftingColumnIndex = -1;

         // go through nonzero entries of the column
         for( int k = colVector.size() - 1; k >= 0; k-- )
         {
            const Rational& value = colVector.value(k);

            if( spxAbs(value) < minValue )
            {
               MSG_DEBUG( std::cout << "   --> nonzero " << k << " has value " << rationalToString(value) << " in row " << colVector.index(k) << "\n" );

               // add new column equal to maxValue times original column
               if( !addedLiftingRow )
               {
                  MSG_DEBUG( std::cout << "            --> adding lifting row\n" );

                  assert(liftingRowVector.size() == 0);

                  liftingColumnIndex = numColsRational();
                  liftingRowVector.add(i, minValue);
                  liftingRowVector.add(liftingColumnIndex, -1);

                  _rationalLP->addRow(LPRowRational(0, liftingRowVector, 0));
                  _realLP->addRow(LPRowReal(0.0, DSVectorReal(liftingRowVector), 0.0));

                  assert(liftingColumnIndex == numColsRational() - 1);
                  assert(liftingColumnIndex == numCols() - 1);

                  _rationalLP->changeBounds(liftingColumnIndex, _rationalNegInfty, _rationalPosInfty);
                  _realLP->changeBounds(liftingColumnIndex, -realParam(SoPlexBase<R>::INFTY), realParam(SoPlexBase<R>::INFTY));

                  liftingRowVector.clear();
                  addedLiftingRow = true;
               }

               // get row index
               int rowIndex = colVector.index(k);
               assert(rowIndex >= 0);
               assert(rowIndex < _beforeLiftRows);
               assert(liftingColumnIndex == numColsRational() - 1);

               MSG_DEBUG( std::cout << "            --> changing matrix\n" );

               // remove nonzero from original column
               _rationalLP->changeElement(rowIndex, i, 0);
               _realLP->changeElement(rowIndex, i, 0.0);

               // add nonzero divided by maxValue to new column
               Rational newValue(value);
               newValue /= minValue;
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

         _basisStatusCols.append(numColsRational() - _beforeLiftCols, SPxSolverBase<R>::BASIC);
         _basisStatusRows.append(numRowsRational() - _beforeLiftRows, SPxSolverBase<R>::FIXED);
         _rationalLUSolver.clear();
      }

      MSG_DEBUG( _realLP->writeFileLPBase("afterLift.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();

      if( numColsRational() > _beforeLiftCols || numRowsRational() > _beforeLiftRows )
      {
         MSG_INFO1( spxout, spxout << "Added " << numColsRational() - _beforeLiftCols << " columns and "
            << numRowsRational() - _beforeLiftRows << " rows to reduce large matrix coefficients\n." );
      }
   }



   /// undoes lifting
  template <class R>
   void SoPlexBase<R>::_project(SolRational& sol)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("beforeProject.lp", 0, 0, 0) );

      assert(numColsRational() >= _beforeLiftCols);
      assert(numRowsRational() >= _beforeLiftRows);

      // shrink rational LP to original size
      _rationalLP->removeColRange(_beforeLiftCols, numColsRational() - 1);
      _rationalLP->removeRowRange(_beforeLiftRows, numRowsRational() - 1);

      // shrink real LP to original size
      _realLP->removeColRange(_beforeLiftCols, numColsReal() - 1);
      _realLP->removeRowRange(_beforeLiftRows, numRowsReal() - 1);

      // adjust solution
      if( sol.isPrimalFeasible() )
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
      const Rational maxValue = realParam(SoPlexBase<R>::LIFTMAXVAL);

      for( int i = _beforeLiftCols; i < numColsRational() && sol._isDualFeasible; i++ )
      {
         if( spxAbs(maxValue * sol._redCost[i]) > _rationalOpttol )
         {
            MSG_INFO1( spxout, spxout << "Warning: lost dual solution during project phase.\n" );
            sol._isDualFeasible = false;
         }
      }

      if( sol.isDualFeasible() )
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
         if( _basisStatusCols[i] != SPxSolverBase<R>::BASIC )
         {
            MSG_INFO1( spxout, spxout << "Warning: lost basis during project phase because of nonbasic lifting column.\n" );
            _hasBasis = false;
            _rationalLUSolver.clear();
         }
      }

      for( int i = _beforeLiftRows; i < numRowsRational() && _hasBasis; i++ )
      {
         if( _basisStatusRows[i] == SPxSolverBase<R>::BASIC )
         {
            MSG_INFO1( spxout, spxout << "Warning: lost basis during project phase because of basic lifting row.\n" );
            _hasBasis = false;
            _rationalLUSolver.clear();
         }
      }

      if( _hasBasis )
      {
         _basisStatusCols.reSize(_beforeLiftCols);
         _basisStatusRows.reSize(_beforeLiftRows);
         _rationalLUSolver.clear();
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("afterProject.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// stores objective, bounds, and sides of real LP
  template <class R>
   void SoPlexBase<R>::_storeLPReal()
   {
#ifndef SOPLEX_MANUAL_ALT
      if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_MANUAL )
      {
         _manualRealLP = *_realLP;
         return;
      }
#endif

      _manualLower = _realLP->lower();
      _manualUpper = _realLP->upper();
      _manualLhs = _realLP->lhs();
      _manualRhs = _realLP->rhs();
      _manualObj.reDim(_realLP->nCols());
      _realLP->getObj(_manualObj);
   }



   /// restores objective, bounds, and sides of real LP
  template <class R>
   void SoPlexBase<R>::_restoreLPReal()
   {
      if( intParam(SoPlexBase<R>::SYNCMODE) == SYNCMODE_MANUAL )
      {
#ifndef SOPLEX_MANUAL_ALT
         _solver.loadLP(_manualRealLP);
#else
         _realLP->changeLower(_manualLower);
         _realLP->changeUpper(_manualUpper);
         _realLP->changeLhs(_manualLhs);
         _realLP->changeRhs(_manualRhs);
         _realLP->changeObj(_manualObj);
#endif

         if( _hasBasis )
         {
            // in manual sync mode, if the right-hand side of an equality constraint is not floating-point
            // representable, the user might have constructed the constraint in the real LP by rounding down the
            // left-hand side and rounding up the right-hand side; if the basis status is fixed, we need to adjust it
            for( int i = 0; i < _solver.nRows(); i++ )
            {
               if( _basisStatusRows[i] == SPxSolverBase<R>::FIXED && _solver.lhs(i) != _solver.rhs(i) )
               {
                  assert(_solver.rhs(i) == spxNextafter(_solver.lhs(i), infinity));
                  if( _hasSolRational && _solRational.isDualFeasible()
                     && ((intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE && _solRational._dual[i] > 0)
                        || (intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MINIMIZE && _solRational._dual[i] < 0)) )
                  {
                     _basisStatusRows[i] = SPxSolverBase<R>::ON_UPPER;
                  }
                  else
                  {
                     _basisStatusRows[i] = SPxSolverBase<R>::ON_LOWER;
                  }
               }
            }

            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
            _hasBasis = (_solver.basis().status() > SPxBasisBase<R>::NO_PROBLEM);
         }
      }
      else
      {
         _realLP->changeLower(_manualLower);
         _realLP->changeUpper(_manualUpper);
         _realLP->changeLhs(_manualLhs);
         _realLP->changeRhs(_manualRhs);
         _realLP->changeObj(_manualObj);
      }
   }



   /// introduces slack variables to transform inequality constraints into equations for both rational and real LP,
   /// which should be in sync
  template <class R>
   void SoPlexBase<R>::_transformEquality()
   {
      MSG_DEBUG( std::cout << "Transforming rows to equation form.\n" );

      // start timing
      _statistics->transformTime->start();

      MSG_DEBUG( _realLP->writeFileLPBase("beforeTransEqu.lp", 0, 0, 0) );

      // clear array of slack columns
      _slackCols.clear();

      // add artificial slack variables to convert inequality to equality constraints
      for( int i = 0; i < numRowsRational(); i++ )
      {
         assert((lhsRational(i) == rhsRational(i)) == (_rowTypes[i] == RANGETYPE_FIXED));
         if( _rowTypes[i] != RANGETYPE_FIXED )
         {
            _slackCols.add(_rationalZero, -rhsRational(i), *_unitVectorRational(i), -lhsRational(i));
            if( _rationalLP->lhs(i) != 0 )
               _rationalLP->changeLhs(i, _rationalZero);
            if( _rationalLP->rhs(i) != 0 )
               _rationalLP->changeRhs(i, _rationalZero);
            assert(_rationalLP->lhs(i) == 0);
            assert(_rationalLP->rhs(i) == 0);
            _realLP->changeRange(i, 0.0, 0.0);
            _colTypes.append(_switchRangeType(_rowTypes[i]));
            _rowTypes[i] = RANGETYPE_FIXED;
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
            case SPxSolverBase<R>::ON_LOWER:
               _basisStatusCols.append(SPxSolverBase<R>::ON_UPPER);
               break;
            case SPxSolverBase<R>::ON_UPPER:
               _basisStatusCols.append(SPxSolverBase<R>::ON_LOWER);
               break;
            case SPxSolverBase<R>::BASIC:
            case SPxSolverBase<R>::FIXED:
            default:
               _basisStatusCols.append(_basisStatusRows[row]);
               break;
            }

            _basisStatusRows[row] = SPxSolverBase<R>::FIXED;
         }

         _rationalLUSolver.clear();
      }

      MSG_DEBUG( _realLP->writeFileLPBase("afterTransEqu.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();

      if( _slackCols.num() > 0 )
      {
         MSG_INFO1( spxout, spxout << "Added " << _slackCols.num() << " slack columns to transform rows to equality form.\n" );
      }
   }



   /// restores original problem
  template <class R>
   void SoPlexBase<R>::_untransformEquality(SolRational& sol)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("beforeUntransEqu.lp", 0, 0, 0) );

      int numCols = numColsRational();
      int numOrigCols = numColsRational() - _slackCols.num();

      // adjust solution
      if( sol.isPrimalFeasible() )
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

      // adjust basis
      if( _hasBasis )
      {
         for( int i = 0; i < _slackCols.num(); i++ )
         {
            int col = numOrigCols + i;
            int row = _slackCols.colVector(i).index(0);

            assert(row >= 0);
            assert(row < numRowsRational());
            assert(_basisStatusRows[row] != SPxSolverBase<R>::UNDEFINED);
            assert(_basisStatusRows[row] != SPxSolverBase<R>::ZERO || lhsRational(row) == 0);
            assert(_basisStatusRows[row] != SPxSolverBase<R>::ZERO || rhsRational(row) == 0);
            assert(_basisStatusRows[row] != SPxSolverBase<R>::BASIC || _basisStatusCols[col] != SPxSolverBase<R>::BASIC);

            MSG_DEBUG( std::cout << "slack column " << col << " for row " << row
               << ": col status=" << _basisStatusCols[col]
               << ", row status=" << _basisStatusRows[row]
               << ", redcost=" << rationalToString(sol._redCost[col])
               << ", dual=" << rationalToString(sol._dual[row]) << "\n" );

            if( _basisStatusRows[row] != SPxSolverBase<R>::BASIC )
            {
               switch( _basisStatusCols[col] )
               {
               case SPxSolverBase<R>::ON_LOWER:
                  _basisStatusRows[row] = SPxSolverBase<R>::ON_UPPER;
                  break;
               case SPxSolverBase<R>::ON_UPPER:
                  _basisStatusRows[row] = SPxSolverBase<R>::ON_LOWER;
                  break;
               case SPxSolverBase<R>::BASIC:
               case SPxSolverBase<R>::FIXED:
               default:
                  _basisStatusRows[row] = _basisStatusCols[col];
                  break;
               }
            }
         }

         _basisStatusCols.reSize(numOrigCols);
         if( _slackCols.num() > 0 )
            _rationalLUSolver.clear();
      }

      // not earlier because of debug message
      if( sol.isDualFeasible() )
      {
         sol._redCost.reDim(numOrigCols);
      }

      // restore sides and remove slack columns
      for( int i = 0; i < _slackCols.num(); i++ )
      {
         int col = numOrigCols + i;
         int row = _slackCols.colVector(i).index(0);

         if( upperRational(col) != 0 )
            _rationalLP->changeLhs(row, -upperRational(col));
         if( lowerRational(col) != 0 )
            _rationalLP->changeRhs(row, -lowerRational(col));
         assert(_rationalLP->lhs(row) == -upperRational(col));
         assert(_rationalLP->rhs(row) == -lowerRational(col));
         _rowTypes[row] = _switchRangeType(_colTypes[col]);
      }

      _rationalLP->removeColRange(numOrigCols, numCols - 1);
      _realLP->removeColRange(numOrigCols, numCols - 1);
      _colTypes.reSize(numOrigCols);

      // objective, bounds, and sides of real LP are restored only after _solveRational()

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("afterUntransEqu.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// transforms LP to unboundedness problem by moving the objective function to the constraints, changing right-hand
   /// side and bounds to zero, and adding an auxiliary variable for the decrease in the objective function
  template <class R>
   void SoPlexBase<R>::_transformUnbounded()
   {
      MSG_INFO1( spxout, spxout << "Setting up LP to compute primal unbounded ray.\n" );

      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("beforeTransUnbounded.lp", 0, 0, 0) );

      // store bounds
      _unboundedLower.reDim(numColsRational());
      _unboundedUpper.reDim(numColsRational());
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         if( _lowerFinite(_colTypes[c]) )
            _unboundedLower[c] = lowerRational(c);
         if( _upperFinite(_colTypes[c]) )
            _unboundedUpper[c] = upperRational(c);
      }

      // store sides
      _unboundedLhs.reDim(numRowsRational());
      _unboundedRhs.reDim(numRowsRational());
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         if( _lowerFinite(_rowTypes[r]) )
            _unboundedLhs[r] = lhsRational(r);
         if( _upperFinite(_rowTypes[r]) )
            _unboundedRhs[r] = rhsRational(r);
      }

      // make right-hand side zero
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
         if( _lowerFinite(_rowTypes[r]) )
         {
            _rationalLP->changeLhs(r, 0);
            _realLP->changeLhs(r, 0.0);
         }
         else if( _realLP->lhs(r) > -realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeLhs(r, -realParam(SoPlexBase<R>::INFTY));

         assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));
         if( _upperFinite(_rowTypes[r]) )
         {
            _rationalLP->changeRhs(r, 0);
            _realLP->changeRhs(r, 0.0);
         }
         else if( _realLP->rhs(r) < realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeRhs(r, realParam(SoPlexBase<R>::INFTY));
      }

      // transform objective function to constraint and add auxiliary variable
      int numOrigCols = numColsRational();
      DSVectorRational obj(numOrigCols + 1);
      ///@todo implement this without copying the objective function
      obj = _rationalLP->maxObj();
      obj.add(numOrigCols, -1);
      _rationalLP->addRow(LPRowRational(0, obj, 0));
      _realLP->addRow(LPRowReal(0, DSVectorReal(obj), 0));
      _rowTypes.append(RANGETYPE_FIXED);

      assert(numColsRational() == numOrigCols + 1);

      // set objective coefficient and bounds for auxiliary variable
      _rationalLP->changeMaxObj(numOrigCols, 1);
      _realLP->changeMaxObj(numOrigCols, 1.0);

      _rationalLP->changeBounds(numOrigCols, _rationalNegInfty, 1);
      _realLP->changeBounds(numOrigCols, -realParam(SoPlexBase<R>::INFTY), 1.0);
      _colTypes.append(RANGETYPE_UPPER);

      // set objective coefficients to zero and adjust bounds for problem variables
      for( int c = numColsRational() - 2; c >= 0; c-- )
      {
         _rationalLP->changeObj(c, 0);
         _realLP->changeObj(c, 0.0);

         assert((lowerRational(c) > _rationalNegInfty) == _lowerFinite(_colTypes[c]));
         if( _lowerFinite(_colTypes[c]) )
         {
            _rationalLP->changeLower(c, 0);
            _realLP->changeLower(c, 0.0);
         }
         else if( _realLP->lower(c) > -realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeLower(c, -realParam(SoPlexBase<R>::INFTY));

         assert((upperRational(c) < _rationalPosInfty) == _upperFinite(_colTypes[c]));
         if( _upperFinite(_colTypes[c]) )
         {
            _rationalLP->changeUpper(c, 0);
            _realLP->changeUpper(c, 0.0);
         }
         else if( _realLP->upper(c) < realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeUpper(c, realParam(SoPlexBase<R>::INFTY));
      }

      // adjust basis
      if( _hasBasis )
      {
         _basisStatusCols.append(SPxSolverBase<R>::ON_UPPER);
         _basisStatusRows.append(SPxSolverBase<R>::BASIC);
         _rationalLUSolver.clear();
      }

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("afterTransUnbounded.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// undoes transformation to unboundedness problem
  template <class R>
   void SoPlexBase<R>::_untransformUnbounded(SolRational& sol, bool unbounded)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("beforeUntransUnbounded.lp", 0, 0, 0) );

      int numOrigCols = numColsRational() - 1;
      int numOrigRows = numRowsRational() - 1;
      const Rational& tau = sol._primal[numOrigCols];

      // adjust solution and basis
      if( unbounded )
      {
         assert(tau >= _rationalPosone);

         sol._isPrimalFeasible= false;
         sol._hasPrimalRay = true;
         sol._isDualFeasible = false;
         sol._hasDualFarkas = false;

         if( tau != 1 )
            sol._primal /= tau;

         sol._primalRay = sol._primal;
         sol._primalRay.reDim(numOrigCols);

         _hasBasis = (_basisStatusCols[numOrigCols] != SPxSolverBase<R>::BASIC && _basisStatusRows[numOrigRows] == SPxSolverBase<R>::BASIC);
         _basisStatusCols.reSize(numOrigCols);
         _basisStatusRows.reSize(numOrigRows);
      }
      else if( boolParam(SoPlexBase<R>::TESTDUALINF) && tau < _rationalFeastol )
      {
         const Rational& alpha = sol._dual[numOrigRows];

         assert(sol._isDualFeasible);
         assert(alpha <= _rationalFeastol - _rationalPosone);

         sol._isPrimalFeasible= false;
         sol._hasPrimalRay = false;
         sol._hasDualFarkas = false;

         if( alpha != -1 )
         {
            sol._dual /= -alpha;
            sol._redCost /= -alpha;
         }
         sol._dual.reDim(numOrigRows);
         sol._redCost.reDim(numOrigCols);
      }
      else
      {
         sol.invalidate();
         _hasBasis = false;
         _basisStatusCols.reSize(numOrigCols);
         _basisStatusCols.reSize(numOrigRows);
      }

      // recover objective function
      const SVectorRational& objRowVector = _rationalLP->rowVector(numOrigRows);
      for( int i = objRowVector.size() - 1; i >= 0; i-- )
      {
         _rationalLP->changeMaxObj(objRowVector.index(i), objRowVector.value(i));
         _realLP->changeMaxObj(objRowVector.index(i), Real(objRowVector.value(i)));
      }

      // remove objective function constraint and auxiliary variable
      _rationalLP->removeRow(numOrigRows);
      _realLP->removeRow(numOrigRows);
      _rowTypes.reSize(numOrigRows);

      _rationalLP->removeCol(numOrigCols);
      _realLP->removeCol(numOrigCols);
      _colTypes.reSize(numOrigCols);

      // restore objective, sides and bounds
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         if( _lowerFinite(_rowTypes[r]) )
         {
            _rationalLP->changeLhs(r, _unboundedLhs[r]);
            _realLP->changeLhs(r, Real(_unboundedLhs[r]));
         }
         if( _upperFinite(_rowTypes[r]) )
         {
            _rationalLP->changeRhs(r, _unboundedRhs[r]);
            _realLP->changeRhs(r, Real(_unboundedRhs[r]));
         }
         assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
         assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));
         assert((lhsReal(r) > -realParam(SoPlexBase<R>::INFTY)) == _lowerFinite(_rowTypes[r]));
         assert((rhsReal(r) < realParam(SoPlexBase<R>::INFTY)) == _upperFinite(_rowTypes[r]));
      }
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         if( _lowerFinite(_colTypes[c]) )
         {
            _rationalLP->changeLower(c, _unboundedLower[c]);
            _realLP->changeLower(c, Real(_unboundedLower[c]));
         }
         if( _upperFinite(_colTypes[c]) )
         {
            _rationalLP->changeUpper(c, _unboundedUpper[c]);
            _realLP->changeUpper(c, Real(_unboundedUpper[c]));
         }
         assert((lowerRational(c) > _rationalNegInfty) == _lowerFinite(_colTypes[c]));
         assert((upperRational(c) < _rationalPosInfty) == _upperFinite(_colTypes[c]));
         assert((lowerReal(c) > -realParam(SoPlexBase<R>::INFTY)) == _lowerFinite(_colTypes[c]));
         assert((upperReal(c) < realParam(SoPlexBase<R>::INFTY)) == _upperFinite(_colTypes[c]));
      }

      // invalidate rational basis factorization
      _rationalLUSolver.clear();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("afterUntransUnbounded.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// store basis
  template <class R>
   void SoPlexBase<R>::_storeBasis()
   {
      assert(!_storedBasis);

      if( _hasBasis )
      {
         _storedBasis = true;
         _storedBasisStatusCols = _basisStatusCols;
         _storedBasisStatusRows = _basisStatusRows;
      }
      else
         _storedBasis = false;
   }



   /// restore basis
  template <class R>
   void SoPlexBase<R>::_restoreBasis()
   {
      if( _storedBasis )
      {
         _hasBasis = true;
         _basisStatusCols = _storedBasisStatusCols;
         _basisStatusRows = _storedBasisStatusRows;
         _storedBasis = false;
      }
   }



   /// transforms LP to feasibility problem by removing the objective function, shifting variables, and homogenizing the
   /// right-hand side
  template <class R>
   void SoPlexBase<R>::_transformFeasibility()
   {
      MSG_INFO1( spxout, spxout << "Setting up LP to test for feasibility.\n" );

      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("beforeTransFeas.lp", 0, 0, 0) );

      // store objective function
      _feasObj.reDim(numColsRational());
      for( int c = numColsRational() - 1; c >= 0; c-- )
         _feasObj[c] = _rationalLP->maxObj(c);

      // store bounds
      _feasLower.reDim(numColsRational());
      _feasUpper.reDim(numColsRational());
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         if( _lowerFinite(_colTypes[c]) )
            _feasLower[c] = lowerRational(c);
         if( _upperFinite(_colTypes[c]) )
            _feasUpper[c] = upperRational(c);
      }

      // store sides
      _feasLhs.reDim(numRowsRational());
      _feasRhs.reDim(numRowsRational());
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         if( _lowerFinite(_rowTypes[r]) )
            _feasLhs[r] = lhsRational(r);
         if( _upperFinite(_rowTypes[r]) )
            _feasRhs[r] = rhsRational(r);
      }

      // set objective coefficients to zero; shift primal space such as to guarantee that the zero solution is within
      // the bounds
      Rational shiftValue;
      Rational shiftValue2;
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         _rationalLP->changeMaxObj(c, 0);
         _realLP->changeMaxObj(c, 0.0);

         if( lowerRational(c) > 0 )
         {
            const SVectorRational& colVector = colVectorRational(c);

            for( int i = 0; i < colVector.size(); i++ )
            {
               shiftValue = colVector.value(i);
               shiftValue *= lowerRational(c);
               int r = colVector.index(i);

               assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
               assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));

               if( _lowerFinite(_rowTypes[r]) && _upperFinite(_rowTypes[r]) )
               {
                  shiftValue2 = lhsRational(r);
                  shiftValue2 -= shiftValue;
                  _rationalLP->changeLhs(r, shiftValue2);
                  _realLP->changeLhs(r, Real(shiftValue2));

                  shiftValue -= rhsRational(r);
                  shiftValue *= -1;
                  _rationalLP->changeRhs(r, shiftValue);
                  _realLP->changeRhs(r, Real(shiftValue));
               }
               else if( _lowerFinite(_rowTypes[r]) )
               {
                  shiftValue -= lhsRational(r);
                  shiftValue *= -1;
                  _rationalLP->changeLhs(r, shiftValue);
                  _realLP->changeLhs(r, Real(shiftValue));
               }
               else if( _upperFinite(_rowTypes[r]) )
               {
                  shiftValue -= rhsRational(r);
                  shiftValue *= -1;
                  _rationalLP->changeRhs(r, shiftValue);
                  _realLP->changeRhs(r, Real(shiftValue));
               }
            }

            assert((upperRational(c) < _rationalPosInfty) == _upperFinite(_colTypes[c]));
            if( _upperFinite(_colTypes[c]) )
            {
               _rationalLP->changeBounds(c, 0, upperRational(c) - lowerRational(c));
               _realLP->changeBounds(c, 0.0, Real(upperRational(c)));
            }
            else if( _realLP->upper(c) < realParam(SoPlexBase<R>::INFTY) )
            {
               _rationalLP->changeLower(c, 0);
               _realLP->changeBounds(c, 0.0, realParam(SoPlexBase<R>::INFTY));
            }
            else
            {
               _rationalLP->changeLower(c, 0);
               _realLP->changeLower(c, 0.0);
            }
         }
         else if( upperRational(c) < 0 )
         {
            const SVectorRational& colVector = colVectorRational(c);

            for( int i = 0; i < colVector.size(); i++ )
            {
               shiftValue = colVector.value(i);
               shiftValue *= upperRational(c);
               int r = colVector.index(i);

               assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
               assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));

               if( _lowerFinite(_rowTypes[r]) && _upperFinite(_rowTypes[r]) )
               {
                  shiftValue2 = lhsRational(r);
                  shiftValue2 -= shiftValue;
                  _rationalLP->changeLhs(r, shiftValue2);
                  _realLP->changeLhs(r, Real(shiftValue2));

                  shiftValue -= rhsRational(r);
                  shiftValue *= -1;
                  _rationalLP->changeRhs(r, shiftValue);
                  _realLP->changeRhs(r, Real(shiftValue));
               }
               else if( _lowerFinite(_rowTypes[r]) )
               {
                  shiftValue -= lhsRational(r);
                  shiftValue *= -1;
                  _rationalLP->changeLhs(r, shiftValue);
                  _realLP->changeLhs(r, Real(shiftValue));
               }
               else if( _upperFinite(_rowTypes[r]) )
               {
                  shiftValue -= rhsRational(r);
                  shiftValue *= -1;
                  _rationalLP->changeRhs(r, shiftValue);
                  _realLP->changeRhs(r, Real(shiftValue));
               }
            }

            assert((lowerRational(c) > _rationalNegInfty) == _lowerFinite(_colTypes[c]));
            if( _lowerFinite(_colTypes[c]) )
            {
               _rationalLP->changeBounds(c, lowerRational(c) - upperRational(c), 0);
               _realLP->changeBounds(c, Real(lowerRational(c)), 0.0);
            }
            else if( _realLP->lower(c) > -realParam(SoPlexBase<R>::INFTY) )
            {
               _rationalLP->changeUpper(c, 0);
               _realLP->changeBounds(c, -realParam(SoPlexBase<R>::INFTY), 0.0);
            }
            else
            {
               _rationalLP->changeUpper(c, 0);
               _realLP->changeUpper(c, 0.0);
            }
         }
         else
         {
            if( _lowerFinite(_colTypes[c]) )
               _realLP->changeLower(c, Real(lowerRational(c)));
            else if( _realLP->lower(c) > -realParam(SoPlexBase<R>::INFTY) )
               _realLP->changeLower(c, -realParam(SoPlexBase<R>::INFTY));
            if( _upperFinite(_colTypes[c]) )
               _realLP->changeUpper(c, Real(upperRational(c)));
            else if( _realLP->upper(c) < realParam(SoPlexBase<R>::INFTY) )
               _realLP->changeUpper(c, realParam(SoPlexBase<R>::INFTY));
         }

         assert(lowerReal(c) <= upperReal(c));
      }

      // homogenize sides
      _tauColVector.clear();
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         if( lhsRational(r) > 0 )
         {
            _tauColVector.add(r, lhsRational(r));
            assert((rhsRational(r) < _rationalPosInfty) == _upperFinite(_rowTypes[r]));
            if( _upperFinite(_rowTypes[r]) )
            {
               _rationalLP->changeRange(r, 0, rhsRational(r) - lhsRational(r));
               _realLP->changeRange(r, 0.0, Real(rhsRational(r)));
            }
            else
            {
               _rationalLP->changeLhs(r, 0);
               _realLP->changeLhs(r, 0.0);
               if( _realLP->rhs(r) < realParam(SoPlexBase<R>::INFTY) )
                  _realLP->changeRhs(r, realParam(SoPlexBase<R>::INFTY));
            }
         }
         else if( rhsRational(r) < 0 )
         {
            _tauColVector.add(r, rhsRational(r));
            assert((lhsRational(r) > _rationalNegInfty) == _lowerFinite(_rowTypes[r]));
            if( _lowerFinite(_rowTypes[r]) )
            {
               _rationalLP->changeRange(r, lhsRational(r) - rhsRational(r), 0);
               _realLP->changeRange(r, Real(lhsRational(r)), 0.0);
            }
            else
            {
               _rationalLP->changeRhs(r, 0);
               _realLP->changeRhs(r, 0.0);
               if( _realLP->lhs(r) > -realParam(SoPlexBase<R>::INFTY) )
                  _realLP->changeLhs(r, -realParam(SoPlexBase<R>::INFTY));
            }
         }
         else
         {
            if( _lowerFinite(_rowTypes[r]) )
               _realLP->changeLhs(r, Real(lhsRational(r)));
            else if( _realLP->lhs(r) > -realParam(SoPlexBase<R>::INFTY) )
               _realLP->changeLhs(r, -realParam(SoPlexBase<R>::INFTY));
            if( _upperFinite(_rowTypes[r]) )
               _realLP->changeRhs(r, Real(rhsRational(r)));
            else if( _realLP->rhs(r) < realParam(SoPlexBase<R>::INFTY) )
               _realLP->changeRhs(r, realParam(SoPlexBase<R>::INFTY));
         }

         assert(rhsReal(r) <= rhsReal(r));
      }

      ///@todo exploit this case by returning without LP solving
      if( _tauColVector.size() == 0 )
      {
         MSG_INFO3( spxout, spxout << "LP is trivially feasible.\n" );
      }

      // add artificial column
      SPxColId id;
      _tauColVector *= -1;
      _rationalLP->addCol(id,
         LPColRational((intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE ? _rationalPosone : _rationalNegone),
            _tauColVector, 1, 0));
      _realLP->addCol(id,
         LPColReal((intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE ? 1.0 : -1.0),
            DSVectorReal(_tauColVector), 1.0, 0.0));
      _colTypes.append(RANGETYPE_BOXED);

      // adjust basis
      if( _hasBasis )
      {
         _basisStatusCols.append(SPxSolverBase<R>::ON_UPPER);
      }

      // invalidate rational basis factorization
      _rationalLUSolver.clear();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("afterTransFeas.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();
   }



   /// undoes transformation to feasibility problem
  template <class R>
   void SoPlexBase<R>::_untransformFeasibility(SolRational& sol, bool infeasible)
   {
      // start timing
      _statistics->transformTime->start();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("beforeUntransFeas.lp", 0, 0, 0) );

      int numOrigCols = numColsRational() - 1;

      // adjust solution and basis
      if( infeasible )
      {
         assert(sol._isDualFeasible);
         assert(sol._primal[numOrigCols] < 1);

         sol._isPrimalFeasible= false;
         sol._hasPrimalRay = false;
         sol._isDualFeasible = false;
         sol._hasDualFarkas = true;

         sol._dualFarkas = sol._dual;

         _hasBasis = false;
         _basisStatusCols.reSize(numOrigCols);
      }
      else if( sol._isPrimalFeasible)
      {
         assert(sol._primal[numOrigCols] >= 1);

         sol._hasPrimalRay = false;
         sol._isDualFeasible = false;
         sol._hasDualFarkas = false;

         if( sol._primal[numOrigCols] != 1 )
         {
            sol._slacks /= sol._primal[numOrigCols];
            for( int i = 0; i < numOrigCols; i++ )
               sol._primal[i] /= sol._primal[numOrigCols];
            sol._primal[numOrigCols] = 1;
         }

         sol._primal.reDim(numOrigCols);
         sol._slacks -= _rationalLP->colVector(numOrigCols);

         _hasBasis = (_basisStatusCols[numOrigCols] != SPxSolverBase<R>::BASIC);
         _basisStatusCols.reSize(numOrigCols);
      }
      else
      {
         _hasBasis = false;
         _basisStatusCols.reSize(numOrigCols);
      }

      // restore right-hand side
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         assert(rhsRational(r) >= _rationalPosInfty || lhsRational(r) <= _rationalNegInfty
            || _feasLhs[r] - lhsRational(r) == _feasRhs[r] - rhsRational(r));

         if( _lowerFinite(_rowTypes[r]) )
         {
            _rationalLP->changeLhs(r, _feasLhs[r]);
            _realLP->changeLhs(r, Real(_feasLhs[r]));
         }
         else if( _realLP->lhs(r) > -realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeLhs(r, -realParam(SoPlexBase<R>::INFTY));
         assert(_lowerFinite(_rowTypes[r]) == (lhsRational(r) > _rationalNegInfty));
         assert(_lowerFinite(_rowTypes[r]) == (lhsReal(r) > -realParam(SoPlexBase<R>::INFTY)));

         if( _upperFinite(_rowTypes[r]) )
         {
            _rationalLP->changeRhs(r, _feasRhs[r]);
            _realLP->changeRhs(r, Real(_feasRhs[r]));
         }
         else if( _realLP->rhs(r) < realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeRhs(r, realParam(SoPlexBase<R>::INFTY));
         assert(_upperFinite(_rowTypes[r]) == (rhsRational(r) < _rationalPosInfty));
         assert(_upperFinite(_rowTypes[r]) == (rhsReal(r) < realParam(SoPlexBase<R>::INFTY)));

         assert(lhsReal(r) <= rhsReal(r));
      }

      // unshift primal space and restore objective coefficients
      Rational shiftValue;
      for( int c = numOrigCols - 1; c >= 0; c-- )
      {
         bool shifted = (_lowerFinite(_colTypes[c]) && _feasLower[c] > 0) || (_upperFinite(_colTypes[c]) && _feasUpper[c] < 0);
         assert(shifted || !_lowerFinite(_colTypes[c]) || _feasLower[c] == lowerRational(c));
         assert(shifted || !_upperFinite(_colTypes[c]) || _feasUpper[c] == upperRational(c));
#ifdef SOPLEX_WITH_GMP
         assert(upperRational(c) >= _rationalPosInfty || lowerRational(c) <= _rationalNegInfty
            || _feasLower[c] - lowerRational(c) == _feasUpper[c] - upperRational(c));
#endif

         if( shifted )
         {
            if( _lowerFinite(_colTypes[c]) )
            {
               shiftValue = _feasLower[c];
               shiftValue -= lowerRational(c);
            }
            else if( _upperFinite(_colTypes[c]) )
            {
               shiftValue = _feasUpper[c];
               shiftValue -= upperRational(c);
            }
            if( sol._isPrimalFeasible)
            {
               sol._primal[c] += shiftValue;
               sol._slacks.multAdd(shiftValue, _rationalLP->colVector(c));
            }
         }

         if( _lowerFinite(_colTypes[c]) )
         {
            if( shifted )
               _rationalLP->changeLower(c, _feasLower[c]);
            _realLP->changeLower(c, Real(_feasLower[c]));
         }
         else if( _realLP->lower(c) > -realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeLower(c, -realParam(SoPlexBase<R>::INFTY));
         assert(_lowerFinite(_colTypes[c]) == (lowerRational(c) > -_rationalPosInfty));
         assert(_lowerFinite(_colTypes[c]) == (lowerReal(c) > -realParam(SoPlexBase<R>::INFTY)));

         if( _upperFinite(_colTypes[c]) )
         {
            if( shifted )
               _rationalLP->changeUpper(c, _feasUpper[c]);
            _realLP->changeUpper(c, Real(upperRational(c)));
         }
         else if( _realLP->upper(c) < realParam(SoPlexBase<R>::INFTY) )
            _realLP->changeUpper(c, realParam(SoPlexBase<R>::INFTY));
         assert(_upperFinite(_colTypes[c]) == (upperRational(c) < _rationalPosInfty));
         assert(_upperFinite(_colTypes[c]) == (upperReal(c) < realParam(SoPlexBase<R>::INFTY)));

         _rationalLP->changeMaxObj(c, _feasObj[c]);
         _realLP->changeMaxObj(c, Real(_feasObj[c]));

         assert(lowerReal(c) <= upperReal(c));
      }

      // remove last column
      _rationalLP->removeCol(numOrigCols);
      _realLP->removeCol(numOrigCols);
      _colTypes.reSize(numOrigCols);

      // invalidate rational basis factorization
      _rationalLUSolver.clear();

      // print LP if in debug mode
      MSG_DEBUG( _realLP->writeFileLPBase("afterUntransFeas.lp", 0, 0, 0) );

      // stop timing
      _statistics->transformTime->stop();

#ifndef NDEBUG
#ifdef SOPLEX_WITH_GMP
      if( sol._isPrimalFeasible)
      {
         DVectorRational activity(numRowsRational());
         _rationalLP->computePrimalActivity(sol._primal, activity);
         assert(sol._slacks == activity);
      }
#endif
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
  template <class R>
   void SoPlexBase<R>::_computeInfeasBox(SolRational& sol, bool transformed)
   {
      assert(sol.hasDualFarkas());

      const VectorRational& lower = transformed ? _feasLower : lowerRational();
      const VectorRational& upper = transformed ? _feasUpper : upperRational();
      const VectorRational& lhs = transformed ? _feasLhs : lhsRational();
      const VectorRational& rhs = transformed ? _feasRhs : rhsRational();
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

      ///@todo this currently works only if all constraints are equations aggregate rows and sides using the multipliers of the Farkas ray
      for( int r = 0; r < numRows; r++ )
      {
         ytransA += y[r] * _rationalLP->rowVector(r);
         ytransb += y[r] * (y[r] > 0 ? lhs[r] : rhs[r]);
      }

      // if we work on the feasibility problem, we ignore the last column
      if( transformed )
         ytransA.reDim(numCols);

      MSG_DEBUG( std::cout << "ytransb = " << rationalToString(ytransb) << "\n" );

      // if we choose minus ytransb as vector of multipliers for the bound constraints on the variables, we obtain an
      // exactly feasible dual solution for the LP with zero objective function; we aggregate the bounds of the
      // variables accordingly and store its negation in temp
      temp = 0;
      bool isTempFinite = true;
      for( int c = 0; c < numCols && isTempFinite; c++ )
      {
         const Rational& minusRedCost = ytransA[c];

         if( minusRedCost > 0 )
         {
            assert((upper[c] < _rationalPosInfty) == _upperFinite(_colTypes[c]));
            if( _upperFinite(_colTypes[c]) )
               temp.addProduct(minusRedCost, upper[c]);
            else
               isTempFinite = false;
         }
         else if( minusRedCost < 0 )
         {
            assert((lower[c] > _rationalNegInfty) == _lowerFinite(_colTypes[c]));
            if( _lowerFinite(_colTypes[c]) )
               temp.addProduct(minusRedCost, lower[c]);
            else
               isTempFinite = false;
         }
      }

      MSG_DEBUG( std::cout << "max ytransA*[x_l,x_u] = " << (isTempFinite ? rationalToString(temp) : "infinite") << "\n" );

      // ytransb - temp is the increase in the dual objective along the Farkas ray; if this is positive, the dual is
      // unbounded and certifies primal infeasibility
      if( isTempFinite && temp < ytransb )
      {
         MSG_INFO1( spxout, spxout << "Farkas infeasibility proof verified exactly. (1)\n" );
         return;
      }

      // ensure that array of nonzero elements in ytransA is available
      assert(ytransA.isSetup());
      ytransA.setup();

      // if ytransb is negative, try to make it zero by including a positive lower bound or a negative upper bound
      if( ytransb < 0 )
      {
         for( int c = 0; c < numCols; c++ )
         {
            if( lower[c] > 0 )
            {
               ytransA.setValue(c, ytransA[c] - ytransb / lower[c]);
               ytransb = 0;
               break;
            }
            else if( upper[c] < 0 )
            {
               ytransA.setValue(c, ytransA[c] - ytransb / upper[c]);
               ytransb = 0;
               break;
            }
         }
      }

      // if ytransb is still zero then the zero solution is inside the bounds and cannot be cut off by the Farkas
      // constraint; in this case, we cannot compute a Farkas box
      if( ytransb < 0 )
      {
         MSG_INFO1( spxout, spxout << "Approximate Farkas proof to weak.  Could not compute Farkas box. (1)\n" );
         return;
      }

      // compute the one norm of ytransA
      temp = 0;
      const int size = ytransA.size();
      for( int n = 0; n < size; n++ )
         temp += spxAbs(ytransA.value(n));

      // if the one norm is zero then ytransA is zero the Farkas proof should have been verified above
      assert(temp != 0);

      // initialize variables in loop: size of Farkas box B, flag whether B has been increased, and number of current
      // nonzero in ytransA
      Rational B = ytransb / temp;
      bool success = false;
      int n = 0;

      // loop through nonzeros of ytransA
      MSG_DEBUG( std::cout << "B = " << rationalToString(B) << "\n" );
      assert(ytransb >= 0);

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
         assert((upper[colIdx] < _rationalPosInfty) == _upperFinite(_colTypes[colIdx]));
         assert((lower[colIdx] > _rationalNegInfty) == _lowerFinite(_colTypes[colIdx]));
         if( minusRedCost < 0 && lower[colIdx] > -B && _lowerFinite(_colTypes[colIdx]) )
         {
            ytransA.clearNum(n);
            ytransb.subProduct(minusRedCost, lower[colIdx]);
            temp += minusRedCost;

            assert(ytransb >= 0);
            assert(temp >= 0);
            assert(temp == 0 || ytransb / temp > B);

            // if ytransA and ytransb are zero, we have 0^T x >= 0 and cannot compute a Farkas box
            if( temp == 0 && ytransb == 0 )
            {
               MSG_INFO1( spxout, spxout << "Approximate Farkas proof to weak.  Could not compute Farkas box. (2)\n" );
               return;
            }
            // if ytransb is positive and ytransA is zero, we have 0^T x > 0, proving infeasibility
            else if( temp == 0 )
            {
               assert(ytransb > 0);
               MSG_INFO1( spxout, spxout << "Farkas infeasibility proof verified exactly. (2)\n" );
               return;
            }
            else
            {
               B = ytransb / temp;
               MSG_DEBUG( std::cout << "B = " << rationalToString(B) << "\n" );
            }

            success = true;
         }
         // if the multiplier is negative we inspect the upper bound: if it is finite and within the Farkas box, we can
         // increase B by including it in the Farkas proof
         else if( minusRedCost > 0 && upper[colIdx] < B && _upperFinite(_colTypes[colIdx]) )
         {
            ytransA.clearNum(n);
            ytransb.subProduct(minusRedCost, upper[colIdx]);
            temp -= minusRedCost;

            assert(ytransb >= 0);
            assert(temp >= 0);
            assert(temp == 0 || ytransb / temp > B);

            // if ytransA and ytransb are zero, we have 0^T x >= 0 and cannot compute a Farkas box
            if( temp == 0 && ytransb == 0 )
            {
               MSG_INFO1( spxout, spxout << "Approximate Farkas proof to weak.  Could not compute Farkas box. (2)\n" );
               return;
            }
            // if ytransb is positive and ytransA is zero, we have 0^T x > 0, proving infeasibility
            else if( temp == 0 )
            {
               assert(ytransb > 0);
               MSG_INFO1( spxout, spxout << "Farkas infeasibility proof verified exactly. (2)\n" );
               return;
            }
            else
            {
               B = ytransb / temp;
               MSG_DEBUG( std::cout << "B = " << rationalToString(B) << "\n" );
            }

            success = true;
         }
         // the multiplier is zero, we can ignore the bound constraints on this variable
         else if( minusRedCost == 0 )
            ytransA.clearNum(n);
         // currently this bound cannot be used to increase B; we will check it again in the next round, because B might
         // have increased by then
         else
            n++;
      }

      if( B > 0 )
      {
         MSG_INFO1( spxout, spxout << "Computed Farkas box: provably no feasible solutions with components less than "
            << rationalToString(B) << " in absolute value.\n" );
      }
   }






   /// computes rational inverse of basis matrix as defined by _rationalLUSolverBind
  template <class R>
   void SoPlexBase<R>::_computeBasisInverseRational()
   {
      assert(_rationalLUSolver.status() == SLinSolverRational::UNLOADED
         || _rationalLUSolver.status() == SLinSolverRational::TIME);

      const int matrixdim = numRowsRational();
      assert(_rationalLUSolverBind.size() == matrixdim);

      DataArray< const SVectorRational* > matrix(matrixdim);
      _rationalLUSolverBind.reSize(matrixdim);

      for( int i = 0; i < matrixdim; i++ )
      {
         if( _rationalLUSolverBind[i] >= 0 )
         {
            assert(_rationalLUSolverBind[i] < numColsRational());
            matrix[i] = &colVectorRational(_rationalLUSolverBind[i]);
         }
         else
         {
            assert(-1 - _rationalLUSolverBind[i] >= 0);
            assert(-1 - _rationalLUSolverBind[i] < numRowsRational());
            matrix[i] = _unitVectorRational(-1 - _rationalLUSolverBind[i]);
         }
      }

      // load and factorize rational basis matrix
      if( realParam(SoPlexBase<R>::TIMELIMIT) < realParam(SoPlexBase<R>::INFTY) )
         _rationalLUSolver.setTimeLimit(realParam(SoPlexBase<R>::TIMELIMIT) - _statistics->solvingTime->time());
      else
         _rationalLUSolver.setTimeLimit(-1.0);
      _rationalLUSolver.load(matrix.get_ptr(), matrixdim);

      // record statistics
      _statistics->luFactorizationTimeRational += _rationalLUSolver.getFactorTime();
      _statistics->luFactorizationsRational += _rationalLUSolver.getFactorCount();
      _rationalLUSolver.resetCounters();

      if( _rationalLUSolver.status() == SLinSolverRational::TIME )
      {
         MSG_INFO2( spxout, spxout << "Rational factorization hit time limit.\n" );
      }
      else if( _rationalLUSolver.status() != SLinSolverRational::OK )
      {
         MSG_INFO1( spxout, spxout << "Error performing rational LU factorization.\n" );
      }

      return;
   }



   /// factorizes rational basis matrix in column representation
  template <class R>
   void SoPlexBase<R>::_factorizeColumnRational(SolRational& sol, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols, bool& stoppedTime, bool& stoppedIter, bool& error, bool& optimal)
   {
      // start rational solving time
      _statistics->rationalTime->start();

      stoppedTime = false;
      stoppedIter = false;
      error = false;
      optimal = false;

      const bool maximizing = (intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE);
      const int matrixdim = numRowsRational();
      bool loadMatrix = (_rationalLUSolver.status() == SLinSolverRational::UNLOADED
         || _rationalLUSolver.status() == SLinSolverRational::TIME);
      int numBasicRows;

      assert(loadMatrix || matrixdim == _rationalLUSolver.dim());
      assert(loadMatrix || matrixdim == _rationalLUSolverBind.size());
      if( !loadMatrix && (matrixdim != _rationalLUSolver.dim() || matrixdim != _rationalLUSolverBind.size()) )
      {
         MSG_WARNING( spxout, spxout << "Warning: dimensioning error in rational matrix factorization (recovered).\n" );
         loadMatrix = true;
      }

      _workSol._primal.reDim(matrixdim);
      _workSol._slacks.reDim(matrixdim);
      _workSol._dual.reDim(matrixdim);
      _workSol._redCost.reDim(numColsRational() > matrixdim ? numColsRational() : matrixdim);
      if( loadMatrix )
         _rationalLUSolverBind.reSize(matrixdim);

      DVectorRational& basicPrimalRhs = _workSol._slacks;
      DVectorRational& basicDualRhs = _workSol._redCost;
      DVectorRational& basicPrimal = _workSol._primal;
      DVectorRational& basicDual = _workSol._dual;

      Rational violation;
      Rational primalViolation;
      Rational dualViolation;
      bool primalFeasible = false;
      bool dualFeasible = false;

      assert(basisStatusCols.size() == numColsRational());
      assert(basisStatusRows.size() == numRowsRational());

      int j = 0;
      for( int i = 0; i < basisStatusRows.size(); i++ )
      {
         if( basisStatusRows[i] == SPxSolverBase<R>::BASIC && j < matrixdim )
         {
            basicPrimalRhs[i] = 0;
            basicDualRhs[j] = 0;
            if( loadMatrix )
               _rationalLUSolverBind[j] = -1 - i;
            j++;
         }
         else if( basisStatusRows[i] == SPxSolverBase<R>::ON_LOWER )
            basicPrimalRhs[i] = lhsRational(i);
         else if( basisStatusRows[i] == SPxSolverBase<R>::ON_UPPER )
            basicPrimalRhs[i] = rhsRational(i);
         else if( basisStatusRows[i] == SPxSolverBase<R>::ZERO )
            basicPrimalRhs[i] = 0;
         else if( basisStatusRows[i] == SPxSolverBase<R>::FIXED )
         {
            assert(lhsRational(i) == rhsRational(i));
            basicPrimalRhs[i] = lhsRational(i);
         }
         else if( basisStatusRows[i] == SPxSolverBase<R>::UNDEFINED )
         {
            MSG_INFO1( spxout, spxout << "Undefined basis status of row in rational factorization.\n" );
            error = true;
            goto TERMINATE;
         }
         else
         {
            assert(basisStatusRows[i] == SPxSolverBase<R>::BASIC);
            MSG_INFO1( spxout, spxout << "Too many basic rows in rational factorization.\n" );
            error = true;
            goto TERMINATE;
         }
      }
      numBasicRows = j;

      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( basisStatusCols[i] == SPxSolverBase<R>::BASIC && j < matrixdim )
         {
            basicDualRhs[j] = objRational(i);
            if( loadMatrix )
               _rationalLUSolverBind[j] = i;
            j++;
         }
         else if( basisStatusCols[i] == SPxSolverBase<R>::ON_LOWER )
            basicPrimalRhs.multAdd(-lowerRational(i), colVectorRational(i));
         else if( basisStatusCols[i] == SPxSolverBase<R>::ON_UPPER )
            basicPrimalRhs.multAdd(-upperRational(i), colVectorRational(i));
         else if( basisStatusCols[i] == SPxSolverBase<R>::ZERO )
         {}
         else if( basisStatusCols[i] == SPxSolverBase<R>::FIXED )
         {
            assert(lowerRational(i) == upperRational(i));
            basicPrimalRhs.multAdd(-lowerRational(i), colVectorRational(i));
         }
         else if( basisStatusCols[i] == SPxSolverBase<R>::UNDEFINED )
         {
            MSG_INFO1( spxout, spxout << "Undefined basis status of column in rational factorization.\n" );
            error = true;
            goto TERMINATE;
         }
         else
         {
            assert(basisStatusCols[i] == SPxSolverBase<R>::BASIC);
            MSG_INFO1( spxout, spxout << "Too many basic columns in rational factorization.\n" );
            error = true;
            goto TERMINATE;
         }
      }

      if( j != matrixdim )
      {
         MSG_INFO1( spxout, spxout << "Too few basic entries in rational factorization.\n" );
         error = true;
         goto TERMINATE;
      }

      // load and factorize rational basis matrix
      if( loadMatrix )
         _computeBasisInverseRational();

      if( _rationalLUSolver.status() == SLinSolverRational::TIME )
      {
         stoppedTime = true;
         return;
      }
      else if( _rationalLUSolver.status() != SLinSolverRational::OK )
      {
         error = true;
         return;
      }
      assert(_rationalLUSolver.status() == SLinSolverRational::OK);

      // solve for primal solution
      if( realParam(SoPlexBase<R>::TIMELIMIT) < realParam(SoPlexBase<R>::INFTY) )
         _rationalLUSolver.setTimeLimit(realParam(SoPlexBase<R>::TIMELIMIT) - _statistics->solvingTime->time());
      else
         _rationalLUSolver.setTimeLimit(-1.0);
      _rationalLUSolver.solveRight(basicPrimal, basicPrimalRhs);

      // record statistics
      _statistics->luSolveTimeRational += _rationalLUSolver.getSolveTime();
      _rationalLUSolver.resetCounters();

      if( _isSolveStopped(stoppedTime, stoppedIter) )
      {
         MSG_INFO2( spxout, spxout << "Rational factorization hit time limit while solving for primal.\n" );
         return;
      }

      // check bound violation on basic rows and columns
      j = 0;
      primalViolation = 0;
      primalFeasible = true;
      for( int i = 0; i < basisStatusRows.size(); i++ )
      {
         if( basisStatusRows[i] == SPxSolverBase<R>::BASIC )
         {
            assert(j < matrixdim);
            assert(_rationalLUSolverBind[j] == -1 - i);
            violation = lhsRational(i);
            violation += basicPrimal[j];
            if( violation > primalViolation )
            {
               primalFeasible = false;
               primalViolation = violation;
            }
            violation = rhsRational(i);
            violation += basicPrimal[j];
            violation *= -1;
            if( violation > primalViolation )
            {
               primalFeasible = false;
               primalViolation = violation;
            }
            j++;
         }
      }
      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( basisStatusCols[i] == SPxSolverBase<R>::BASIC )
         {
            assert(j < matrixdim);
            assert(_rationalLUSolverBind[j] == i);
            if( basicPrimal[j] < lowerRational(i) )
            {
               violation = lowerRational(i);
               violation -= basicPrimal[j];
               if( violation > primalViolation )
               {
                  primalFeasible = false;
                  primalViolation = violation;
               }
            }
            if( basicPrimal[j] > upperRational(i) )
            {
               violation = basicPrimal[j];
               violation -= upperRational(i);
               if( violation > primalViolation )
               {
                  primalFeasible = false;
                  primalViolation = violation;
               }
            }
            j++;
         }
      }

      if( !primalFeasible )
      {
         MSG_INFO1( spxout, spxout << "Rational solution primal infeasible.\n" );
      }

      // solve for dual solution
      if( realParam(SoPlexBase<R>::TIMELIMIT) < realParam(SoPlexBase<R>::INFTY) )
         _rationalLUSolver.setTimeLimit(realParam(SoPlexBase<R>::TIMELIMIT) - _statistics->solvingTime->time());
      else
         _rationalLUSolver.setTimeLimit(-1.0);
      _rationalLUSolver.solveLeft(basicDual, basicDualRhs);

      // record statistics
      _statistics->luSolveTimeRational += _rationalLUSolver.getSolveTime();
      _rationalLUSolver.resetCounters();

      if( _isSolveStopped(stoppedTime, stoppedIter) )
      {
         MSG_INFO2( spxout, spxout << "Rational factorization hit time limit while solving for dual.\n" );
         return;
      }

      // check dual violation on nonbasic rows
      dualViolation = 0;
      dualFeasible = true;
      for( int i = 0; i < basisStatusRows.size(); i++ )
      {
         if( _rowTypes[i] == RANGETYPE_FIXED
            && (basisStatusRows[i] == SPxSolverBase<R>::ON_LOWER || basisStatusRows[i] == SPxSolverBase<R>::ON_UPPER) )
         {
            assert(lhsRational(i) == rhsRational(i));
            basisStatusRows[i] = SPxSolverBase<R>::FIXED;
         }

         assert(basisStatusRows[i] != SPxSolverBase<R>::BASIC || basicDual[i] == 0);
         if( basisStatusRows[i] == SPxSolverBase<R>::BASIC || basisStatusRows[i] == SPxSolverBase<R>::FIXED )
            continue;
         else if( basicDual[i] < 0 )
         {
            if( ((maximizing && basisStatusRows[i] != SPxSolverBase<R>::ON_LOWER) || (!maximizing && basisStatusRows[i] != SPxSolverBase<R>::ON_UPPER))
               && (basisStatusRows[i] != SPxSolverBase<R>::ZERO || rhsRational(i) != 0) )
            {
               dualFeasible = false;
               violation = -basicDual[i];
               if( violation > dualViolation )
                  dualViolation = violation;
               MSG_DEBUG( spxout << "negative dual multliplier for row " << i
                  << " with dual " << rationalToString(basicDual[i])
                  << " and status " << basisStatusRows[i]
                  << " and [lhs,rhs] = [" << rationalToString(lhsRational(i)) << "," << rationalToString(rhsRational(i)) << "]"
                  << "\n" );
            }
         }
         else if( basicDual[i] > 0 )
         {
            if( ((maximizing && basisStatusRows[i] != SPxSolverBase<R>::ON_UPPER) || (!maximizing && basisStatusRows[i] != SPxSolverBase<R>::ON_LOWER))
               && (basisStatusRows[i] != SPxSolverBase<R>::ZERO || lhsRational(i) == 0) )
            {
               dualFeasible = false;
               if( basicDual[i] > dualViolation )
                  dualViolation = basicDual[i];
               MSG_DEBUG( spxout << "positive dual multliplier for row " << i
                  << " with dual " << rationalToString(basicDual[i])
                  << " and status " << basisStatusRows[i]
                  << " and [lhs,rhs] = [" << rationalToString(lhsRational(i)) << "," << rationalToString(rhsRational(i)) << "]"
                  << "\n" );
            }
         }
      }

      // check reduced cost violation on nonbasic columns
      for( int i = 0; i < basisStatusCols.size(); i++ )
      {
         if( _colTypes[i] == RANGETYPE_FIXED
            && (basisStatusCols[i] == SPxSolverBase<R>::ON_LOWER || basisStatusCols[i] == SPxSolverBase<R>::ON_UPPER) )
         {
            assert(lowerRational(i) == upperRational(i));
            basisStatusCols[i] = SPxSolverBase<R>::FIXED;
         }

#ifdef SOPLEX_WITH_GMP
         assert(basisStatusCols[i] != SPxSolverBase<R>::BASIC || basicDual * colVectorRational(i) == objRational(i));
#endif
         if( basisStatusCols[i] == SPxSolverBase<R>::BASIC || basisStatusCols[i] == SPxSolverBase<R>::FIXED )
            continue;
         else
         {
            _workSol._redCost[i] = basicDual * colVectorRational(i);
            _workSol._redCost[i] -= objRational(i);
            if( _workSol._redCost[i] > 0 )
            {
               if( ((maximizing && basisStatusCols[i] != SPxSolverBase<R>::ON_LOWER) || (!maximizing && basisStatusCols[i] != SPxSolverBase<R>::ON_UPPER))
                  && (basisStatusCols[i] != SPxSolverBase<R>::ZERO || upperRational(i) != 0) )
               {
                  dualFeasible = false;
                  if( _workSol._redCost[i] > dualViolation )
                     dualViolation = _workSol._redCost[i];
               }
               _workSol._redCost[i] *= -1;
            }
            else if( _workSol._redCost[i] < 0 )
            {
               _workSol._redCost[i] *= -1;
               if( ((maximizing && basisStatusCols[i] != SPxSolverBase<R>::ON_UPPER) || (!maximizing && basisStatusCols[i] != SPxSolverBase<R>::ON_LOWER))
                  && (basisStatusCols[i] != SPxSolverBase<R>::ZERO || lowerRational(i) != 0) )
               {
                  dualFeasible = false;
                  if( _workSol._redCost[i] > dualViolation )
                     dualViolation = _workSol._redCost[i];
               }
            }
            else
               _workSol._redCost[i] *= -1;
         }
      }

      if( !dualFeasible )
      {
         MSG_INFO1( spxout, spxout << "Rational solution dual infeasible.\n" );
      }

      // store solution
      optimal = primalFeasible && dualFeasible;
      if( optimal || boolParam(SoPlexBase<R>::RATFACJUMP) )
      {
         _hasBasis = true;
         if( &basisStatusRows != &_basisStatusRows )
            _basisStatusRows = basisStatusRows;
         if( &basisStatusCols != &_basisStatusCols )
            _basisStatusCols = basisStatusCols;

         sol._primal.reDim(numColsRational());
         j = numBasicRows;
         for( int i = 0; i < basisStatusCols.size(); i++ )
         {
            if( basisStatusCols[i] == SPxSolverBase<R>::BASIC )
            {
               assert(j < matrixdim);
               assert(_rationalLUSolverBind[j] == i);
               sol._primal[i] = basicPrimal[j];
               j++;
            }
            else if( basisStatusCols[i] == SPxSolverBase<R>::ON_LOWER )
               sol._primal[i] = lowerRational(i);
            else if( basisStatusCols[i] == SPxSolverBase<R>::ON_UPPER )
               sol._primal[i] = upperRational(i);
            else if( basisStatusCols[i] == SPxSolverBase<R>::ZERO )
               sol._primal[i] = 0;
            else if( basisStatusCols[i] == SPxSolverBase<R>::FIXED )
            {
               assert(lowerRational(i) == upperRational(i));
               sol._primal[i] = lowerRational(i);
            }
            else
            {
               assert(basisStatusCols[i] == SPxSolverBase<R>::UNDEFINED);
               MSG_INFO1( spxout, spxout << "Undefined basis status of column in rational factorization.\n" );
               error = true;
               goto TERMINATE;
            }
         }
         sol._slacks.reDim(numRowsRational());
         _rationalLP->computePrimalActivity(sol._primal, sol._slacks);
         sol._isPrimalFeasible= true;

         sol._dual = basicDual;
         for( int i = 0; i < numColsRational(); i++ )
         {
            if( basisStatusCols[i] == SPxSolverBase<R>::BASIC )
               sol._redCost[i] = 0;
            else if( basisStatusCols[i] == SPxSolverBase<R>::FIXED )
            {
               sol._redCost[i] = basicDual * colVectorRational(i);
               sol._redCost[i] -= objRational(i);
               sol._redCost[i] *= -1;
            }
            else
               sol._redCost[i] = _workSol._redCost[i];
         }
         sol._isDualFeasible  = true;
      }

   TERMINATE:
      // stop rational solving time
      _statistics->rationalTime->stop();
      return;
   }

   /// attempts rational reconstruction of primal-dual solution
  template <class R>
   bool SoPlexBase<R>::_reconstructSolutionRational(SolRational& sol, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusRows, DataArray< typename SPxSolverBase<R>::VarStatus >& basisStatusCols, const Rational& denomBoundSquared)
   {
      bool success;
      bool isSolBasic;
      DIdxSet basicIndices(numColsRational());

      success = false;
      isSolBasic = true;

      if( !sol.isPrimalFeasible() || !sol.isDualFeasible() )
         return success;

      // start timing and increment statistics counter
      _statistics->reconstructionTime->start();
      _statistics->rationalReconstructions++;

      // reconstruct primal vector
      _workSol._primal = sol._primal;
      for( int j = 0; j < numColsRational(); ++j )
      {
         if( basisStatusCols[j] == SPxSolverBase<R>::BASIC )
            basicIndices.addIdx(j);
      }

      success = reconstructVector(_workSol._primal, denomBoundSquared, &basicIndices);

      if( !success )
      {
         MSG_INFO1( spxout, spxout << "Rational reconstruction of primal solution failed.\n" );
         _statistics->reconstructionTime->stop();
         return success;
      }

      MSG_DEBUG( spxout << "Rational reconstruction of primal solution successful.\n" );

      // check violation of bounds
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         // we want to notify the user whether the reconstructed solution is basic; otherwise, this would be redundant
         typename SPxSolverBase<R>::VarStatus& basisStatusCol = _basisStatusCols[c];
         if( (basisStatusCol == SPxSolverBase<R>::FIXED && _workSol._primal[c] != lowerRational(c))
            || (basisStatusCol == SPxSolverBase<R>::ON_LOWER && _workSol._primal[c] != lowerRational(c))
            || (basisStatusCol == SPxSolverBase<R>::ON_UPPER && _workSol._primal[c] != upperRational(c))
            || (basisStatusCol == SPxSolverBase<R>::ZERO && _workSol._primal[c] != 0)
            || (basisStatusCol == SPxSolverBase<R>::UNDEFINED) )
         {
            isSolBasic = false;
         }

         if( _lowerFinite(_colTypes[c]) && _workSol._primal[c] < lowerRational(c) )
         {
            MSG_DEBUG( std::cout << "Lower bound of variable " << c << " violated by " << rationalToString(lowerRational(c) - _workSol._primal[c]) << "\n" );
            MSG_INFO1( spxout, spxout << "Reconstructed solution primal infeasible (1).\n" );
            _statistics->reconstructionTime->stop();
            return false;
         }

         if( _upperFinite(_colTypes[c]) && _workSol._primal[c] > upperRational(c) )
         {
            MSG_DEBUG( std::cout << "Upper bound of variable " << c << " violated by " << rationalToString(_workSol._primal[c] - upperRational(c)) << "\n" );
            MSG_INFO1( spxout, spxout << "Reconstructed solution primal infeasible (2).\n" );
            _statistics->reconstructionTime->stop();
            return false;
         }
      }

      // compute slacks
      ///@todo we should compute them one by one so we can abort when encountering an infeasibility
      _workSol._slacks.reDim(numRowsRational());
      _rationalLP->computePrimalActivity(_workSol._primal, _workSol._slacks);

      // check violation of sides
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         // we want to notify the user whether the reconstructed solution is basic; otherwise, this would be redundant
         typename SPxSolverBase<R>::VarStatus& basisStatusRow = _basisStatusRows[r];
         if( (basisStatusRow == SPxSolverBase<R>::FIXED && _workSol._slacks[r] != lhsRational(r))
            || (basisStatusRow == SPxSolverBase<R>::ON_LOWER && _workSol._slacks[r] != lhsRational(r))
            || (basisStatusRow == SPxSolverBase<R>::ON_UPPER && _workSol._slacks[r] != rhsRational(r))
            || (basisStatusRow == SPxSolverBase<R>::ZERO && _workSol._slacks[r] != 0)
            || (basisStatusRow == SPxSolverBase<R>::UNDEFINED) )
         {
            isSolBasic = false;
         }

         if( _lowerFinite(_rowTypes[r]) && _workSol._slacks[r] < lhsRational(r) )
         {
            MSG_DEBUG( std::cout << "Lhs of row " << r << " violated by " << rationalToString(lhsRational(r) - _workSol._slacks[r]) << "\n" );
            MSG_INFO1( spxout, spxout << "Reconstructed solution primal infeasible (3).\n" );
            _statistics->reconstructionTime->stop();
            return false;
         }

         if( _upperFinite(_rowTypes[r]) && _workSol._slacks[r] > rhsRational(r) )
         {
            MSG_DEBUG( std::cout << "Rhs of row " << r << " violated by " << rationalToString(_workSol._slacks[r] - rhsRational(r)) << "\n" );
            MSG_INFO1( spxout, spxout << "Reconstructed solution primal infeasible (4).\n" );
            _statistics->reconstructionTime->stop();
            return false;
         }
      }

      // reconstruct dual vector
      _workSol._dual = sol._dual;

      success = reconstructVector(_workSol._dual, denomBoundSquared);

      if( !success )
      {
         MSG_INFO1( spxout, spxout << "Rational reconstruction of dual solution failed.\n" );
         _statistics->reconstructionTime->stop();
         return success;
      }

      MSG_DEBUG( spxout << "Rational reconstruction of dual vector successful.\n" );

      // check dual multipliers before reduced costs because this check is faster since it does not require the
      // computation of reduced costs
      const bool maximizing = (intParam(SoPlexBase<R>::OBJSENSE) == SoPlexBase<R>::OBJSENSE_MAXIMIZE);
      for( int r = numRowsRational() - 1; r >= 0; r-- )
      {
         int sig = sign(_workSol._dual[r]);

         if( (!maximizing && sig > 0) || (maximizing && sig < 0) )
         {
            if( !_lowerFinite(_rowTypes[r]) || _workSol._slacks[r] > lhsRational(r) )
            {
               MSG_DEBUG( std::cout << "complementary slackness violated by row " << r
                  << " with dual " << rationalToString(_workSol._dual[r])
                  << " and slack " << rationalToString(_workSol._slacks[r])
                  << " not at lhs " << rationalToString(lhsRational(r))
                  << "\n" );
               MSG_INFO1( spxout, spxout << "Reconstructed solution dual infeasible (1).\n" );
               _statistics->reconstructionTime->stop();
               return false;
            }

            if( _basisStatusRows[r] != SPxSolverBase<R>::ON_LOWER && _basisStatusRows[r] != SPxSolverBase<R>::FIXED )
            {
               if( _basisStatusRows[r] == SPxSolverBase<R>::BASIC || _basisStatusRows[r] == SPxSolverBase<R>::UNDEFINED )
                  isSolBasic = false;
               else
                  _basisStatusRows[r] = SPxSolverBase<R>::ON_LOWER;
            }
         }
         else if( (!maximizing && sig < 0) || (maximizing && sig > 0) )
         {
            if( !_upperFinite(_rowTypes[r]) || _workSol._slacks[r] < rhsRational(r) )
            {
               MSG_DEBUG( std::cout << "complementary slackness violated by row " << r
                  << " with dual " << rationalToString(_workSol._dual[r])
                  << " and slack " << rationalToString(_workSol._slacks[r])
                  << " not at rhs " << rationalToString(rhsRational(r))
                  << "\n" );
               MSG_INFO1( spxout, spxout << "Reconstructed solution dual infeasible (2).\n" );
               _statistics->reconstructionTime->stop();
               return false;
            }

            if( _basisStatusRows[r] != SPxSolverBase<R>::ON_UPPER && _basisStatusRows[r] != SPxSolverBase<R>::FIXED )
            {
               if( _basisStatusRows[r] == SPxSolverBase<R>::BASIC || _basisStatusRows[r] == SPxSolverBase<R>::UNDEFINED )
                  isSolBasic = false;
               else
                  _basisStatusRows[r] = SPxSolverBase<R>::ON_UPPER;
            }
         }
      }

      // compute reduced cost vector; we assume that the objective function vector has less nonzeros than the reduced
      // cost vector, and so multiplying with -1 first and subtracting the dual activity should be faster than adding
      // the dual activity and negating afterwards
      ///@todo we should compute them one by one so we can abort when encountering an infeasibility
      _workSol._redCost.reDim(numColsRational());
      _rationalLP->getObj(_workSol._redCost);
      _rationalLP->subDualActivity(_workSol._dual, _workSol._redCost);

      // check reduced cost violation
      for( int c = numColsRational() - 1; c >= 0; c-- )
      {
         int sig = sign(_workSol._redCost[c]);

         if( (!maximizing && sig > 0) || (maximizing && sig < 0) )
         {
            if( !_lowerFinite(_colTypes[c]) || _workSol._primal[c] > lowerRational(c) )
            {
               MSG_DEBUG( std::cout << "complementary slackness violated by column " << c
                  << " with reduced cost " << rationalToString(_workSol._redCost[c])
                  << " and value " << rationalToString(_workSol._primal[c])
                  << " not at lower bound " << rationalToString(lowerRational(c))
                  << "\n" );
               MSG_INFO1( spxout, spxout << "Reconstructed solution dual infeasible (3).\n" );
               _statistics->reconstructionTime->stop();
               return false;
            }

            if( _basisStatusCols[c] != SPxSolverBase<R>::ON_LOWER && _basisStatusCols[c] != SPxSolverBase<R>::FIXED )
            {
               if( _basisStatusCols[c] == SPxSolverBase<R>::BASIC || _basisStatusCols[c] == SPxSolverBase<R>::UNDEFINED )
                  isSolBasic = false;
               else
                  _basisStatusCols[c] = SPxSolverBase<R>::ON_LOWER;
            }
         }
         else if( (!maximizing && sig < 0) || (maximizing && sig > 0) )
         {
            if( !_upperFinite(_colTypes[c]) || _workSol._primal[c] < upperRational(c) )
            {
               MSG_DEBUG( std::cout << "complementary slackness violated by column " << c
                  << " with reduced cost " << rationalToString(_workSol._redCost[c])
                  << " and value " << rationalToString(_workSol._primal[c])
                  << " not at upper bound " << rationalToString(upperRational(c))
                  << "\n" );
               MSG_INFO1( spxout, spxout << "Reconstructed solution dual infeasible (4).\n" );
               _statistics->reconstructionTime->stop();
               return false;
            }

            if( _basisStatusCols[c] != SPxSolverBase<R>::ON_UPPER && _basisStatusCols[c] != SPxSolverBase<R>::FIXED )
            {
               if( _basisStatusCols[c] == SPxSolverBase<R>::BASIC || _basisStatusCols[c] == SPxSolverBase<R>::UNDEFINED )
                  isSolBasic = false;
               else
                  _basisStatusCols[c] = SPxSolverBase<R>::ON_UPPER;
            }
         }
      }

      // update solution
      sol._primal = _workSol._primal;
      sol._slacks = _workSol._slacks;
      sol._dual = _workSol._dual;
      sol._redCost = _workSol._redCost;

      if( !isSolBasic )
      {
         MSG_WARNING( spxout, spxout << "Warning: Reconstructed solution not basic.\n" );
         _hasBasis = false;
      }

      // stop timing
      _statistics->reconstructionTime->stop();

      return success;
   }
} // namespace soplex
