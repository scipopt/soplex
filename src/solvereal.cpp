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
   void SoPlex2::_solveReal()
   {
      // start timing
      _statistics->solvingTime.start();

      try
      {
         _preprocessAndSolveReal(!_hasBasis);
         _storeSolutionReal();
      }
      catch( ... )
      {
         MSG_ERROR( spxout << "Exception thrown while solving real LP.\n" );
         _status = SPxSolver::ERROR;
      }

      // stop timing
      _statistics->solvingTime.stop();
   }



   /// checks result of the solving process and solves again without preprocessing if necessary
   void SoPlex2::_evaluateSolutionReal(SPxSimplifier::Result simplificationStatus)
   {
      if( simplificationStatus == SPxSimplifier::INFEASIBLE )
         _status = SPxSolver::INFEASIBLE;
      else if( simplificationStatus == SPxSimplifier::DUAL_INFEASIBLE )
         _status = SPxSolver::INForUNBD;
      else if( simplificationStatus == SPxSimplifier::UNBOUNDED )
         _status = SPxSolver::UNBOUNDED;
      else if( simplificationStatus == SPxSimplifier::VANISHED )
         _status = SPxSolver::OPTIMAL;
      else if( simplificationStatus == SPxSimplifier::OKAY )
         _status = _solver.status();

      // process result
      switch( _status )
      {
      case SPxSolver::OPTIMAL:
         if( !_isRealLPLoaded )
         {
            _resolveWithoutPreprocessing(simplificationStatus);
            return;
         }
         else
            _hasBasis = true;
         break;

      case SPxSolver::UNBOUNDED:
      case SPxSolver::INFEASIBLE:
      case SPxSolver::INForUNBD:
         // in case of infeasibility or unboundedness, we currently can not unsimplify, but have to solve the original LP again
         if( !_isRealLPLoaded )
         {
            _preprocessAndSolveReal(false);
            return;
         }
         else
            _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
         break;

      case SPxSolver::SINGULAR:
         // if preprocessing was applied, try to run again without to avoid singularity
         if( !_isRealLPLoaded )
         {
            _preprocessAndSolveReal(false);
            return;
         }
         // if there was a regular starting basis and the original problem is in the solver, load the basis
         else if( _hasBasis )
         {
            assert(_basisStatusRows.size() == _solver.nRows());
            assert(_basisStatusCols.size() == _solver.nCols());
            _solver.setBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
         }
         break;

      case SPxSolver::ABORT_CYCLING:
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::REGULAR:
      case SPxSolver::RUNNING:
         // store regular basis if there is no simplifier and the original problem is not in the solver because of
         // scaling; non-optimal bases should currently not be unsimplified
         if( _simplifier == 0 && !_isRealLPLoaded && _solver.basis().status() > SPxBasis::NO_PROBLEM )
         {
            _basisStatusRows.reSize(numRowsReal());
            _basisStatusCols.reSize(numColsReal());
            assert(_basisStatusRows.size() == _solver.nRows());
            assert(_basisStatusCols.size() == _solver.nCols());

            _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
            _hasBasis = true;
         }
         else
            _hasBasis = false;
         break;

      default:
         _hasBasis = false;
         break;
      }
   }



   /// solves real LP with/without preprocessing
   void SoPlex2::_preprocessAndSolveReal(bool applyPreprocessing)
   {
      if( applyPreprocessing )
      {
         _enableSimplifierAndScaler();
         _solver.setTerminationValue(realParam(SoPlex2::INFTY));
      }
      else
      {
         _disableSimplifierAndScaler();
         ///@todo implement for both objective senses
         _solver.setTerminationValue(intParam(SoPlex2::OBJSENSE) == SoPlex2::OBJSENSE_MINIMIZE
            ? realParam(SoPlex2::OBJLIMIT_UPPER) : realParam(SoPlex2::INFTY));
      }

      applyPreprocessing = (_simplifier != 0 || _scaler != 0);

      if( _isRealLPLoaded )
      {
         assert(_realLP == &_solver);

         // preprocessing is always applied to the LP in the solver; hence we have to create a copy of the original LP
         // if preprocessing is turned on
         if( applyPreprocessing )
         {
            _realLP = 0;
            spx_alloc(_realLP);
            _realLP = new (_realLP) SPxLPReal(_solver);
            _isRealLPLoaded = false;
         }
      }
      else
      {
         assert(_realLP != &_solver);

         // ensure that the solver has the original problem
         _solver.loadLP(*_realLP);

         // load basis if available
         if( _hasBasis )
         {
            assert(_basisStatusRows.size() == numRowsReal());
            assert(_basisStatusCols.size() == numColsReal());

            ///@todo this should not fail even if the basis is invalid (wrong dimension or wrong number of basic
            ///      entries); fix either in SPxSolver or in SPxBasis
            _solver.setBasis(_basisStatusRows.get_const_ptr(), _basisStatusCols.get_const_ptr());
         }

         // if there is no preprocessing, then the original and the transformed problem are identical and it is more
         // memory-efficient to keep only the problem in the solver
         if( !applyPreprocessing )
         {
            _realLP->~SPxLPReal();
            spx_free(_realLP);
            _realLP = &_solver;
            _isRealLPLoaded = true;
         }
      }

      // assert that we have two problems if and only if we apply preprocessing
      assert(_realLP == &_solver || applyPreprocessing);
      assert(_realLP != &_solver || !applyPreprocessing);

      // apply problem simplification
      SPxSimplifier::Result simplificationStatus = SPxSimplifier::OKAY;
      if( _simplifier != 0 )
      {
         simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlex2::EPSILON_ZERO), Real(rationalParam(SoPlex2::FEASTOL)), Real(rationalParam(SoPlex2::OPTTOL)));
      }

      // run the simplex method if problem has not been solved by the simplifier
      if( simplificationStatus == SPxSimplifier::OKAY )
      {
         if( _scaler != 0 )
            _scaler->scale(_solver);

         _solveRealLPAndRecordStatistics();
      }

      // check the result and run again without preprocessing if necessary
      _evaluateSolutionReal(simplificationStatus);
   }



   /// loads original problem into solver and solves again after it has been solved to optimality with preprocessing
   void SoPlex2::_resolveWithoutPreprocessing(SPxSimplifier::Result simplificationStatus)
   {
      assert(!_isRealLPLoaded);
      assert(_simplifier != 0 || _scaler != 0);
      assert(simplificationStatus == SPxSimplifier::VANISHED
         || (simplificationStatus == SPxSimplifier::OKAY && _solver.status() == SPxSolver::OPTIMAL));

      // if simplifier is active and LP is solved in presolving or to optimality, then we unsimplify to get the basis
      if( _simplifier != 0 )
      {
         assert(!_simplifier->isUnsimplified());

         bool vanished = simplificationStatus == SPxSimplifier::VANISHED;

         // get solution vectors for transformed problem
         DVectorReal primal(vanished ? 0 : _solver.nCols());
         DVectorReal slacks(vanished ? 0 : _solver.nRows());
         DVectorReal dual(vanished ? 0 : _solver.nRows());
         DVectorReal redCost(vanished ? 0 : _solver.nCols());

         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());
         assert(vanished || _basisStatusRows.size() >= _solver.nRows());
         assert(vanished || _basisStatusCols.size() >= _solver.nCols());

         if( !vanished )
         {
            assert(_solver.status() == SPxSolver::OPTIMAL);

            _solver.getPrimal(primal);
            _solver.getSlacks(slacks);
            _solver.getDual(dual);
            _solver.getRedCost(redCost);

            // unscale vectors
            if( _scaler != 0 )
            {
               _scaler->unscalePrimal(primal);
               _scaler->unscaleSlacks(slacks);
               _scaler->unscaleDual(dual);
               _scaler->unscaleRedCost(redCost);
            }

            // get basis of transformed problem
            _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
         }

         try
         {
            _simplifier->unsimplify(primal, dual, slacks, redCost, _basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());

            // store basis for original problem
            _simplifier->getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
            _hasBasis = true;
         }
         catch( ... )
         {
            MSG_INFO1( spxout << "Exception thrown during unsimplification. Resolving without simplifier and scaler.\n" );
         }
      }
      // if the original problem is not in the solver because of scaling, we also need to store the basis
      else if( _scaler != 0 )
      {
         _basisStatusRows.reSize(numRowsReal());
         _basisStatusCols.reSize(numColsReal());
         assert(_basisStatusRows.size() == _solver.nRows());
         assert(_basisStatusCols.size() == _solver.nCols());

         _solver.getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
         _hasBasis = true;
      }

      // resolve the original problem
      _preprocessAndSolveReal(false);
      return;
   }



   /// stores solution of the real LP; before calling this, the real LP must be loaded in the solver and solved (again)
   void SoPlex2::_storeSolutionReal()
   {
      assert(_isRealLPLoaded);

      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::ERROR);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_RATIOTESTER);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_PRICER);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_SOLVER);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NOT_INIT);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::SINGULAR);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::NO_PROBLEM);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::UNBOUNDED);
      assert(_solver.basis().status() != SPxBasis::PRIMAL || status() != SPxSolver::INFEASIBLE);
      assert(_solver.basis().status() != SPxBasis::UNBOUNDED || status() == SPxSolver::UNBOUNDED);
      assert(_solver.basis().status() == SPxBasis::UNBOUNDED || _solver.basis().status() == SPxBasis::NO_PROBLEM || status() != SPxSolver::UNBOUNDED);

      _solReal._hasPrimal = (status() == SPxSolver::OPTIMAL
         || ((_solver.basis().status() == SPxBasis::PRIMAL || _solver.basis().status() == SPxBasis::UNBOUNDED)
            && _solver.shift() < 10.0 * realParam(SoPlex2::EPSILON_ZERO)));
      if( _solReal._hasPrimal )
      {
         _solReal._primal.reDim(_solver.nCols());
         _solReal._slacks.reDim(_solver.nRows());
         _solver.getPrimal(_solReal._primal);
         _solver.getSlacks(_solReal._slacks);
      }

      _solReal._hasPrimalRay = (status() == SPxSolver::UNBOUNDED);
      if( _solReal._hasPrimalRay )
      {
         _solReal._primalRay.reDim(_solver.nCols());
         _solver.getPrimalray(_solReal._primalRay);
      }

      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::ERROR);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_RATIOTESTER);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_PRICER);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_SOLVER);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NOT_INIT);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::SINGULAR);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::NO_PROBLEM);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::UNBOUNDED);
      assert(_solver.basis().status() != SPxBasis::DUAL || status() != SPxSolver::INFEASIBLE);
      assert(_solver.basis().status() != SPxBasis::INFEASIBLE || status() == SPxSolver::INFEASIBLE);
      assert(_solver.basis().status() == SPxBasis::INFEASIBLE || _solver.basis().status() == SPxBasis::NO_PROBLEM || status() != SPxSolver::INFEASIBLE);

      _solReal._hasDual = (status() == SPxSolver::OPTIMAL
         || ((_solver.basis().status() == SPxBasis::DUAL || _solver.basis().status() == SPxBasis::INFEASIBLE)
            && _solver.shift() < 10.0 * realParam(SoPlex2::EPSILON_ZERO)));
      if( _solReal._hasDual )
      {
         _solReal._dual.reDim(_solver.nRows());
         _solReal._redCost.reDim(_solver.nCols());
         _solver.getDual(_solReal._dual);
         _solver.getRedCost(_solReal._redCost);
      }

      _solReal._hasDualFarkas = (status() == SPxSolver::INFEASIBLE && _simplifier == 0);
      if( _solReal._hasDualFarkas )
      {
         _solReal._dualFarkas.reDim(_solver.nRows());
         _solver.getDualfarkas(_solReal._dualFarkas);
      }

      _hasSolReal = true;
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