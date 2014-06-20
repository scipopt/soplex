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

/* This file contains the private functions for the Improved Dual Simplex (IDS)
 * 
 * An important note about the implementation of the IDS is the reliance on the row representation of the basis matrix.
 * The forming of the reduced and complementary problems is involves identifying rows in the row-form of the basis
 * matrix that have zero dual multipliers.
 *
 * Ideally, this work will be extended such that the IDS can be executed using the column-form of the basis matrix. */

namespace soplex
{
   /// solves LP using the improved dual simplex
   void SoPlex::_solveImprovedDualSimplex()
   {
      assert(_hasBasis == true);
      assert(_solver.rep() == SPxSolver::ROW);

      // start timing
      _statistics->solvingTime.start();

      //@todo Need to check this. Where will the solving of the IDS take place? In the SoPlex::solve, _solvereal or
      //_solverational?
      // remember that last solve was in floating-point
      _lastSolveMode = SOLVEMODE_REAL;

      // creating the initial reduced problem from the basis information
      _formIdsReducedProblem();


      // stop timing
      _statistics->solvingTime.stop();
   }

   /// forms the reduced problem
   void SoPlex::_formIdsReducedProblem()
   {
      int* bind = 0;
      int* nonposind = 0;
      int* compatind = 0;
      int nnonposind = 0;
      int ncompatind = 0;


      // the improved dual simplex modifies the LP to form a reduced problem and complementary problem. To achieve this
      // it is necessary to copy the current LP and then perform the modifications on the copy.
      // copying the current LP information
      _idsLP = 0;
      spx_alloc(_idsLP);
      _idsLP = new (_idsLP) SPxLPIds(_solver);


      // get thie indices of the rows with positive dual multipliers and columns with positive reduced costs.
      spx_alloc(bind, numColsReal());
      spx_alloc(nonposind, numColsReal());
      _getNonPositiveDualMultiplierInds(nonposind, bind, &nnonposind);


      // get the compatible columns from the constraint matrix w.r.t the current basis matrix
      spx_alloc(compatind, numRowsReal());
      _getCompatibleColumns(nonposind, compatind, nnonposind, &ncompatind);


   }

   void SoPlex::_getNonPositiveDualMultiplierInds(int* nonposind, int* bind, int* nnonposind) const
   {
      nnonposind = 0;

      // iterating over all columns in the basis matrix
      // this identifies the basis indices and the indicies that are positive.
      for( int i = 0; i < numColsReal(); ++i )
      {
         if( baseId(i).isSPxRowId() )
         {
            bind[i] = -1 - number(SPxRowId(baseId(i)));

            //@todo need to check this regarding min and max problems
            if( fVec()[i] <= 0 )
            {
               nonposind[k] = bind[i];
               nnonposind++;
            }
         }
         else if( baseId(i).isSPxColId() )
         {
            bind[i] = number(SPxColId(baseId(i)));

            if (spxSense() == SPxLP::MINIMIZE)
            {
               //@todo need to check this regarding min and max problems
               if( -fVec()[i] <= 0 )
               {
                  nonposind[k] = bind[i];
                  nnonposind++;
               }
            }
            else
            {
               //@todo need to check this regarding min and max problems
               if( fVec()[i] <= 0 )
               {
                  nonposind[k] = bind[i];
                  nnonposind++;
               }
            }
         }
      }
   }


   /// retrieves the compatible columns from the constraint matrix
   void SoPlex::_getCompatibleColumns(int* nonposind, int* compatind, int nnonposind, int* ncompatind)
   {
      bool compatible;
      SSVector y(numColsReal());

      ncompatind  = 0;

      for( int i = 0; i < numRowsReal(); ++i )
      {
         // the rhs of this calculation are the rows of the constraint matrix
         // so we are solving y B = A_{i,.}
         try
         {
            _solver.basis().solve(y, _solver.thevectors[i]);
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> while computing compatability.\n" );
            return false;
         }

         compatible = true;
         // a compatible row is given by zeros in all columns related to the nonpositive indices
         for( int j = 0; j < nnonposind; ++j )
         {
            if( y[nnonposind[j]] != 0 )
            {
               compatible = false;
               break;
            }
         }

         if( compatible )
         {
            compatid[k] = i;
            ncompatind++;
         }
      }

   }


   /// checks result of the solving process and solves again without preprocessing if necessary
   void SoPlex::_evaluateSolutionReal(SPxSimplifier::Result simplificationStatus)
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
            MSG_INFO1( spxout << " --- transforming basis into original space" << std::endl; )
            _solver.changeObjOffset(0.0);
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
            _solver.changeObjOffset(0.0);
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
            _solver.changeObjOffset(0.0);
            _preprocessAndSolveReal(false);
            return;
         }
         break;

      case SPxSolver::ABORT_CYCLING:
         // if preprocessing was applied, try to run again without to avoid cycling
         if( !_isRealLPLoaded )
         {
            _solver.changeObjOffset(0.0);
            _preprocessAndSolveReal(false);
            return;
         }
         // else fallthrough
      case SPxSolver::ABORT_TIME:
      case SPxSolver::ABORT_ITER:
      case SPxSolver::ABORT_VALUE:
      case SPxSolver::REGULAR:
      case SPxSolver::RUNNING:
         if( _isRealLPLoaded )
            _hasBasis = (_solver.basis().status() > SPxBasis::NO_PROBLEM);
         // store regular basis if there is no simplifier and the original problem is not in the solver because of
         // scaling; non-optimal bases should currently not be unsimplified
         else if( _simplifier == 0 && _solver.basis().status() > SPxBasis::NO_PROBLEM )
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
   void SoPlex::_preprocessAndSolveReal(bool applyPreprocessing)
   {
      _statistics->preprocessingTime.start();

      if( applyPreprocessing )
      {
         _enableSimplifierAndScaler();
         _solver.setTerminationValue(realParam(SoPlex::INFTY));
      }
      else
      {
         _disableSimplifierAndScaler();
         ///@todo implement for both objective senses
         _solver.setTerminationValue(intParam(SoPlex::OBJSENSE) == SoPlex::OBJSENSE_MINIMIZE
            ? realParam(SoPlex::OBJLIMIT_UPPER) : realParam(SoPlex::INFTY));
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
         simplificationStatus = _simplifier->simplify(_solver, realParam(SoPlex::EPSILON_ZERO), realParam(SoPlex::FEASTOL), realParam(SoPlex::OPTTOL));
         _solver.changeObjOffset(_simplifier->getObjoffset());
      }

      _statistics->preprocessingTime.stop();

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
   void SoPlex::_resolveWithoutPreprocessing(SPxSimplifier::Result simplificationStatus)
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
            _simplifier->getBasis(_basisStatusRows.get_ptr(), _basisStatusCols.get_ptr());
            _hasBasis = true;
         }
         catch( SPxException E )
         {
            MSG_ERROR( spxout << "Caught exception <" << E.what() << "> during unsimplification. Resolving without simplifier and scaler.\n" );
         }
         catch( ... )
         {
            MSG_ERROR( spxout << "Caught unknown exception during unsimplification. Resolving without simplifier and scaler.\n" );
            _status = SPxSolver::ERROR;
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
   void SoPlex::_storeSolutionReal()
   {
      assert(status() != SPxSolver::OPTIMAL || _isRealLPLoaded);

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
            && _solver.shift() < 10.0 * realParam(SoPlex::EPSILON_ZERO))) && _isRealLPLoaded;
      if( _solReal._hasPrimal )
      {
         _solReal._primal.reDim(_solver.nCols());
         _solReal._slacks.reDim(_solver.nRows());
         _solver.getPrimal(_solReal._primal);
         _solver.getSlacks(_solReal._slacks);
         _solReal._primalObjVal = _solver.objValue();
      }

      _solReal._hasPrimalRay = (status() == SPxSolver::UNBOUNDED && _isRealLPLoaded);
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
            && _solver.shift() < 10.0 * realParam(SoPlex::EPSILON_ZERO))) && _isRealLPLoaded;
      if( _solReal._hasDual )
      {
         _solReal._dual.reDim(_solver.nRows());
         _solReal._redCost.reDim(_solver.nCols());
         _solver.getDual(_solReal._dual);
         _solver.getRedCost(_solReal._redCost);
         _solReal._dualObjVal = ( _solReal._hasPrimal ? _solReal._primalObjVal : _solver.objValue() );
      }

      _solReal._hasDualFarkas = (status() == SPxSolver::INFEASIBLE && _isRealLPLoaded);
      if( _solReal._hasDualFarkas )
      {
         _solReal._dualFarkas.reDim(_solver.nRows());
         _solver.getDualfarkas(_solReal._dualFarkas);
      }

      _hasSolReal = true;
   }
} // namespace soplex
#endif
